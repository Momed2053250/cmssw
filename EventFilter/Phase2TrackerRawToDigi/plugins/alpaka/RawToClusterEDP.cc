// CMSSW includes
#include "DataFormats/FEDRawData/interface/StripPixelHostCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropHostCollection.h"
#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropDeviceCollection.h"

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/alpaka/ClusterPropSoACollection.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "CondFormats/SiPhase2TrackerObjects/interface/TrackerDetToDTCELinkCablingMap.h"
#include "CondFormats/SiPhase2TrackerObjects/interface/DTCELinkId.h"
#include "CondFormats/DataRecord/interface/TrackerDetToDTCELinkCablingMapRcd.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include <unordered_map>

#include "EventFilter/Phase2TrackerRawToDigi/interface/TrackerHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/ChannelsOffset.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerSpecifications.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2DAQFormatSpecification.h"
//#include "EventFilter/Phase2TrackerRawToDigi/plugins/alpaka/RawToCluster.dev.cc"
#include "EventFilter/Phase2TrackerRawToDigi/plugins/alpaka/RawToClusterAlgo.h"
#include <iomanip> // for std::setw
#include <future>
#include "FWCore/Framework/interface/ESWatcher.h"

//From the CPU based code 
using namespace Phase2TrackerSpecifications;
using namespace Phase2DAQFormatSpecification;
using namespace Phase2RawToCluster;


namespace ALPAKA_ACCELERATOR_NAMESPACE {

	using namespace cms::alpakatools;

	class Phase2RawToClusterProducer : public stream::EDProducer<> {
		public:
			explicit Phase2RawToClusterProducer(const edm::ParameterSet&);
			static void fillDescriptions(edm::ConfigurationDescriptions&);
			// enumaration declaration for the module types  
			enum WhichModule:int {undef, TwoS ,PS };
		private:
			void produce(device::Event&, device::EventSetup const&) override;
			void beginRun(edm::Run const&, edm::EventSetup const&) override;

			// Tokens
			const edm::EDGetTokenT<FEDRawDataCollection> fedRawDataToken_;
			const edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> cablingMapToken_;
			const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
			const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
			//New output token
			device::EDPutToken<Phase2RawToCluster::ClusterPropSoACollection> outputToken_;

			// cached ES pointers
			const TrackerDetToDTCELinkCablingMap* cablingMap_ = nullptr;
			const TrackerGeometry* trackerGeometry_ = nullptr;
			const TrackerTopology* trackerTopology_ = nullptr;
			std::map<int, std::pair<int,int>> stackMap_;
			// make the host buffers 
			// DTC*slink *channel -> DetIdx
			// Make new queue to access a queue before the produce method 
			Queue myqueue;
			//Get the available device for memory allocation 
			Device devAcc = alpaka::getDevByIdx(Platform{}, 0u);

			// TODO:: check if the buffer can be changed to SoAs 
			cms::alpakatools::host_buffer<int[]> detIdxModuleTypeMap_;
			cms::alpakatools::device_buffer<Device, int[]> detIdxModuleTypeDevice_;

	};

	Phase2RawToClusterProducer::Phase2RawToClusterProducer(const edm::ParameterSet& iConfig)
		: stream::EDProducer<>(iConfig),
		fedRawDataToken_(consumes<FEDRawDataCollection>(iConfig.getParameter<edm::InputTag>("fedRawDataCollection"))),
		cablingMapToken_(esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd, edm::Transition::BeginRun>()),
		trackerGeometryToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
		trackerTopologyToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
		outputToken_{ produces() },
		myqueue(devAcc),
		// Make the Global host and device buffer for each Module Type witht he given size = dtcId * number of slinks * number of channels 
		detIdxModuleTypeMap_{cms::alpakatools::make_host_buffer<int[], Platform>((MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC * CICs_PER_SLINK )},	
		detIdxModuleTypeDevice_{cms::alpakatools::make_device_buffer<int[]>(myqueue,(MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC * CICs_PER_SLINK )}
	{
	}

	void Phase2RawToClusterProducer::beginRun(
			edm::Run const&, edm::EventSetup const& iSetup) {
		cablingMap_     = &iSetup.getData(cablingMapToken_);
		trackerGeometry_ = &iSetup.getData(trackerGeometryToken_);
		trackerTopology_ = &iSetup.getData(trackerTopologyToken_);
		// Build the stack map 
		/*
		stackMap_.clear();
		for (auto const& detUnit : trackerGeometry_->detUnits()) {
			uint32_t rawId = detUnit->geographicalId().rawId();
			DetId detId(rawId);
			if (detId.det() != DetId::Detector::Tracker) continue;
			int stackIdx = trackerTopology_->stack(detId);
			if (trackerTopology_->isLower(detId))
				stackMap_[stackIdx].first = rawId;
			if (trackerTopology_->isUpper(detId))
				stackMap_[stackIdx].second = rawId;
		}*/
		stackMap_.clear();
for (auto iu = trackerGeometry_->detUnits().begin(); iu != trackerGeometry_->detUnits().end(); ++iu) {
    unsigned int detId_raw = (*iu)->geographicalId().rawId();
    DetId detId = DetId(detId_raw);
    if (detId.det() == DetId::Detector::Tracker) {
        int stackIdx = trackerTopology_->stack(detId);
        if (trackerTopology_->isLower(detId))
            stackMap_[stackIdx].first = detId;
        if (trackerTopology_->isUpper(detId))
            stackMap_[stackIdx].second = detId;
    }
}

		// Step 1) read for the module type and store the information in a global buffer: we do this in the begin run to just do it once not per event and we store this information in a buffer that can be used later 

		// Read one entire DTC (#dtcID), as per the producer logic

		for (int dtcID = MIN_DTC_ID; dtcID < MAX_DTC_ID + 1; dtcID++){

			// read the 4 slinks

			for (unsigned int iSlink = 0; iSlink < SLINKS_PER_DTC; iSlink++)
			{
				// as defined in the DAQProducer code
				// unsigned totID = iSlink + SLINKS_PER_DTC * (dtcID - 1) + CMSSW_TRACKER_ID ;


				// now read the payload (channel header + clusters)
				// all channel headers should be there, even if 0 clusters are found
				// the loop is not on the actual channel number, as in the ClusterToRaw conversion each channel is split by CIC0_CIC1
				// in order to get all the clusters from the same lpGBT and fill them once at the end
				// Loop over the Channels 
				for (unsigned int iChannel = 0; iChannel < CICs_PER_SLINK; iChannel++)
				{
					// retrieve the module type:
					// first we need to construct the DTCElinkId object ## dtc_id, gbtlink_id, elink_id
					// to get the gbt_id we should reverse what is done in the packer function,
					// where clusters from channel X are split into 2*i and 2*i+1 based on being from CIC0 or CIC1

					unsigned int gbt_id = iSlink * MODULES_PER_SLINK + std::div(iChannel, 2).quot;
					DTCELinkId thisDTCElinkId(dtcID, gbt_id, 0);

					int thisDetId = -1;
					bool is2SModule = false;

					// Define a detctor index as following : DetIdx = Channel Id + number of channels * sLink Id + number of channels * number of sLinks * DTCId 
					// To offset the DTCId to zero we do dtcID - MIN_DTC_ID
					auto DetIdx = iChannel + (CICs_PER_SLINK * iSlink) + (CICs_PER_SLINK * SLINKS_PER_DTC * (dtcID - MIN_DTC_ID));

					// then pass it to the map to get the detid
	/*				if (cablingMap_->knowsDTCELinkId(thisDTCElinkId))
					{
						auto possibleDetIds = cablingMap_->dtcELinkIdToDetId(thisDTCElinkId); // this returns a pair, detid will be an uint32_t (not a DetId)
						thisDetId = possibleDetIds->second;
						// Adding a printout for an intermediate check 
						std::cout << "slink is :" << iSlink << "\n" << "DtcID is:" << unsigned(dtcID) << "\n" << " detId is:" << thisDetId << "\n";
						//
						LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID) << "\tiGBT: " << unsigned(gbt_id)
							<< "\tielink: " << unsigned(0) << "\t -> detId:" << thisDetId;
						// check if it is 2S or PS
						is2SModule = trackerGeometry_->getDetectorType( stackMap_[thisDetId].first) == TrackerGeometry::ModuleType::Ph2SS;
						// If the slinks iDTC and iGBT are connected then check which module it is 
						detIdxModuleTypeMap_[DetIdx] = is2SModule ? WhichModule::TwoS : WhichModule::PS;		
					}
					else {
						LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID) << "\tiGBT: " << unsigned(gbt_id)
							<< " -> not connected? " ;
						// If the slinks iDTC and iGBT are not connected then return an undefined module 
						detIdxModuleTypeMap_[DetIdx] = WhichModule::undef;
						continue;
					}*/
					if (cablingMap_->knowsDTCELinkId(thisDTCElinkId)) {
    auto possibleDetIds = cablingMap_->dtcELinkIdToDetId(thisDTCElinkId); // returns a pair
    thisDetId = possibleDetIds->second;

    std::cout << "slink is :" << iSlink << "\n"
              << "DtcID is:" << unsigned(dtcID) << "\n"
              << "detId is:" << thisDetId << "\n";

    LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID)
                                     << "\tiGBT: " << unsigned(gbt_id)
                                     << "\tielink: " << unsigned(0)
                                     << "\t -> detId:" << thisDetId;

    // Get stack index and safely access stackMap_
    int stackIdx = trackerTopology_->stack(DetId(thisDetId));
    auto it = stackMap_.find(stackIdx);
    if (it != stackMap_.end()) {
        is2SModule = trackerGeometry_->getDetectorType(it->second.first) == TrackerGeometry::ModuleType::Ph2SS;
        detIdxModuleTypeMap_[DetIdx] = is2SModule ? WhichModule::TwoS : WhichModule::PS;
    } else {
        edm::LogWarning("RawToClusterProducer") << "Warning: stackIdx = " << stackIdx
                                                << " (from detId " << thisDetId << ") not found in stackMap_!";
        detIdxModuleTypeMap_[DetIdx] = WhichModule::undef;
    }
} else {
    LogTrace("RawToClusterProducer") << "slink: " << iSlink
                                     << "\tiDTC: " << unsigned(dtcID)
                                     << "\tiGBT: " << unsigned(gbt_id)
                                     << " -> not connected?";
    detIdxModuleTypeMap_[DetIdx] = WhichModule::undef;
    continue;
}

				} // channel 
			} // slink 
		} // det id 

		// Copy the information to the memory storages we created from the host read veiw to the device buffer store 
		alpaka::memcpy(
				myqueue,				
				detIdxModuleTypeDevice_,        // device destination pointer
				detIdxModuleTypeMap_, // host source pointer
				static_cast< unsigned int>(((MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC * CICs_PER_SLINK ) * sizeof(int))  // total bytes to copy
			      );
	} // Begin run 

	// Produce 
	void Phase2RawToClusterProducer::produce(
			device::Event& iEvent, device::EventSetup const&) {

		auto queue = iEvent.queue();


		// -------The commented section might be used to change from buffer approach to SoAs ------------- For now kept -------	  
		/*  
		// Create and fill host digi collection 
		// Construct new containers to hold data in host 
		// auto hostDigis = HGCalDigiHostCollection(newSize, queue);
		//	int size = some token.view().metadata().size(); ///will decide later 
		//auto hostStripPixel = StripPixelHostCollection(size, queue);
		for(int i=0; i<size;i++){
		//    hostStripPixel.view()[i].stripClusterWords() = some token.view()[i].stripClusterWords();
		//    hostStripPixel.view()[i].pixelClusterWords() = some token.view()[i].pixelClusterWords();

		}
		*/

		// 1) Build the flat rawword array:
		//FED Raw Collection as rowColl holds data fragments from each SLINK for every DTC
		auto const& rawColl = iEvent.get(fedRawDataToken_);
		// assert that the size is equal to number of slinks * dtcId 
		// TODO: Check if this assert is correct or there also has to be the * CICs_PER_SLINK also 
		assert(rawColl.size() == SLINKS_PER_DTC * (MAX_DTC_ID - MIN_DTC_ID));
		//store all rawColl in three vecotrs where these will store its data size and an offset
		//Data: raw byte
		std::vector<unsigned char> linearData;
		//size:size of each fragment, and offset: the starting index for each fragment in linearData  
		// TODO change from size_t to uint_32 : OPTIMIZATION  
		std::vector<size_t> size(rawColl.size());
		std::vector<size_t> offset(rawColl.size());
		// Fill size with the length of each FEDData fragment
		for (auto j = 0u; j < rawColl.size(); ++j) {
			auto& data = rawColl.FEDData(j);
			size.push_back(data.size());
		}
		// Compute offsets via exclusive scan:
		// offset[i] = sum of size[0] through size[i-1]
		std::exclusive_scan(size.begin(), size.end(), offset.begin(), 0);
		// Reserve total capacity for linearData: last offset + last fragment size
		linearData.reserve(offset[offset.size()-1] + size[size.size() - 1]);
		// Make a raw pointer to the beggining of the LinearData
		unsigned char* start = linearData.data();
		// Copy each fragment into linearData at its computed offset
		for (auto j = 0u; j < rawColl.size(); ++j) {
			auto& data = rawColl.FEDData(j);
			std::memcpy(start + offset[j], data.data(), data.size() );
		}
		// Make memory allocations to veiw these data from the CPU and copy them inot a buffer in GPU 
		auto linearData_HostView = cms::alpakatools::make_host_view<unsigned char>(linearData.data(), static_cast<long unsigned int>(linearData.size()));
		auto linearData_DevBuffer = cms::alpakatools::make_device_buffer<unsigned char[]>(queue, static_cast<long unsigned int>(linearData.size()));
		alpaka::memcpy(
				queue,
				linearData_DevBuffer,        // device destination pointer
				linearData_HostView, // host source pointer
				static_cast< unsigned int>( linearData.size() )  // total bytes to copy
			      ); 
		auto size_HostView = cms::alpakatools::make_host_view<size_t>(size.data(), static_cast<long unsigned int>(size.size()));
		auto size_DevBuffer = cms::alpakatools::make_device_buffer<size_t[]>(queue, static_cast<long unsigned int>(size.size()));
		alpaka::memcpy(
				queue,
				size_DevBuffer,        // device destination pointer
				size_HostView , // host source pointer
				static_cast< unsigned int>( size.size() )  // total bytes to copy
			      );
		auto offset_HostView = cms::alpakatools::make_host_view<size_t>(offset.data(), static_cast<long unsigned int>(offset.size()));
		auto offset_DevBuffer = cms::alpakatools::make_device_buffer<size_t[]>(queue, static_cast<long unsigned int>(offset.size()));
		alpaka::memcpy(
				queue,
				offset_DevBuffer,        // device destination pointer
				offset_HostView , // host source pointer
				static_cast< unsigned int>( offset.size() )  // total bytes to copy
			      );
		// wait for the copy to finish before launching kernels
		alpaka::wait(queue);

		// Launch the kernals 
		launchUnpacker(queue, linearData_DevBuffer, size_DevBuffer, offset_DevBuffer ,detIdxModuleTypeDevice_);
		// ----------------------------Commented section kept for later stages -------------------------//
		/*
		// 2) Fill host input SOA
		auto hostStripPixel = StripPixelHostCollection(totalWords, queue);
		for (size_t i = 0; i < totalWords; ++i) {
		hostStripPixel.view()[i].stripClustersWords() = rawWords[i];
		}

		// 3) Allocate device input buffer and memcpy the host SOA into it
		Phase2RawToCluster::StripPixelDeviceCollection devStripPixel(totalWords, queue);
		//aloocating a buffer 
		//	const uint32_t wordCounter = 0;
		//	      auto buffer = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, wordCounter);
		alpaka::memcpy(
		queue,
		devStripPixel.buffer(),        // device destination pointer
		hostStripPixel.buffer() //, // host source pointer
		//  totalWords * sizeof(Phase2RawToCluster::StripPixelSoA::SoALayout)  // total bytes to copy
		);
		// once you know devStripPixel is filled, allocate the device output SOA
		Phase2RawToCluster::ClusterPropDeviceCollection devClusterProp(totalWords, queue);
		// wait for the copy to finish before launching your kernel
		//cms::alpakatools::CopyToDevice(queue, devStripPixel, hostStripPixel);
		alpaka::wait(queue);

		// 5) Allocate a host-side SOA for the output and memcpy the result back
		//	  ClusterPropHostCollection 
		Phase2RawToCluster::ClusterPropSoACollection  hostClusterProp(totalWords, queue);
		alpaka::memcpy(
		queue,
		hostClusterProp.buffer(),      // host destination pointer
		devClusterProp.const_buffer() //, // device source pointer
		//   totalWords * sizeof(Phase2RawToCluster::ClusterPropSoA::SoALayout)
		);
		// wait for the device→host copy to complete
		alpaka::wait(queue);


		// -------------------------------------------------- iEvent.emplace -----------------------------------//
		// finally, put it into the event
		iEvent.emplace(outputToken_, std::move(hostClusterProp));
		// call the unoacker kernal 


*/
	} //produce 


	void Phase2RawToClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
		edm::ParameterSetDescription desc;
		desc.add<edm::InputTag>("fedRawDataCollection");
		descriptions.addWithDefaultLabel(desc);
	}



}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// define this as a plug-in
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(Phase2RawToClusterProducer);
