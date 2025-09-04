// EDProducer 
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
#include "EventFilter/Phase2TrackerRawToDigi/plugins/alpaka/RawToClusterAlgo.h"
#include <iomanip> // for std::setw
#include <future>
#include "FWCore/Framework/interface/ESWatcher.h"

//From the CPU based code 
using namespace Phase2TrackerSpecifications;
using namespace Phase2DAQFormatSpecification;
using namespace Phase2RawToCluster;

// debug flag
//#define Debug_CPU


namespace ALPAKA_ACCELERATOR_NAMESPACE {

	using namespace cms::alpakatools;

	class Phase2RawToClusterProducer : public stream::EDProducer<> {
		public:
			explicit Phase2RawToClusterProducer(const edm::ParameterSet&);
			static void fillDescriptions(edm::ConfigurationDescriptions&);
			void beginRun(edm::Run const&, edm::EventSetup const&) override;
			// enumaration declaration for the module types  
			enum WhichModule:int {undef, TwoS ,PS };
		private:
			void produce(device::Event&, device::EventSetup const&) override;

			// Tokens for aquiring the RAW data 
			const edm::EDGetTokenT<FEDRawDataCollection> fedRawDataToken_;
			const edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> cablingMapToken_;
			const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
			const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
			//New output token
			device::EDPutToken<Phase2RawToCluster::ClusterPropSoACollection> outputToken_;
			// For the converted one
			//edm::EDPutTokenT<Phase2TrackerCluster1DCollectionNew> clusterCollectionToken_;
			// cached ES pointers
			const TrackerDetToDTCELinkCablingMap* cablingMap_ = nullptr;
			const TrackerGeometry* trackerGeometry_ = nullptr;
			const TrackerTopology* trackerTopology_ = nullptr;
			std::map<int, std::pair<int,int>> stackMap_;
			// make the host buffers 
			// DTC*slink *channel -> detIdx
			//Get the available device for memory allocation 
			Device devAcc = alpaka::getDevByIdx(Platform{}, 0u);
			// Make new queue to access a queue before the produce method 
			Queue myqueue;
			// TODO:: check if the buffer can be changed to SoAs 
			cms::alpakatools::host_buffer<int[]> detIdxModuleTypeMap_;
			cms::alpakatools::device_buffer<Device, int[]> detIdxModuleTypeDevice_;
			// TODO: check if the det id is needed at the end or this can be deleted later (depends on the output data grouping format)
			// To use the detid from the stack map store the stack map
			// a. per process  
			cms::alpakatools::host_buffer<std::pair<int, int>[]> stackMapHost_;
			cms::alpakatools::device_buffer<Device, std::pair<int, int>[]> stackMapDevice_;
			// b. per event 
			cms::alpakatools::host_buffer<int[]> detIdMapHost_;
			cms::alpakatools::device_buffer<Device, int[]> detIdMapDevice_;

	};

	Phase2RawToClusterProducer::Phase2RawToClusterProducer(const edm::ParameterSet& iConfig)
		: stream::EDProducer<>(iConfig),
		fedRawDataToken_(consumes<FEDRawDataCollection>(iConfig.getParameter<edm::InputTag>("fedRawDataCollection"))),
		cablingMapToken_(esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd, edm::Transition::BeginRun>()),
		trackerGeometryToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
		trackerTopologyToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
		outputToken_{ produces() },
		//clusterCollectionToken_{produces<edmNew::DetSetVector<Phase2TrackerCluster1D>>()} // for the converted ones 
		myqueue(devAcc),
		// Make the Global host and device buffer for each Module Type with the given size = dtcId * number of slinks * number of channels -> detIdx 
		//// include +1 because DTC IDs run from MIN_DTC_ID through MAX_DTC_ID inclusive this solves the buffer overflow issue in run time
		detIdxModuleTypeMap_{cms::alpakatools::make_host_buffer<int[], Platform>((MAX_DTC_ID - MIN_DTC_ID +1) * SLINKS_PER_DTC * CICs_PER_SLINK )},	
		detIdxModuleTypeDevice_{cms::alpakatools::make_device_buffer<int[]>(myqueue,(MAX_DTC_ID - MIN_DTC_ID +1) * SLINKS_PER_DTC * CICs_PER_SLINK )},
		// initialize the stack map host and device buffers
		// In Phase2RawToClusterProducer constructor
		stackMapHost_{cms::alpakatools::make_host_buffer<std::pair<int, int>[], Platform>(stackMap_.size())},
		stackMapDevice_{cms::alpakatools::make_device_buffer<std::pair<int, int>[]>(myqueue, stackMap_.size())},
		//stackMapHost_{},
		//stackMapDevice_{},
		detIdMapHost_{cms::alpakatools::make_host_buffer<int[], Platform>((MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC * CICs_PER_SLINK)},
		detIdMapDevice_{cms::alpakatools::make_device_buffer<int[]>(myqueue, (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC * CICs_PER_SLINK)}
		{
	}

	void Phase2RawToClusterProducer::beginRun(
			edm::Run const&, edm::EventSetup const& iSetup) {
		cablingMap_     = &iSetup.getData(cablingMapToken_);
		trackerGeometry_ = &iSetup.getData(trackerGeometryToken_);
		trackerTopology_ = &iSetup.getData(trackerTopologyToken_);
		// Build the stack map 
		stackMap_.clear();
		for (auto iu = trackerGeometry_->detUnits().begin(); iu != trackerGeometry_->detUnits().end(); ++iu) {
			unsigned int detId_raw = (*iu)->geographicalId().rawId();
			DetId detId = DetId(detId_raw);
			if (detId.det() == DetId::Detector::Tracker) {
				// build map of upper and lower for each module
				if (trackerTopology_->isLower(detId) != 0) {
					stackMap_[trackerTopology_->stack(detId)].first = detId;
				}
				if (trackerTopology_->isUpper(detId) != 0) {
					stackMap_[trackerTopology_->stack(detId)].second = detId;
				}
			}
		} 

// Copy the stack map to the device buffer
const unsigned int sM = static_cast<unsigned int>(stackMap_.size());
// Now (re)allocate the host/device buffers with correct size
stackMapHost_ = cms::alpakatools::make_host_buffer<std::pair<int, int>[], Platform>(sM);
stackMapDevice_ = cms::alpakatools::make_device_buffer<std::pair<int, int>[]>(myqueue, sM);
// Fill host buffer
size_t stackIdx = 0;
for (const auto& entry : stackMap_) {
  stackMapHost_[stackIdx++] = entry.second;
}
//alpaka::memset(myqueue, stackMapDevice_, 0x00);
alpaka::memcpy(myqueue, stackMapDevice_, stackMapHost_, sM);
alpaka::wait(myqueue);

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
					if (cablingMap_->knowsDTCELinkId(thisDTCElinkId)) {
						auto possibleDetIds = cablingMap_->dtcELinkIdToDetId(thisDTCElinkId); // returns a pair
						thisDetId = possibleDetIds->second;
						// After `thisDetId = possibleDetIds->second; // TODO: check  
						detIdMapHost_[DetIdx] = thisDetId;

#ifdef Debug_CPU
						std::cout << "slink is :" << iSlink << "\n"
							<< "DtcID is:" << unsigned(dtcID) << "\n"
							<< "detId is:" << thisDetId << "\n";
#endif

						LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID)
							<< "\tiGBT: " << unsigned(gbt_id)
							<< "\tielink: " << unsigned(0)
							<< "\t -> detId:" << thisDetId;
						// check is 2S or PS
						is2SModule = trackerGeometry_->getDetectorType( stackMap_[thisDetId].first) == TrackerGeometry::ModuleType::Ph2SS;
						//std::cout << "is2SModule is:" << is2SModule << "\n"; // here it is correctly calculated 
						detIdxModuleTypeMap_[DetIdx] = is2SModule ? WhichModule::TwoS : WhichModule::PS;	
						
#ifdef Debug_CPU
						std::cout << "Mapped dtcID=" << dtcID << ", iSlink=" << iSlink << ", iChannel=" << iChannel 
							<< ", gbt_id=" << gbt_id << ", detId=" << thisDetId << ", stackIdx=" << stackIdx 
							<< //", moduleType=" << (is2SModule ? "TwoS" : "PS") <<", DetIdx=" << DetIdx << 
							"\n";
#endif

					}
					else {
						LogTrace("RawToClusterProducer")
							<< "slink: " << iSlink
							<< "\tiDTC: " << unsigned(dtcID)
							<< "\tiGBT: " << unsigned(gbt_id)
							<< " -> not connected?";
						detIdxModuleTypeMap_[DetIdx] = WhichModule::undef;
						continue;
					} 
				} // channel 
			} // slink 
		} // det id 
// copy detidmap to the device buffer
//constexpr unsigned int M = 4097;
const unsigned int M =
    (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC * CICs_PER_SLINK;
//alpaka::memset( myqueue, detIdMapDevice_, 0x00 );
alpaka::memcpy(
	myqueue, 
	detIdMapDevice_, 
	detIdMapHost_, 
	M
);
alpaka::wait(myqueue);
		// Copy the information to the memory storages we created from the host read veiw to the device buffer store 
		// TODO: check if N is correct or needs to be the same as detidx 
		constexpr unsigned N = (MAX_DTC_ID - MIN_DTC_ID +1)*SLINKS_PER_DTC*CICs_PER_SLINK;
		//alpaka::memset( myqueue,detIdxModuleTypeDevice_, 0x00 );
		alpaka::memcpy(
				myqueue,
				detIdxModuleTypeDevice_,
				detIdxModuleTypeMap_,
				N    //   number of int elements
			      );
		// wait for the copy to finish 
		alpaka::wait(myqueue);
		
	} // Begin run 

	// Produce 
	// Produce 
void Phase2RawToClusterProducer::produce(
    device::Event& iEvent, device::EventSetup const&) {

  auto queue = iEvent.queue();
  // maximum total clusters 
  static constexpr size_t MaxTotalClusters = (N_CLUSTER_MASK + 1) * CICs_PER_SLINK * (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;
  // -------The commented section might be used to change from buffer approach to SoAs ------------- For now kept -------	  
  
  /*   Example 
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
  // CHANGED: Compute the number of possible slinks to match kernel loop range and CPU logic
  const size_t numSlinks = (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;
  // as defined in the DAQProducer code
  //unsigned totID = iSlink + SLINKS_PER_DTC * (dtcID - 1) + CMSSW_TRACKER_ID ;   // not used in the port  
  // assert that the size is equal to number of slinks * dtcId 
  // TODO: Check if this assert is correct or there also has to be the * CICs_PER_SLINK also 
  //	assert(rawColl.size() == SLINKS_PER_DTC * (MAX_DTC_ID - MIN_DTC_ID));
  //store all rawColl in three vecotrs where these will store its data size and an offset
  //Data: raw byte
  std::vector<unsigned char> linearData;
  //size:size of each fragment, and offset: the starting index for each fragment in linearData  
  // TODO change from size_t to uint_32 : OPTIMIZATION  
  // CHANGED: Resize size and offset to numSlinks (instead of rawColl.size()) to only handle relevant tracker slinks
  std::vector<size_t> size(numSlinks);          // potential zero hsunting feild this is intialized with all zeros 
  std::vector<size_t> offset(numSlinks);
  // CHANGED: Add totIDs vector to store computed FED IDs for each slinkIdx
  std::vector<unsigned int> totIDs(numSlinks);
  // CHANGED: Loop over dtcID and iSlink like original CPU code to compute totID and fetch only relevant FED data
  size_t slinkIdx = 0;
  for (int dtcID = MIN_DTC_ID; dtcID < MAX_DTC_ID + 1; dtcID++) {
    for (unsigned int iSlink = 0; iSlink < SLINKS_PER_DTC; iSlink++) {
      unsigned totID = iSlink + SLINKS_PER_DTC * (dtcID - 1) + CMSSW_TRACKER_ID;
      totIDs[slinkIdx] = totID;
      const FEDRawData& fedData = rawColl.FEDData(totID);
      size[slinkIdx] = fedData.size();
      #ifdef Debug_CPU
      // print the all the elements of the fed raw data to check if zeros are there 
      // CHANGED: Print with totID instead of j for consistency with CPU debug
      std::cout << "FEDRawDataCollection[" << totID << "] size: " << size[slinkIdx] << "\n";
      #endif
      slinkIdx++;
    }
  }
  //TODO:: Check this section if needed to move back to reserve from resize :: Perhaps not 
  // Compute offsets via exclusive scan:
  // offset[i] = sum of size[0] through size[i-1]
  std::exclusive_scan(size.begin(), size.end(), offset.begin(), 0);
  // Reserve total capacity for linearData: last offset + last fragment size
  //linearData.reserve(offset[offset.size()-1] + size[size.size() - 1]);
  //linearData.reserve(offset.back() + size.back());
  //changing reserve to resize in attempt to solve the illegal memory access runtime erro 
  size_t totalBytes = offset.back() + size.back();
  linearData.resize(totalBytes);
  // commenting/adding to solve the init problem
  // ====== NEW: EXPLICITLY INITIALIZE HOST MEMORY ======
  std::fill(linearData.begin(), linearData.end(), 0);  // ZERO-INITIALIZE ENTIRE BUFFER 
  // Now linearData.data() points to a buffer of length = totalBytes
  // Make a raw pointer to the beggining of the LinearData
  unsigned char* start = linearData.data();
  // Copy each fragment into linearData at its computed offset
  // CHANGED: Loop over numSlinks (0 to numSlinks-1), use totIDs[idx] to fetch data, skip if size[idx]==0
  for (size_t idx = 0; idx < numSlinks; ++idx) {
    if (size[idx] == 0) continue; // Skip empty fragments
    const FEDRawData& data = rawColl.FEDData(totIDs[idx]);
    // ====== NEW: ADDED BOUNDS CHECK FOR SAFETY ======
    if (offset[idx] + size[idx] > totalBytes) {
      throw std::runtime_error("BUFFER OVERFLOW DETECTED IN RAW DATA COPYING");
    }
    //if (data.size() == 0) continue; // Skip empty fragments make no diffrrence in the zero peak 
    std::memcpy(start + offset[idx], data.data(), size[idx]);
  }
  // Make memory allocations to veiw these data from the CPU and copy them inot a buffer in GPU 
  auto linearData_HostView = cms::alpakatools::make_host_view<unsigned char>(linearData.data(), static_cast<long unsigned int>(linearData.size()));
  auto linearData_DevBuffer = cms::alpakatools::make_device_buffer<unsigned char[]>(queue, static_cast<long unsigned int>(linearData.size()));
  // ====== NEW: INITIALIZE DEVICE BUFFER BEFORE COPY ======
  // CHANGED: Uncommented alpaka::memset for linearData_DevBuffer to ensure no garbage data
  alpaka::memset(queue, linearData_DevBuffer, 0x00);  // ZERO DEVICE MEMORY FIRST
  alpaka::memcpy(
      queue,
      linearData_DevBuffer,        // device destination pointer
      linearData_HostView, // host source pointer
      static_cast< unsigned int>( linearData.size() )  // total bytes to copy
        ); 
      alpaka::wait(queue);
  auto size_HostView = cms::alpakatools::make_host_view<size_t>(size.data(), static_cast<long unsigned int>(size.size()));
  auto size_DevBuffer = cms::alpakatools::make_device_buffer<size_t[]>(queue, static_cast<long unsigned int>(size.size()));
  // ====== NEW: PROPER BYTE COUNT FOR SIZE_T BUFFER ======
  // CHANGED: Uncommented alpaka::memset for size_DevBuffer to ensure no garbage data
  alpaka::memset(queue, size_DevBuffer, 0x00);
  alpaka::memcpy(
      queue,
      size_DevBuffer,        // device destination pointer
      size_HostView , // host source pointer
      static_cast< unsigned int>( size.size() )  // total bytes to copy
        );
      alpaka::wait(queue);
  auto offset_HostView = cms::alpakatools::make_host_view<size_t>(offset.data(), static_cast<long unsigned int>(offset.size()));
  auto offset_DevBuffer = cms::alpakatools::make_device_buffer<size_t[]>(queue, static_cast<long unsigned int>(offset.size()));
  // ====== NEW: SAME FIXES FOR OFFSET BUFFER ======
  // CHANGED: Uncommented alpaka::memset for offset_DevBuffer to ensure no garbage data
  alpaka::memset(queue, offset_DevBuffer, 0x00);
  alpaka::memcpy(
    queue,
    offset_DevBuffer,        // device destination pointer
    offset_HostView ,        // host source pointer
    static_cast<unsigned int>(offset.size())  // total bytes to copy
  );
  alpaka::wait(queue);

  // Check this part 
  // Allocate output SoA and global counter
  auto devClusterProp = Phase2RawToCluster::ClusterPropDeviceCollection(MaxTotalClusters, queue);
  auto&& devClusterPropBuffer = devClusterProp.buffer();  // Use forwarding reference
  // CHANGED: Uncommented alpaka::memset for devClusterPropBuffer to ensure no garbage data
  alpaka::memset(queue, devClusterPropBuffer, 0x00);  // Zero-initialize the device buffer
  auto globalCounter = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, 1u);
  // ✅ initialize the counter so we don’t read garbage later
  alpaka::memset(queue, globalCounter, 0u);
  alpaka::wait(queue);

  // wait for the copy to finish before launching kernels
  //alpaka::wait(queue);

  // Launch the kernals 
  //launchUnpacker(queue, linearData_DevBuffer, size_DevBuffer, offset_DevBuffer, detIdxModuleTypeDevice_,
  //              devClusterProp.view(), globalCounter.data()); 
  //New: launch the unpacker with the stackmapinfor for retainig the detId
  // Modify the launchUnpacker call
  //launchUnpacker(queue, linearData_DevBuffer, size_DevBuffer, offset_DevBuffer, detIdxModuleTypeDevice_, detIdMapDevice_, stackMapDevice_, devClusterProp.view(), globalCounter.data());        
  //auto hostClusterProp = ClusterPropHostCollection(MaxTotalClusters, queue);
  uint32_t stackMapSize = stackMap_.size();
  launchUnpacker(
  queue,
  linearData_DevBuffer,
  size_DevBuffer,
  offset_DevBuffer,
  detIdxModuleTypeDevice_,
  detIdMapDevice_,
  stackMapDevice_,
  stackMapSize,  // Pass the size
  devClusterProp.view(),
  globalCounter.data()
  );
  //alpaka::memcpy(queue, hostClusterProp.buffer(), devClusterProp.const_buffer());
  //alpaka::wait(queue);

  // Copy output back to host as ClusterPropSoACollection
  Phase2RawToCluster::ClusterPropSoACollection hostClusterPropSoA(MaxTotalClusters, queue);
  auto hostBuf = hostClusterPropSoA.buffer();
  // CHANGED: Uncommented alpaka::memset for hostBuf to ensure no garbage data (though host-side, it's safer)
  alpaka::memset(queue, hostBuf, 0x00);
  alpaka::memcpy(queue, hostBuf, devClusterProp.const_buffer());
  alpaka::wait(queue);


// Print the size column of the SoA
// needs check returs a segmentation violation in runtime  
//#ifdef Debug_CPU
/*
// 1. Allocate host buffer to read back the counter
cms::alpakatools::host_buffer<uint32_t[]> hostCounter =
    cms::alpakatools::make_host_buffer<uint32_t[], Platform>(1u);

// 2. Copy from device to host
alpaka::memcpy(queue, hostCounter, globalCounter, 1u);
alpaka::wait(queue);  // make sure copy completes

// 3. Read number of clusters
size_t numElements = hostCounter[0];
auto viewSoA = hostClusterPropSoA.view();

std::cout << "Number of clusters: " << numElements << std::endl;

// 4. Clamp to allocated capacity
const size_t cap = static_cast<size_t>(viewSoA.metadata().size());
if (numElements > cap) {
    std::cerr << "[WARN] globalCounter (" << numElements
              << ") exceeds SoA capacity (" << cap << "). Clamping.\n";
}
const size_t safeN = std::min(numElements, cap);

// 5. Print first few entries (size only)
const size_t toPrint = std::min<size_t>(safeN, 64);
for (size_t i = 0; i < toPrint; ++i) {
    std::cout << "Element " << i << ": size = " << viewSoA[i].size() << std::endl;
}
if (safeN > toPrint) {
    std::cout << "... (" << (safeN - toPrint) << " more not shown)\n";
}
*/
//#endif

// converter SoA->DetSetVector  

//  Convert SOA to DetSetVector for compatibility with original output
/*
auto outputClusterCollection = std::make_unique<Phase2TrackerCluster1DCollectionNew>();
auto view = hostClusterPropSoA.view();
std::map<uint32_t, edmNew::DetSet<Phase2TrackerCluster1D>> detSets;
size_t sizeC = view.metadata().size();
for (size_t i = 0; i < sizeC; ++i) {
	uint32_t detId = view[i].detId();
	Phase2TrackerDigi firstDigi(view[i].strip(), view[i].row());
	Phase2TrackerCluster1D cluster(firstDigi, view[i].size(), view[i].threshold());
	detSets[detId].data().push_back(cluster);
}

for (auto& [detId, detSet] : detSets) {
	edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller ff(*outputClusterCollection, detId);
	for (const auto& cluster : detSet.data) {
		ff.push_back(cluster);
	}
} 

auto outputClusterCollection = std::make_unique<Phase2TrackerCluster1DCollectionNew>();
auto view = hostClusterPropSoA.view();

size_t sizeC = view.metadata().size();
std::unordered_map<uint32_t, edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller> fillers;

for (size_t i = 0; i < sizeC; ++i) {
    uint32_t detId = view[i].detId();
    Phase2TrackerDigi firstDigi(view[i].strip(), view[i].row());
    Phase2TrackerCluster1D cluster(firstDigi, view[i].size(), view[i].threshold());

    // Only create a FastFiller once per detId
    auto it = fillers.find(detId);
    if (it == fillers.end()) {
        auto [newIt, success] = fillers.emplace(
            detId,
            edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller(*outputClusterCollection, detId));
        it = newIt;
    }

    it->second.push_back(cluster);
}

// END ADDED

// ADDED: Put the DetSetVector into the event
//iEvent.put(std::move(outputClusterCollection));
//iEvent.emplace(clusterCollectionToken_, std::move(outputClusterCollection));
*/
// Emplace the SoA collection into the event
iEvent.emplace(outputToken_, std::move(hostClusterPropSoA));

		// ----------------------------Commented section kept for later stages -------------------------//
		/*
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
		*/
		//finally, put it into the event
		//iEvent.emplace(outputToken_, Phase2RawToCluster::ClusterPropSoACollection(std::move(hostClusterProp)));


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
