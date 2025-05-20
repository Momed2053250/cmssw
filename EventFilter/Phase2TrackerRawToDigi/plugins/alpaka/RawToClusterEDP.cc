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
// just a peaceful timer 
template<class T> double duration(T t0,T t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

inline std::chrono::time_point<std::chrono::steady_clock> now()
{
  return std::chrono::steady_clock::now();
}

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  class Phase2RawToClusterProducer : public stream::EDProducer<> {
  public:
    explicit Phase2RawToClusterProducer(const edm::ParameterSet&);
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    void produce(device::Event&, device::EventSetup const&) override;
    void beginRun(edm::Run const&, edm::EventSetup const&) override;
   //  edm::ESWatcher<HGCalCondSerializableModuleInfoRcd> calibWatcher_;
  //  edm::ESWatcher<HGCalCondSerializableConfigRcd> configWatcher_;
  // Token
    const edm::EDGetTokenT<FEDRawDataCollection> fedRawDataToken_;
    const edm::ESGetToken<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd> cablingMapToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
   // edm::EDGetTokenT<Phase2TrackerCluster1DCollectionNew> OutputClusterCollectionToken_;
    // our new output token
    //cms::alpaka
	   device::EDPutToken<Phase2RawToCluster::ClusterPropSoACollection> outputToken_;

       // cached ES pointers
    const TrackerDetToDTCELinkCablingMap* cablingMap_ = nullptr;
    const TrackerGeometry* trackerGeometry_ = nullptr;
    const TrackerTopology* trackerTopology_ = nullptr;
    std::map<int, std::pair<int,int>> stackMap_;
  };

    Phase2RawToClusterProducer::Phase2RawToClusterProducer(const edm::ParameterSet& iConfig)
    : stream::EDProducer<>(iConfig),
     fedRawDataToken_(consumes<FEDRawDataCollection>(iConfig.getParameter<edm::InputTag>("fedRawDataCollection"))),
      cablingMapToken_(esConsumes<TrackerDetToDTCELinkCablingMap, TrackerDetToDTCELinkCablingMapRcd, edm::Transition::BeginRun>()),
      trackerGeometryToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
      trackerTopologyToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
      
	// outputToken_(produces())
	//outputToken_{produces<Phase2RawToCluster::ClusterPropHostCollection>()}
	outputToken_{ produces() //
         }
       //		.template deviceProduces<Phase2RawToCluster::ClusterPropDeviceCollection, Phase2RawToCluster::ClusterPropHostCollection>() }
// outputToken_{ produces<Phase2RawToCluster::ClusterPropHostCollection>() }
	{}

//Phase2RawToClusterProducer::~Phase2RawToClusterProducer() 
//{
//}

  void Phase2RawToClusterProducer::beginRun(
      edm::Run const&, edm::EventSetup const& iSetup) {
    cablingMap_     = &iSetup.getData(cablingMapToken_);
    trackerGeometry_ = &iSetup.getData(trackerGeometryToken_);
    trackerTopology_ = &iSetup.getData(trackerTopologyToken_);

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
    }
  }

  void Phase2RawToClusterProducer::produce(
      device::Event& iEvent, device::EventSetup const&) {
    auto queue = iEvent.queue();
/*
    // Read digis
   // auto const& deviceCalibParamProvider = iSetup.getData(calibToken_);
  //  auto const& deviceConfigParamProvider = iSetup.getData(configToken_);
  //  auto const& hostDigisIn = iEvent.get(digisToken_);


    
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
/*
    	// 1) Load raw FED words into host vector
  	auto const& rawColl = iEvent.get(fedRawDataToken_);
  	size_t totalWords = 0;
  	for (auto const& fedIt : rawColl) {
    	totalWords += fedIt.second.size() / sizeof(uint32_t);
  	}
  	std::vector<uint32_t> rawWords;
  	rawWords.reserve(totalWords);
  	for (auto const& fedIt : rawColl) {
    	auto ptr32 = reinterpret_cast<const uint32_t*>(fedIt.second.data());
    	auto nWords = fedIt.second.size() / sizeof(uint32_t);
    	rawWords.insert(rawWords.end(), ptr32, ptr32 + nWords);
  	}
*/
	// 1) Build the flat raw‐word array from exactly those FEDs the cabling map knows:

	auto const& rawColl = iEvent.get(fedRawDataToken_);

	// (a) Ask the map for every DTCELinkId it knows
	auto dtcLinkIds = cablingMap_->getKnownDTCELinkIds();

	// (b) Turn each DTCELinkId into a unique FED ID
	std::vector<int> fedIds;
	fedIds.reserve(dtcLinkIds.size());
	for (auto const& linkId : dtcLinkIds) {
	  // extract the DTC number (1…MAX_DTC_ID) and the S-link index (0…SLINKS_PER_DTC-1)
	  unsigned dtcId     = linkId.dtc_id();
	  unsigned slinkIdx  = linkId.gbtlink_id();

	  int fedId = slinkIdx
	            + SLINKS_PER_DTC * (dtcId - 1)
	            + CMSSW_TRACKER_ID;
	  fedIds.push_back(fedId);
	}

	// sort+unique so we don’t double-count any
	std::sort(fedIds.begin(), fedIds.end());
	fedIds.erase(std::unique(fedIds.begin(), fedIds.end()), fedIds.end());

	// (c) First pass: count total 32-bit words
	size_t totalWords = 0;
	for (int fedId : fedIds) {
	  auto const& fed = rawColl.FEDData(fedId);
	  totalWords += fed.size() / sizeof(uint32_t);
	}

	// (d) Reserve and fill the flat word array
	std::vector<uint32_t> rawWords;
	rawWords.reserve(totalWords);

	for (int fedId : fedIds) {
	  auto const& fed = rawColl.FEDData(fedId);
	  auto ptr32  = reinterpret_cast<const uint32_t*>(fed.data());
	  auto nWords = fed.size() / sizeof(uint32_t);
	  rawWords.insert(rawWords.end(), ptr32, ptr32 + nWords);
	}


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

	// 4) Launch 2S unpacker
   	launchS2UnpackerKernel(queue, devStripPixel, devClusterProp, totalWords);
	alpaka::wait(queue);
/*	const auto threadsPerBlock = 512;
StripPixelHos	auto blocks = cms::alpakatools::divide_up_by(static_cast<uint32_t>(totalWords), threadsPerBlock);
	auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);
  	alpaka::exec<Acc1D>(queue, workDiv, S2UnpackerKernel{}, devStripPixel.view(), devClusterProp.view());
  	alpaka::wait(queue);
*/

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
	  // finally, put it into the event
	  iEvent.emplace(outputToken_, std::move(hostClusterProp));

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
