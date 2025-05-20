// alpaka-related imports
#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/alpaka/ClusterPropDeviceCollection.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/SensorHybrid.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerSpecifications.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2DAQFormatSpecification.h"

// Calibration algorithms header
//#include "RecoLocalCalo/HGCalRecAlgos/plugins/alpaka/HGCalRecHitCalibrationAlgorithms.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  using namespace Phase2RawToCluster;
  using namespace Phase2TrackerSpecifications;
  using namespace Phase2DAQFormatSpecification;
// a class kernal example 
/*  class CalibrationKernel_digisToRecHits {
  public:

    template <typename TAcc>

    ALPAKA_FN_ACC void operator()(TAcc const& acc, StripPixelDeviceCollection::View strippixel, StripPixelOutputDeviceCollection::View output ) const {
      // dummy example
      for (auto index : elements_with_stride(acc, strippixel.metadata().size())) {
        output[index].stripClustersWords() = strippixel[index].stripClusterWords();
	output[index].pixelClustersWords() = strippixel[index].pixelClusterWords();
      }
    }
  };
  */

  // Kernel for 2S modules: unpack only stripClustersWords into (x,y,width)
	struct S2UnpackerKernel {
	  template <typename Acc, typename InView, typename OutView>
	  ALPAKA_FN_ACC void operator()(Acc const& acc,
	                                InView const& in,
	                                OutView out) const {
	    auto n = in.metadata().size();
	    for (auto idx : uniform_elements(acc, n)) {
	      uint32_t word = in[idx].stripClustersWords();
	      uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
	      uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_ONLY_BITS_2S)) & SCLUSTER_ADDRESS_MASK;
	      bool seed = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_2S)) & 1;
	      uint32_t w = word & WIDTH_MAX_VALUE;
	      if (w == 0) w = 8;
	      out[idx].x()        = STRIPS_PER_CBC * chip + addr;
	      out[idx].y()        = seed ? 1u : 0u;
	      out[idx].width()    = w;
	      out[idx].seedFlag() = seed;
	    }
	  }
	};	
        // 4th step  of the EDP.cc) Launch 2S unpacker
      void launchS2UnpackerKernel(
	    Queue& queue,
	    const Phase2RawToCluster::StripPixelDeviceCollection& devStripPixel,
	    Phase2RawToCluster::ClusterPropDeviceCollection& devClusterProp,
	    size_t totalWords) {
	  const uint32_t threadsPerBlock = 512;
	  const uint32_t blocks = cms::alpakatools::divide_up_by(totalWords, threadsPerBlock);
	  auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);
	  alpaka::exec<Acc1D>(queue, workDiv, S2UnpackerKernel{}, devStripPixel.view(), devClusterProp.view());
	}

  /* NEEED TO BE RECHECKED ONCE THE 2S MODEL WORKS THIS WAY

// Kernel for PS modules: unpack both strip **and** pixel words
  struct PSUnpackerKernel {
    template <typename Acc>
    ALPAKA_FN_ACC void operator()(Acc const & acc,
                                  StripPixelDeviceCollection::View in,
                                  ClusterPropDeviceCollection::View out) const {
      auto n = in.metadata().size();
      // for now I assume the first nStrip entries are stripClusters, rest are pixels.
      // later we  must pass nStrip in metadata or as template param; here we assume fixed
      uint32_t nStrip = in.metadata().size() - in.metadata().pixelCount;
      for (auto idx : elements_with_stride(acc, n)) {
        uint32_t word, x,y,w;
        bool seed = false;
        if (idx < nStrip) {
          // strip-on-PS
          word = in[idx].stripClustersWords();
          uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
          uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS))
                          & SCLUSTER_ADDRESS_PS_MAX_VALUE;
          w = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS))
              & WIDTH_MAX_VALUE;
          if (w==0) w = 8;
          seed = false;
          x = STRIPS_PER_SSA * chip + addr;
          y = (idx%2==0 ? 0u : 1u);
        } else {
          // pixel-on-PS
          word = in[idx].pixelClustersWords();
          uint32_t chip = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
          uint32_t addr = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS))
                          & SCLUSTER_ADDRESS_PS_MAX_VALUE;
          w = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS))
              & WIDTH_MAX_VALUE;
          if (w==0) w = 8;
          seed = true;
          x = STRIPS_PER_SSA * chip + addr;
          uint32_t z = word & PS_Z_BITS_MASK;
          y = (idx%2==0 ? z : (z + 16));
        }
        out[idx].x()        = x;
        out[idx].y()        = y;
        out[idx].width()    = w;
        out[idx].seedFlag() = seed;
      }
    }
  };
*/

  //setter fucntion example. jEREMI'S CODE 
 /*
  void HGCalRecHitCalibrationAlgorithms::loadCalibParams(CalibParams& newCalibParams) {
    LogDebug("HGCalRecHitCalibrationAlgorithms") << "\nINFO -- HGCalRecHitCalibrationAlgorithms::loadCalibParams: " << newCalibParams.size() << " elements" << std::endl;
    calibParams = newCalibParams;
  }
*/

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE


