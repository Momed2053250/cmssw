#ifndef EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h
#define EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h" 

#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/alpaka/ClusterPropDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  void launchS2UnpackerKernel(
      Queue& queue,
      const Phase2RawToCluster::StripPixelDeviceCollection& devStripPixel,
      Phase2RawToCluster::ClusterPropDeviceCollection& devClusterProp,
      size_t totalWords);

  // Future kernels:
  // void launchPSUnpackerKernel(...);
  // void launchClusterizer(...);

}

#endif  // EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h

