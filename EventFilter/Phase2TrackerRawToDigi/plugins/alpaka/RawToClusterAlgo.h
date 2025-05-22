#ifndef EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h
#define EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h" 

#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  void launchS2UnpackerKernel(
      Queue& queue,
      const Phase2RawToCluster::StripPixelDeviceCollection& devStripPixel,
      Phase2RawToCluster::ClusterPropDeviceCollection& devClusterProp,
      size_t totalWords);

  // Future kernels:
  // void launchPSUnpackerKernel(...);
  // void launchClusterizer(...);
  void launchUnpacker(
        Queue& queue,
        cms::alpakatools::device_buffer<Device, unsigned char[]> rawdatabuff,
        cms::alpakatools::device_buffer<Device, size_t[]> sizedatabuff,
        cms::alpakatools::device_buffer<Device, size_t[]> offsetdatabuff,
	 cms::alpakatools::device_buffer<Device, int[]> inmap);
       	}

#endif  // EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h

