#ifndef EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h
#define EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropDeviceCollection.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerSpecifications.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2DAQFormatSpecification.h"

using namespace Phase2TrackerSpecifications;
using namespace Phase2DAQFormatSpecification;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

    void launchUnpacker(
            Queue& queue,
            cms::alpakatools::device_buffer<Device, unsigned char[]> rawdatabuff,
            cms::alpakatools::device_buffer<Device, size_t[]> sizedatabuff,
            cms::alpakatools::device_buffer<Device, size_t[]> offsetdatabuff,
            cms::alpakatools::device_buffer<Device, int[]> inmap,
            cms::alpakatools::device_buffer<Device, int[]> detIdMap, // Add detIdMap
            cms::alpakatools::device_buffer<Device, std::pair<int, int>[]> stackMap, // add stackMap
            uint32_t stackMapSize,  // Add size parameter
            Phase2RawToCluster::ClusterPropDeviceCollection::View out,
            uint32_t* globalCounter
    );

} // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif // EventFilter_Phase2TrackerRawToDigi_RawToClusterAlgo_h