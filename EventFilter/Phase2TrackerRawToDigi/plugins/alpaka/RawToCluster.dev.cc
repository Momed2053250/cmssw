// ================================ GPU kernels File ================================
// alpaka-related imports

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropDeviceCollection.h"

#include "EventFilter/Phase2TrackerRawToDigi/interface/SensorHybrid.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerSpecifications.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2DAQFormatSpecification.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/TrackerHeader.h"

using namespace cms::alpakatools;
using namespace Phase2RawToCluster;
using namespace Phase2TrackerSpecifications;
using namespace Phase2DAQFormatSpecification;
using namespace ALPAKA_ACCELERATOR_NAMESPACE;

// Debug flag
//#define Debug_GPU

// ------------1) Define max local scratch sizes ------------------
// Per-thread scratch arrays (safe on GPU backends)
static constexpr int MaxOffsetWords    = (OFFSET_BITS * CICs_PER_SLINK) / N_BITS_PER_WORD;
static constexpr int MaxStripClusters  = N_CLUSTER_MASK + 1;   // 128
static constexpr int MaxPixelClusters  = N_CLUSTER_MASK + 1;   // 128
static constexpr int MaxPayloadLines =
  ((MaxStripClusters * SS_CLUSTER_BITS + MaxPixelClusters * PX_CLUSTER_BITS) / N_BITS_PER_WORD) + 1;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // create masking
  ALPAKA_FN_ACC inline int createMask(int nBits) { return (1 << nBits) - 1; }

  // Read a 32bit word from a byte buffer
  ALPAKA_FN_ACC inline uint32_t readLine(const unsigned char* dataPtr, int byteIdx) {
    return (static_cast<uint32_t>(dataPtr[byteIdx])     << 24) |
           (static_cast<uint32_t>(dataPtr[byteIdx + 1]) << 16) |
           (static_cast<uint32_t>(dataPtr[byteIdx + 2]) << 8)  |
            static_cast<uint32_t>(dataPtr[byteIdx + 3]);
  }

  // Compute byte offset within payload: skip header and channel offset table
  ALPAKA_FN_ACC inline int getLineIndex(int channelIdx, unsigned int iline) {
    return channelIdx + N_BYTES_PER_WORD + iline * N_BYTES_PER_WORD;
  }

  // Extract cluster words across multiple lines with bit-packing
  ALPAKA_FN_ACC inline void readPayload(
      uint32_t* clusterWords,         // output array for extracted words
      const uint32_t* lines,          // input buffer of 32-bit words
      int numClusters,
      int& nAvailableBits,            // bits left in current "line"
      int& iLine,                     // current line index
      int& bitsToRead,                // leftover bits to read if cluster spans words
      int& nFullClusters,             // full clusters read in current line
      const int clusterBits,
      const int clusterWordMask,      // mask to isolate cluster bits
      const bool isPixelCluster,
      int nFullClustersStrips = 0     // count of strip clusters in PS module
  ) {
    for (int icluster = 0; icluster < numClusters; ++icluster) {
      if (nAvailableBits >= clusterBits) {
        // calculate the shift to align bits for extraction as in the CPU code
        int shift = N_BITS_PER_WORD - bitsToRead - (nFullClusters + 1) * clusterBits;
        // account for bits used by last strip cluster (PS only)
        if (icluster == 0 && isPixelCluster) shift -= (nFullClustersStrips) * SS_CLUSTER_BITS;
        nFullClustersStrips = 0; // reset

        // mask, and save cluster word
        clusterWords[icluster] = (lines[iLine] >> shift) & clusterWordMask;

        // update available bits and number of full clusters from this line
        nAvailableBits -= clusterBits;
        nFullClusters++;

        // Advance to next word if we consumed all bits
        if (nAvailableBits == 0) {
          ++iLine;
          nAvailableBits = N_BITS_PER_WORD;
          nFullClusters = 0;
          bitsToRead = 0;
        }
      } else {
        // cluster spans word boundary
        const int nMask = createMask(nAvailableBits);
        const uint16_t wordLeft = static_cast<uint16_t>(lines[iLine] & nMask);

        bitsToRead = clusterBits - nAvailableBits;
        const int nextMask = createMask(bitsToRead);
        const uint16_t wordRight = static_cast<uint16_t>((lines[iLine + 1] >> (N_BITS_PER_WORD - bitsToRead)) & nextMask);

        clusterWords[icluster] = (static_cast<uint32_t>(wordLeft) << bitsToRead) | wordRight;

        // prepare for next read
        nAvailableBits = N_BITS_PER_WORD - bitsToRead;
        ++iLine;
        nFullClusters = 0;
      }
    }
  }

  // maximum total clusters
  static constexpr size_t MaxTotalClusters =
      (N_CLUSTER_MASK + 1) * CICs_PER_SLINK * (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;

  // --------------------------- Unpacker kernel ---------------------------
  struct Unpacker {
    template <
      typename Acc,
      typename RawBufView,
      typename SizeBufView,
      typename OffBufView,
      typename ModuleTypeView,
      typename InnerDetIdView,   // uint32_t view
      typename OuterDetIdView,   // uint32_t view
      typename OutView
    >
    ALPAKA_FN_ACC void operator()(Acc const& acc,
                                  RawBufView in,
                                  SizeBufView sizes,
                                  OffBufView offsets,
                                  ModuleTypeView const& detIdxModuleTypeMap,
                                  InnerDetIdView const& innerDetIdForFlatIdx,
                                  OuterDetIdView const& outerDetIdForFlatIdx,
                                  OutView out,
                                  uint32_t* globalCounter) const {
      // Per-thread scratch (safe on GPU)
      uint32_t offsetWords[MaxOffsetWords];
      uint32_t lines[MaxPayloadLines];
      uint32_t stripClusterWords[MaxStripClusters];
      uint32_t pixelClusterWords[MaxPixelClusters];

      // Global linear thread id and stride
      const uint32_t gtid = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
      const uint32_t gdim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc)[0u];

      // Number of SLINK fragments we actually process
      const uint32_t NSlinks = (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;

      // Each thread processes multiple SLINKs in a grid-stride loop
      for (uint32_t frdId = gtid; frdId < NSlinks; frdId += gdim) {
        // Skip empty fragments
        const uint32_t fragSizeBytes = static_cast<uint32_t>(sizes[frdId]);
        if (fragSizeBytes == 0) continue;

        const unsigned char* dataPtr = in + offsets[frdId];

        // 0) fragment minimal size check for header+offsets
        const uint32_t minHdrOff = static_cast<uint32_t>((HEADER_N_LINES + MODULES_PER_SLINK) * N_BYTES_PER_WORD);
        if (fragSizeBytes < minHdrOff) continue;

        // 1) Offsets start after the fixed-size header (we don't need header contents for unpacking)
        const size_t nOffsetsLines = MaxOffsetWords; // (OFFSET_BITS * CICs_PER_SLINK) / 32
        const size_t initByte      = HEADER_N_LINES * N_BYTES_PER_WORD;

        // Ensure last offset word is within fragment
        const size_t lastOffByte = initByte + (nOffsetsLines - 1) * N_BYTES_PER_WORD + (N_BYTES_PER_WORD - 1);
        if (lastOffByte >= fragSizeBytes) continue;

        for (size_t k = 0; k < nOffsetsLines; ++k) {
          const int byteIdx = static_cast<int>(initByte + k * N_BYTES_PER_WORD);
          offsetWords[k] = readLine(dataPtr, byteIdx);
        }

        // 2) Unpack each channel (same order and logic as CPU)
        for (unsigned int iChannel = 0; iChannel < CICs_PER_SLINK; ++iChannel) {
          // Build flatIdx = frdId * CICs + iChannel
          const unsigned flatIdx = frdId * CICs_PER_SLINK + iChannel;

          // Retrieve module type (0=undef, 1=TwoS, 2=PS)
          const int moduleType = detIdxModuleTypeMap[flatIdx];
          if (moduleType == 0) continue; // skip unconnected

          const bool is2SModule = (moduleType == 1);

          // Compute byte index of channel header
          const size_t offsetTableStart = (HEADER_N_LINES + MODULES_PER_SLINK) * N_BYTES_PER_WORD;

          // Read 16-bit offset for this channel (same as CPU)
          const int wordIdx = static_cast<int>(iChannel / 2);
          const uint16_t channelOffset16 = (iChannel % 2 == 0)
            ? static_cast<uint16_t>(offsetWords[wordIdx] & 0xFFFFu)
            : static_cast<uint16_t>(offsetWords[wordIdx] >> 16);

          const uint32_t idx = static_cast<uint32_t>(offsetTableStart + channelOffset16 * N_BYTES_PER_WORD);

          // header is 4 bytes at idx..idx+3
          if (idx + (N_BYTES_PER_WORD - 1) >= fragSizeBytes) {
            // bogus channel offset -> skip channel
            continue;
          }
          const uint32_t chHeaderWord = readLine(dataPtr, static_cast<int>(idx));

          const unsigned int numStripClusters =
            (chHeaderWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS - N_STRIP_CLUSTER_BITS)) & N_CLUSTER_MASK;
          const unsigned int numPixelClusters = chHeaderWord & N_CLUSTER_MASK;

          // Define number of payload lines
          unsigned int nLines = 0;
          if (numStripClusters + numPixelClusters > 0) {
            const unsigned int neededBits =
              numStripClusters * SS_CLUSTER_BITS + numPixelClusters * PX_CLUSTER_BITS;
            nLines = static_cast<unsigned int>(neededBits / N_BITS_PER_WORD) + 1u;
          }
          if (nLines > MaxPayloadLines) nLines = MaxPayloadLines;

          // ensure the last payload word fits
          if (nLines > 0) {
            const uint32_t lastPayloadByte = static_cast<uint32_t>(
              getLineIndex(static_cast<int>(idx), nLines - 1) + (N_BYTES_PER_WORD - 1));
            if (lastPayloadByte >= fragSizeBytes) {
              // malformed payload -> skip channel
              continue;
            }
          }

          // Retrieve payload lines
          for (unsigned int k = 0; k < nLines; ++k) {
            const int byteIdx = getLineIndex(static_cast<int>(idx), k);
            lines[k] = readLine(dataPtr, byteIdx);
          }

          // Read payloads (bit-unpack into per-thread scratch)
          int nAvailableBits = N_BITS_PER_WORD;
          int iLine = 0;
          int bitsToRead = 0;
          int nFullClustersStrip = 0;
          int nFullClustersPix = 0;

          const unsigned int useStrip = (numStripClusters <= static_cast<unsigned int>(MaxStripClusters)) ? numStripClusters : static_cast<unsigned int>(MaxStripClusters);
          const unsigned int usePixel = (numPixelClusters <= static_cast<unsigned int>(MaxPixelClusters)) ? numPixelClusters : static_cast<unsigned int>(MaxPixelClusters);

          // 2S or PS: strips first
          if (useStrip > 0) {
            readPayload(stripClusterWords, lines, static_cast<int>(useStrip),
                        nAvailableBits, iLine, bitsToRead, nFullClustersStrip,
                        SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false);
          }

          // PS only: then pixels
          if (!is2SModule && usePixel > 0) {
            readPayload(pixelClusterWords, lines, static_cast<int>(usePixel),
                        nAvailableBits, iLine, bitsToRead, nFullClustersPix,
                        PX_CLUSTER_BITS, PX_CLUSTER_WORD_MASK, true, nFullClustersStrip);
          }

          // Reserve output slots; guard against capacity overflow
          const uint32_t want = is2SModule ? useStrip : (useStrip + usePixel);
          if (want == 0) continue;

          const uint32_t base = alpaka::atomicAdd(acc, globalCounter, want);
          const uint32_t cap  = static_cast<uint32_t>(MaxTotalClusters);
          if (base >= cap) {
            // counter exceeded capacity; host will clamp anyway
            continue;
          }
          const uint32_t room = cap - base;

          // how many can we actually write
          const uint32_t takeStrip = is2SModule
            ? (useStrip > room ? room : useStrip)
            : (useStrip > room ? room : useStrip);
          const uint32_t takePix   = (!is2SModule && takeStrip < room)
            ? (usePixel > (room - takeStrip) ? (room - takeStrip) : usePixel)
            : 0u;

          // Unpack to output SoA using precomputed inner/outer detIds
          const uint32_t innerDet = innerDetIdForFlatIdx[flatIdx];
          const uint32_t outerDet = outerDetIdForFlatIdx[flatIdx];

          // 2S: strips only
          if (is2SModule) {
            for (uint32_t ic = 0; ic < takeStrip; ++ic) {
              const uint32_t word = stripClusterWords[ic];
              const uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
              const uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_ONLY_BITS_2S)) & SCLUSTER_ADDRESS_MASK;
              const bool     seed = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_2S)) & IS_SEED_SENSOR_MASK;
              uint32_t       w    = word & WIDTH_MAX_VALUE;
              if (w == 0) w = 8;

              const uint32_t outIdx = base + ic;
              if (outIdx >= cap) break;

              out[outIdx].strip()     = STRIPS_PER_CBC * chip + addr;           // x
              out[outIdx].row()       = (iChannel % 2 == 0) ? 0u : 1u;          // y
              out[outIdx].size()      = w;                                      // width
              out[outIdx].threshold() = seed;                                   // seedFlag
              out[outIdx].mipBit()    = 0;
              out[outIdx].column()    = 0;
              out[outIdx].edge()      = 0;
              out[outIdx].detId()     = seed ? innerDet : outerDet;             // inner=seed, outer=corr
            }
          } else {
            // PS: strips (correlated sensor  outer)
            for (uint32_t ic = 0; ic < takeStrip; ++ic) {
              const uint32_t word = stripClusterWords[ic];
              const uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
              const uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS)) & SCLUSTER_ADDRESS_PS_MAX_VALUE;
              uint32_t       w    = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS)) & WIDTH_MAX_VALUE;
              const uint32_t mip  = word & MIP_BITS_MASK;
              if (w == 0) w = 8;

              const uint32_t outIdx = base + ic;
              if (outIdx >= cap) break;

              out[outIdx].strip()     = STRIPS_PER_SSA * chip + addr;           // x
              out[outIdx].row()       = (iChannel % 2 == 0) ? 0u : 1u;          // y
              out[outIdx].size()      = w;
              out[outIdx].threshold() = false;                                  // strip on PS is correlated sensor
              out[outIdx].mipBit()    = mip;
              out[outIdx].column()    = 0;
              out[outIdx].edge()      = 0;
              out[outIdx].detId()     = outerDet;                               // outer (correlated)
            }

            // PS: pixels (seed  inner); placed immediately after strips
            for (uint32_t ic = 0; ic < takePix; ++ic) {
              const uint32_t word = pixelClusterWords[ic];
              const uint32_t chip = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
              const uint32_t addr = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS)) & SCLUSTER_ADDRESS_PS_MAX_VALUE;
              uint32_t       w    = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS)) & WIDTH_MAX_VALUE;
              const uint32_t z    = word & PS_Z_BITS_MASK;
              if (w == 0) w = 8;

              const uint32_t outIdx = base + takeStrip + ic;
              if (outIdx >= cap) break;

              out[outIdx].strip()     = STRIPS_PER_SSA * chip + addr;           // x
              out[outIdx].row()       = (iChannel % 2 == 0) ? z : (z + 16);     // y
              out[outIdx].size()      = w;
              out[outIdx].threshold() = true;                                   // pixel on PS is seed
              out[outIdx].mipBit()    = 0;
              out[outIdx].column()    = z;                                      // keep as in your mapping
              out[outIdx].edge()      = 0;
              out[outIdx].detId()     = innerDet;                               // inner (seed)
            }
          }
        } // channels
      }   // frdId stride loop
    }     // operator()
  };

  // Launch the generic Unpacker kernel (header + payload) on device
  void launchUnpacker(
      Queue& queue,
      cms::alpakatools::device_buffer<Device, unsigned char[]> rawdatabuff,
      cms::alpakatools::device_buffer<Device, size_t[]>        sizedatabuff,
      cms::alpakatools::device_buffer<Device, size_t[]>        offsetdatabuff,
      cms::alpakatools::device_buffer<Device, int[]>           detIdxModuleTypeDevice,
      cms::alpakatools::device_buffer<Device, uint32_t[]>      innerDetIdDevice,  // uint32_t
      cms::alpakatools::device_buffer<Device, uint32_t[]>      outerDetIdDevice,  // uint32_t
      Phase2RawToCluster::ClusterPropDeviceCollection::View out,
      uint32_t* globalCounter) {

    const uint32_t NSlinks = (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;
    const uint32_t threadsPerBlock = 128;
    const uint32_t blocks = (NSlinks + threadsPerBlock - 1) / threadsPerBlock;

    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

    alpaka::exec<Acc1D>(
      queue,
      workDiv,
      Unpacker{},
      rawdatabuff.data(),
      sizedatabuff.data(),
      offsetdatabuff.data(),
      detIdxModuleTypeDevice.data(),
      innerDetIdDevice.data(),
      outerDetIdDevice.data(),
      out,
      globalCounter
    );
  }

} // namespace ALPAKA_ACCELERATOR_NAMESPACE

// No dynamic shared memory is required
namespace alpaka::trait {
  template<>
  struct BlockSharedMemDynSizeBytes<Unpacker, Acc1D> {
    template<
      typename RawBufView,
      typename SizeBufView,
      typename OffBufView,
      typename ModuleTypeView,
      typename InnerDetIdView,
      typename OuterDetIdView,
      typename OutView
    >
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
      Unpacker const&,
      Vec1D /*threads*/,
      Vec1D /*elements*/,
      RawBufView, SizeBufView, OffBufView,
      ModuleTypeView, InnerDetIdView, OuterDetIdView,
      OutView,
      uint32_t* /*globalCounter*/
    ) {
      return 0u;
    }
  };
}