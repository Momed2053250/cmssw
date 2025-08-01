// GPU kernals File 
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

// ------------1) Define max sharedmem sizes ------------------
static constexpr int MaxHeaderWords    = HEADER_N_LINES;
static constexpr int MaxOffsetWords    = (OFFSET_BITS * CICs_PER_SLINK) / N_BITS_PER_WORD;
static constexpr int MaxStripClusters  = N_CLUSTER_MASK + 1;   // 128
static constexpr int MaxPixelClusters  = N_CLUSTER_MASK + 1;   // 128
static constexpr int MaxPayloadLines =
((MaxStripClusters * SS_CLUSTER_BITS +
  MaxPixelClusters * PX_CLUSTER_BITS)
 / N_BITS_PER_WORD) + 1;
// total 32-bit words we need in shared mem:
static constexpr int MaxTotalSharedWords =
MaxHeaderWords
+ MaxOffsetWords
+ MaxPayloadLines
+ MaxStripClusters
+ MaxPixelClusters;



namespace ALPAKA_ACCELERATOR_NAMESPACE {

	// create masking 
	ALPAKA_FN_ACC int createMask(int nBits) {
		return (1 << nBits) - 1;
	}	
	// Read a 32bit word from a byte buffer
	ALPAKA_FN_ACC uint32_t readLine(const unsigned char* dataPtr, int lineIdx){						
		uint32_t line = (static_cast<uint32_t>(dataPtr[lineIdx]) << 24) | 
			(static_cast<uint32_t>(dataPtr[lineIdx + 1]) << 16) | 
			(static_cast<uint32_t>(dataPtr[lineIdx + 2]) << 8) | 
			(static_cast<uint32_t>(dataPtr[lineIdx + 3]));

		return line;                                
	}
	// Compute byte offset within payload: skip header and channel offset table
	ALPAKA_FN_ACC int getLineIndex(int channelIdx, unsigned int iline){
		return channelIdx + N_BYTES_PER_WORD + iline * N_BYTES_PER_WORD; 
	}

	// Extract cluster words across multiple lines with bit-packing
	ALPAKA_FN_ACC void readPayload(uint32_t* clusterWords,   // output array for extracted words
			uint32_t* lines,                         // input buffer of 32-bit words
			int numClusters,
			int& nAvailableBits,                    // bits left in current "line"
			int& iLine,                             // current line index
			int& bitsToRead,                         // leftover bits to read if cluster spans words
			int& nFullClusters,                     // full clusters read in current line
			int clusterBits,
			int clusterWordMask,                    // mask to isolate cluster bits
			bool isPixelCluster,
			int nFullClustersStrips = 0             // count of strip clusters in PS module !!Perhaps 
			)
	{
		for (int icluster = 0; icluster < numClusters; icluster++) {
			if (nAvailableBits >= clusterBits) {
				// calculate the shift to align bits for extraction as in the code (CPU based) 
				int shift = N_BITS_PER_WORD - bitsToRead - (nFullClusters + 1) * clusterBits;
				// take into account bits already used for the last strip cluster
				if (icluster == 0 && isPixelCluster) 
					// adjust for prior strip clusters in PS modules
					shift -= (nFullClustersStrips)* SS_CLUSTER_BITS;
				nFullClustersStrips = 0; // reset

				// mask, and save cluster word
				clusterWords[icluster] = (lines[iLine] >> shift) & clusterWordMask;
				// update available bits and number of full clusters from this line
				nAvailableBits -= clusterBits;
				nFullClusters++;

				// Advance to next "line" if we've consumed all bits
				if (nAvailableBits == 0) {
					iLine++;
					nAvailableBits = N_BITS_PER_WORD;
					nFullClusters = 0;
					bitsToRead = 0;
				}
			} else {

				//Handle clusters spanning across two 32bit words 
				// get the remaining bits from the current line. first create the mask, then mask
				int nMask = createMask(nAvailableBits);
				uint16_t wordLeft = lines[iLine] & nMask;

				// create mask for next line
				bitsToRead = clusterBits - nAvailableBits;
				int nextMask = createMask(bitsToRead);
				// shift and mask
				uint16_t wordRight = (lines[iLine + 1] >> (N_BITS_PER_WORD - bitsToRead)) & nextMask;

				// compose the full cluster word
				clusterWords[icluster] = (wordLeft << bitsToRead) | wordRight;

				// re-set n available bits
				nAvailableBits = N_BITS_PER_WORD - bitsToRead;
				// advance by one line and re-init the number of complete clusters read from the current line
				iLine++;
				nFullClusters = 0;

			}
		}
	}

	// Read 16-bit offset for a given channel from packed offsetWords array
	ALPAKA_FN_ACC uint16_t getOffsetForChannel(unsigned int iChannel, uint32_t* offsetWords) {
		if (iChannel >= CICs_PER_SLINK) {
			printf("Error: iChannel %u too high\n", iChannel);
			return 0;
		}
		//TODO:: Optimize
		// Even channel: lower 16 bits of word iChannel/2
		int wordIdx = iChannel / 2;
		if (iChannel % 2 == 0) {
			return static_cast<uint16_t>(offsetWords[wordIdx] & 0xFFFF);
		} else {
			// extract the upper 16 bits by shifting right by 16
			return static_cast<uint16_t>(offsetWords[wordIdx] >> 16); 
		}
	}

	// maximum total clusters 
	static constexpr size_t MaxTotalClusters = (N_CLUSTER_MASK + 1) * CICs_PER_SLINK * (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC;
	// Unpacker kernel: top-level device loop over FED fragments
	struct Unpacker {
		template <
			typename Acc,
				 typename RawBufView,
				 typename SizeBufView,
				 typename OffBufView,
				 typename InMapView,
				 typename DetIdMapView, // Add DetIdMapView
				 typename StackMapView, // add stack map 
				 typename OutView
					 >
					 ALPAKA_FN_ACC void operator()(
							 Acc const& acc,
							 RawBufView in,
							 SizeBufView sizes,
							 OffBufView offsets,
							 InMapView const& detIdxModuleTypeMap,
							 DetIdMapView const& detIdMap, // Add detIdMap
							 StackMapView const& stackMap, // add stack map 
							 uint32_t stackMapSize,  // Add stackMapSize parameter
							 OutView out,
							 uint32_t* globalCounter
							 ) const {
						 // (A) Allocate ONE contiguous chunk of dynamic shared memory:
						 uint8_t* smemBytes = alpaka::getDynSharedMem<uint8_t>(acc);
						 uint32_t* smemWords = reinterpret_cast<uint32_t*>(smemBytes);
#ifdef Debug_GPU
						 if (!smemBytes) {
							 printf("Error: smemBytes is null\n");
							 return;
						 }
#endif
						 // (B) Slice that chunk into disjoint regions:
						 uint32_t* headerWords = smemWords;
						 uint32_t* offsetWords = headerWords + MaxHeaderWords;
						 uint32_t* lines = offsetWords + MaxOffsetWords;
						 uint32_t* stripClusterWords = lines + MaxPayloadLines;
						 uint32_t* pixelClusterWords = stripClusterWords + MaxStripClusters;

						 // Track starting index for each channel pair
						 uint32_t channelPairStartIdx = 0;
						 if (cms::alpakatools::once_per_block(acc)) {
							 channelPairStartIdx = alpaka::atomicAdd(acc, globalCounter, static_cast<uint32_t>(MaxStripClusters + MaxPixelClusters));
						 }
						 alpaka::syncBlockThreads(acc);

						 // Local counter for clusters within this channel pair
						 uint32_t localClusterIdx = 0;

						 // Iterate over each FED fragment ID in parallel
						 for (auto frdId : cms::alpakatools::independent_groups(acc, (MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC)) {
							 if (sizes[frdId]  == 0 ) continue;  // Skip empty fragments
							 // {
								 const unsigned char* dataPtr = in + offsets[frdId];

								 // 1) Read the header
								 size_t nHeaderLines = HEADER_N_LINES;
								 for (auto k : cms::alpakatools::independent_group_elements(acc, nHeaderLines)) {
									 auto byteIdx = k * N_BYTES_PER_WORD;
									 headerWords[k] = readLine(dataPtr, byteIdx);
								 }
								 alpaka::syncBlockThreads(acc);
#ifdef Debug_GPU
								 printf("headerWords[0] = %u\n", headerWords[0]);
#endif

								 // 2) Read offset words
								 size_t nOffsetsLines = (OFFSET_BITS * CICs_PER_SLINK) / N_BITS_PER_WORD;
								 size_t initByte = HEADER_N_LINES * N_BYTES_PER_WORD;
								 for (auto k : cms::alpakatools::independent_group_elements(acc, nOffsetsLines)) {
									 int byteIdx = static_cast<int>(initByte + k * N_BYTES_PER_WORD);
									 offsetWords[k] = readLine(dataPtr, byteIdx);
								 }
								 alpaka::syncBlockThreads(acc);
#ifdef Debug_GPU
								 printf("offsetWords[0] = %u\n", offsetWords[0]);
#endif

								 // Unpack each channel
								 for (unsigned int iChannel = 0; iChannel < CICs_PER_SLINK; iChannel++) {
									 // Reset local cluster index for even channels
									 if (iChannel % 2 == 0) {
										 localClusterIdx = 0;
									 }

									 // Retrieve module type
									 const unsigned CICs = CICs_PER_SLINK;
									 unsigned flatIdx = frdId * CICs + iChannel;
									 int thisDetId = detIdMap[flatIdx]; // Get the detId from the map
									 int is2SModule = detIdxModuleTypeMap[flatIdx] == 1 ? 1 : 0;
#ifdef Debug_GPU
									 if (is2SModule != 0) {
										 printf("is2SModule is: %d\n", is2SModule);
									 }
#endif
									 // Compute byte index of channel header
									 size_t offsetTableStart = (HEADER_N_LINES + MODULES_PER_SLINK) * N_BYTES_PER_WORD;
									 int channelOffset16 = static_cast<int>(getOffsetForChannel(iChannel, offsetWords));
#ifdef Debug_GPU
									 printf("ChannelOffset16 is: %u\n", channelOffset16);
#endif
									 int idx = static_cast<int>(offsetTableStart + channelOffset16 * N_BYTES_PER_WORD);
#ifdef Debug_GPU
									 printf("idx is: %u\n", idx);
#endif

									 // Read channel header and extract cluster counts
									 uint32_t chHeaderWord = readLine(dataPtr, idx);
									 unsigned int numStripClusters =
										 (chHeaderWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS - N_STRIP_CLUSTER_BITS)) & N_CLUSTER_MASK;
									 unsigned int numPixelClusters = chHeaderWord & N_CLUSTER_MASK;

									 // Define number of payload lines
									 unsigned int nLines = (numStripClusters + numPixelClusters > 0) ?
										 int((numStripClusters * SS_CLUSTER_BITS + numPixelClusters * PX_CLUSTER_BITS) / N_BITS_PER_WORD) + 1 : 0;
#ifdef Debug_GPU
									 printf("n strip clusters are: %u\n", numStripClusters);
									 printf("n pixel clusters are: %u\n", numPixelClusters);
#endif

									 // Retrieve payload lines
									 for (auto k : cms::alpakatools::independent_group_elements(acc, nLines)) {
										 int byteIdx = getLineIndex(idx, k);
										 lines[k] = readLine(dataPtr, byteIdx);  //RACE DETECTED  
										 // print the lines
#ifdef Debug_GPU
										 if (k == 0)  // Match CPU: only print the first line
											 printf("Lines[0] = %u\n", lines[0]);
#endif
									 }
									 alpaka::syncBlockThreads(acc);
									 // Read payloads
									 int nAvailableBits = N_BITS_PER_WORD;
									 int iLine = 0;
									 int bitsToRead = 0;
									 int nFullClustersStrip = 0;
									 int nFullClustersPix = 0;

									 if (is2SModule) {
										 if (cms::alpakatools::once_per_block(acc)) {
											 readPayload(stripClusterWords, lines, numStripClusters, nAvailableBits, iLine, bitsToRead,
													 nFullClustersStrip, SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false);
										 }
									 } else {
										 if (cms::alpakatools::once_per_block(acc)) {
											 readPayload(stripClusterWords, lines, numStripClusters, nAvailableBits, iLine, bitsToRead,
													 nFullClustersStrip, SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false);
											 // print out the strip cluster words
#ifdef Debug_GPU
											 printf("Strip Cluster words: \n");
											 for (unsigned int i = 0; i < numStripClusters; ++i) {
												 printf("%u ", stripClusterWords[i]);
											 }
											 printf("\n");
#endif
											 readPayload(pixelClusterWords, lines, numPixelClusters, nAvailableBits, iLine, bitsToRead,
													 nFullClustersPix, PX_CLUSTER_BITS, PX_CLUSTER_WORD_MASK, true, nFullClustersStrip);
											 // print out the pixel cluster words
#ifdef Debug_GPU
											 printf("Pixel Cluster words: \n");
											 for (unsigned int i = 0; i < numPixelClusters; ++i) {
												 printf("%u ", pixelClusterWords[i]);		
											 }
											 printf("\n");
#endif
										 }
									 }
									 alpaka::syncBlockThreads(acc);

									 // Unpack clusters and store in output SoA
									 if (is2SModule) {
										#ifdef Debug_GPU
											printf("Unpacking for 2S module\n");
										#endif
											for (auto icluster : cms::alpakatools::independent_group_elements(acc, numStripClusters)) {
												uint32_t word = stripClusterWords[icluster];
												uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
												uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_ONLY_BITS_2S)) & SCLUSTER_ADDRESS_MASK;
												bool seed = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_2S)) & IS_SEED_SENSOR_MASK;
												uint32_t w = word & WIDTH_MAX_VALUE;
												if (w == 0) w = 8;
										
												uint32_t outIdx = channelPairStartIdx + localClusterIdx;
												if (outIdx < MaxTotalClusters) {
													out[outIdx].strip() = STRIPS_PER_CBC * chip + addr;  // Maps to x
													out[outIdx].row() = iChannel % 2 == 0 ? 0u : 1u;     // Maps to y
													out[outIdx].size() = w;                              // Maps to width
													out[outIdx].threshold() = seed;                      // Maps to seedFlag
													out[outIdx].mipBit() = 0;                            // Unchanged
													out[outIdx].column() = 0;                            // Initialize to 0 (unused)
													out[outIdx].edge() = 0;                              // Initialize to 0 (unused)
													if (thisDetId >= 0 && thisDetId < static_cast<int>(stackMapSize)){
														out[outIdx].detId() = seed ? stackMap[thisDetId].first : stackMap[thisDetId].second;
													}
														#ifdef Debug_GPU
													printf("Unpacked values 2S: chipID = %u, addr = %u, size = %u, threshold = %d, strip = %u, row = %u, column = %u, edge = %u\n",
														   chip, addr, w, seed, out[outIdx].strip(), out[outIdx].row(), out[outIdx].column(), out[outIdx].edge());
										#endif
												}
												localClusterIdx++;
											}
										} else {
											// PS strip clusters
											for (auto icluster : cms::alpakatools::independent_group_elements(acc, numStripClusters)) {
												uint32_t word = stripClusterWords[icluster];
												uint32_t chip = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
												uint32_t addr = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS)) & SCLUSTER_ADDRESS_PS_MAX_VALUE;
												uint32_t w = (word >> (SS_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS)) & WIDTH_MAX_VALUE;
												uint32_t mipBit = word & MIP_BITS_MASK;
												if (w == 0) w = 8;
										
												uint32_t outIdx = channelPairStartIdx + localClusterIdx;
												if (outIdx < MaxTotalClusters) {
													out[outIdx].strip() = STRIPS_PER_SSA * chip + addr;  // Maps to x
													out[outIdx].row() = iChannel % 2 == 0 ? 0u : 1u;     // Maps to y
													out[outIdx].size() = w;                              // Maps to width
													out[outIdx].threshold() = false;                     // Maps to seedFlag
													out[outIdx].mipBit() = mipBit;                       // Unchanged
													out[outIdx].column() = 0;                            // Initialize to 0 (unused)
													out[outIdx].edge() = 0;                              // Initialize to 0 (unused)
													//out[outIdx].detId() = stackMap[thisDetId].second; // outer (correlated sensor)
													if (thisDetId >= 0 && thisDetId < static_cast<int>(stackMapSize)) {
														out[outIdx].detId() = stackMap[thisDetId].second;
													}
													//else {
													//	out[outIdx].detId() = -1; // Invalid detId
													//}
										#ifdef Debug_GPU
													printf("Unpacked values S on PS: chipID = %u, addr = %u, size = %u, mipBit = %u, strip = %u, row = %u, column = %u, edge = %u\n",
														   chip, addr, w, mipBit, out[outIdx].strip(), out[outIdx].row(), out[outIdx].column(), out[outIdx].edge());
										#endif
												}
												localClusterIdx++;
											}
											// PS pixel clusters
											for (auto icluster : cms::alpakatools::independent_group_elements(acc, numPixelClusters)) {
												uint32_t word = pixelClusterWords[icluster];
												uint32_t chip = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS)) & CHIP_ID_MAX_VALUE;
												uint32_t addr = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS)) & SCLUSTER_ADDRESS_PS_MAX_VALUE;
												uint32_t w = (word >> (PX_CLUSTER_BITS - CHIP_ID_BITS - SCLUSTER_ADDRESS_BITS_PS - WIDTH_BITS)) & WIDTH_MAX_VALUE;
												uint32_t z = word & PS_Z_BITS_MASK;
												if (w == 0) w = 8;
										
												uint32_t outIdx = channelPairStartIdx + localClusterIdx;
												if (outIdx < MaxTotalClusters) {
													out[outIdx].strip() = STRIPS_PER_SSA * chip + addr;  // Maps to x
													out[outIdx].row() = iChannel % 2 == 0 ? z : (z + 16); // Maps to y
													out[outIdx].size() = w;                              // Maps to width
													out[outIdx].threshold() = true;                      // Maps to seedFlag
													out[outIdx].mipBit() = 0;                            // Unchanged
													out[outIdx].column() = z;                            // Check if set to z is correct for the PS Modules 
													out[outIdx].edge() = 0;                              // Initialize to 0 (unused)
													//out[outIdx].detId() = stackMap[thisDetId].first; // inner (seed sensor)										
													if (thisDetId >= 0 && thisDetId < static_cast<int>(stackMapSize)) {
														out[outIdx].detId() = stackMap[thisDetId].first;
													}
													#ifdef Debug_GPU
													printf("Unpacked values P on PS: chipID = %u, addr = %u, size = %u, z = %u, strip = %u, row = %u, column = %u, edge = %u\n",
														   chip, addr, w, z, out[outIdx].strip(), out[outIdx].row(), out[outIdx].column(), out[outIdx].edge());
										#endif
												}
												localClusterIdx++;
											}
										}
										alpaka::syncBlockThreads(acc);
								 } // end loop on channels for this dtc
							 //} // end fed data size > 0
						 } // independatn group elements 
						 alpaka::syncBlockThreads(acc);
					 } // call operator  
	};


	// Kernel for 2S modules: unpack only stripClustersWords into (x,y,width) TODO :: chabnge the unpackers from kernal to functions (remove the operator and like the functions on top )
	// 2. move uniform elements outside the function now 
	// 3. in the previous todo we added the reserved places before the unpacking and this needs to write the output to the reserved the places : 
	// Launch the generic Unpacker kernel (header + payload) on device
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
			uint32_t* globalCounter) {
		const uint32_t threadsPerBlock = 128;
		// +1 added for normalization of the 3D indexing
		//const uint32_t blocks = (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC ;
		// removing +1
		const uint32_t blocks = (MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC ;
		//Adjust the work division to account for channel-level parallelism:
/*		const uint32_t threadsPerBlock = 128;
		const uint32_t blocks = ((MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC * CICs_PER_SLINK + threadsPerBlock - 1) / threadsPerBlock;
 */
		auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

		alpaka::exec<Acc1D>(
				queue,
				workDiv,
				Unpacker{},
				rawdatabuff.data(),
				sizedatabuff.data(),
				offsetdatabuff.data(),
				inmap.data(),
				detIdMap.data(),
				stackMap.data(), // add stackMap
				static_cast<uint32_t>(stackMapSize),  // Pass the size
				out, 
				globalCounter
				);
	}
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// Specialize trait to tell Alpaka how much to allocate 
// Specialization of BlockSharedMemDynSizeBytes for Unpacker kernel
// This specialization tells Alpaka how much dynamic shared memory to allocate for the Unpacker kernel
// This is needed to avoid illegal memory access errors
// The size is the total number of 32-bit words needed, multiplied by sizeof(uint32_t) to get bytes
namespace alpaka::trait {
	template<>
		struct BlockSharedMemDynSizeBytes<Unpacker, Acc1D> {
			template<
				typename RawBufView,
					 typename SizeBufView,
					 typename OffBufView,
					 typename InMapView,
					 typename DetIdMapView, // Add DetIdMapView
    				 typename StackMapView, // Add StackMapView
					 typename OutView
						 >
						 ALPAKA_FN_HOST_ACC static std::size_t
						 getBlockSharedMemDynSizeBytes(
								 Unpacker const & /*kernel*/,
								 Vec1D threads, 
								 Vec1D elements,
								 RawBufView, //const* /*in*/,
								 SizeBufView, //const* /*sizes*/,
								 OffBufView, //const* /*offsets*/,
								 DetIdMapView, // Add detIdMap
    							 StackMapView, // Add stackMap
								 InMapView, //const* /*detIdxMap*/
								 uint32_t stackMapSize,  // Add stackMapSize parameter
								 OutView,
								 uint32_t* globalCounter
								 ) {
							 return static_cast<std::size_t>(MaxTotalSharedWords) * sizeof(uint32_t);
						 }
		};
} // namespace alpaka::trait
