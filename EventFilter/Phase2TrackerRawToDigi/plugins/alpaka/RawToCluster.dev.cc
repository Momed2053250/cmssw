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
	ALPAKA_FN_ACC uint16_t getOffsetForChannel (unsigned int iChannel, uint32_t* offsetWords){ 
		//TODO:: Optimize
		// Even channel: lower 16 bits of word iChannel/2
		if (iChannel % 2 == 0) {
			// extract the lower 16 bits by masking with 0xFFFF
			return static_cast<uint16_t>(offsetWords[iChannel/2] & 0xFFFF);
		}
		else {
			// Odd channel: upper 16 bits of word (iChannel-1)/2
			// extract the upper 16 bits by shifting right by 16
			return static_cast<uint16_t>(offsetWords[(iChannel -1)/2] >> 16) ;
		}	
	}

	// writing the kernals 
	//
	// Unpacker kernel: top-level device loop over FED fragments
	struct Unpacker {
		template <
			typename Acc,
				 typename RawBufView,
				 typename SizeBufView,
				 typename OffBufView,
				 typename InMapView
					 >
					 ALPAKA_FN_ACC void operator()(
							 Acc const& acc,
							 RawBufView         in,    // device_buffer<Device,unsigned char[]>
							 SizeBufView        sizes,   // device_buffer<Device,size_t[]>
							 OffBufView         offsets, // device_buffer<Device,size_t[]>
							 InMapView          const& detIdxModuleTypeMap //,          // device_buffer<Device,int[]>
							 ) const {
						 // (A) Allocate ONE contiguous chunk of dynamic shared memory:
						 uint8_t*  smemBytes = alpaka::getDynSharedMem<uint8_t>(acc);
						 uint32_t*	smemWords = reinterpret_cast<uint32_t*>(smemBytes);
						 if (!smemBytes) {
							 printf("Error: smemBytes is null\n");
							 return;
						 }

						 // (B) Slice that chunk into disjoint regions:
						 uint32_t* headerWords       = smemWords;                           // [0 .. MaxHeaderWords-1]
						 uint32_t* offsetWords       = headerWords       + MaxHeaderWords;  // [MaxHeaderWords .. ]
						 uint32_t* lines             = offsetWords       + MaxOffsetWords;  // [ ... ]
						 uint32_t* stripClusterWords = lines             + MaxPayloadLines;  // [ ... ]
						 uint32_t* pixelClusterWords = stripClusterWords + MaxStripClusters; // [ ... ]


						 // Iterate over each FED fragment ID in parallel with indepndant group tool 
						 for (auto frdId : cms::alpakatools::independent_groups(acc,
									 (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC)) {

							 // process only if fedData.size() > 0
							 if (sizes[frdId] > 0) 
							 {

								 // Pointer to the start of fragment data
								 //const unsigned char* dataPtr = in;
								 const unsigned char* dataPtr = in + offsets[frdId];

								 // 1) Read the header (HEADER_N_LINES words) into headerWords[]:
								 size_t nHeaderLines =  HEADER_N_LINES;
								 for (auto k : cms::alpakatools::independent_group_elements(acc, nHeaderLines)) {
									 // int byteIdx = static_cast<int>(k) * N_BYTES_PER_WORD;
									 auto byteIdx = k * N_BYTES_PER_WORD;
									 //    printf("Before %u \n", k);
									 if (k >= MaxHeaderWords) {
										 printf("Error: k=%u exceeds MaxHeaderWords=%d\n", k, MaxHeaderWords);
										 return;
									 }
									 headerWords[k] = readLine(dataPtr, byteIdx);
								 }
								 alpaka::syncBlockThreads(acc);
								 // read the offsets into shared memory: each 32 bit word contains two offset words of 16 bit each
								 // 2) Read offset words into offsetWords[]:
								 size_t nOffsetsLines = (OFFSET_BITS * CICs_PER_SLINK) / N_BITS_PER_WORD;
								 size_t initByte      = HEADER_N_LINES * N_BYTES_PER_WORD;
								 for (auto k : cms::alpakatools::independent_group_elements(acc, nOffsetsLines)) {
									 int byteIdx = static_cast<int>(initByte + k * N_BYTES_PER_WORD);
									 offsetWords[k] = readLine(dataPtr, byteIdx);
								 }
								 alpaka::syncBlockThreads(acc);
								 printf("offsetWords[0] = %u\n", offsetWords[0]);


								 // now read the payload (channel header + clusters)
								 // all channel headers should be there, even if 0 clusters are found
								 // the loop is not on the actual channel number, as in the ClusterToRaw conversion each channel is split by CIC0_CIC1
								 // in order to get all the clusters from the same lpGBT and fill them once at the end
								 //---------- Unpack each channel ----------

								 // Debug print out counter 
								 static int printoutCount = 0;
								 for (unsigned int iChannel = 0; iChannel < CICs_PER_SLINK; iChannel++)
								 {
									 // clear the collection if iChannel is even   // See if this need to be a TODO
									 if (iChannel%2==0){
										 //     thisChannel1DSeedClusters.clear();
										 //     thisChannel1DCorrClusters.clear();
									 }  

									 // retrieve the module type:  	TODO: 	check the logic if debug 
									 //	int is2SModule = detIdxModuleTypeMap[frdId+iChannel];
									 // first, compute your flat index into the [dtc][slink][channel] array:
									 //
									 const unsigned CICs = CICs_PER_SLINK;
									 unsigned flatIdx = frdId * CICs + iChannel;
									 // then use 
									 int is2SModule = detIdxModuleTypeMap[ flatIdx ];

									 // Compute byteindex idx of this channels header:
									 size_t offsetTableStart = HEADER_N_LINES * N_BYTES_PER_WORD;
									 //  must match exactly where your offset table begins 
									 int channelOffset16 = static_cast<int>( getOffsetForChannel(iChannel, offsetWords) );
									 int idx = static_cast<int>(offsetTableStart + channelOffset16 * N_BYTES_PER_WORD);
									 // ensure that channelOffset16 times 4 truly points to the channel header 

									 // Read channel header word extract cluster counts:
									 uint32_t chHeaderWord = readLine(dataPtr, idx);
									 unsigned int numStripClusters =
										 (chHeaderWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS - N_STRIP_CLUSTER_BITS))
										 & N_CLUSTER_MASK;  // 7 bits
									 unsigned int numPixelClusters = chHeaderWord & N_CLUSTER_MASK;  // 7 bits

									 /*						// find the channel offset
															int initial_offset = (HEADER_N_LINES + MODULES_PER_SLINK) * N_BYTES_PER_WORD;
															int idx = initial_offset + getOffsetForChannel(iChannel, offsetWords) * N_BYTES_PER_WORD;

									 // get the channel header and unpack it to get the cluster count 
									 uint32_t headerWord = readLine(dataPtr, idx);
									 unsigned int numStripClusters = (headerWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS - N_STRIP_CLUSTER_BITS)) & N_CLUSTER_MASK; // 7-bit field
									 unsigned int numPixelClusters = (headerWord) & N_CLUSTER_MASK; // 7-bit field
									 */

									 // define the number of lines of the payload
									 unsigned int nLines = (numStripClusters + numPixelClusters > 0) ? 
										 int((numStripClusters * SS_CLUSTER_BITS + numPixelClusters * PX_CLUSTER_BITS)/ N_BITS_PER_WORD) + 1 : 0;
									 // Print cluster counts with limit of 10
									 if (printoutCount < 100 && numStripClusters > 0) {
										 printf("n strip clusters are: %u\n", numStripClusters);
										 printoutCount++;
									 }
									 if (printoutCount < 100 && numPixelClusters > 0) {
										 printf("n pixel clusters are: %u\n", numPixelClusters);
										 printoutCount++;
									 }
									 /*          
										     if (numStripClusters + numPixelClusters > 0 ){
										     LogTrace("RawToClusterProducer") << "\tchannel " << iChannel << "\theader: " << std::bitset<N_BITS_PER_WORD>(headerWord) 
										     << "\tn strip clusters = " << numStripClusters 
										     << "\tn pixel clusters = " << numPixelClusters 
										     << " (n lines = " << nLines << ")";
										     }  
										     */         
									 // first retrieve all lines filled with clusters
									 for (auto k : cms::alpakatools::independent_group_elements(acc, nLines)) {
										 int byteIdx = idx + static_cast<int>(k) * N_BYTES_PER_WORD;
										 lines[k] = readLine(dataPtr, byteIdx);
									 }
									 alpaka::syncBlockThreads(acc);

									 /*   if ( lines.size() != nLines) {
									      edm::LogError("RawtoClusterProducer") << "Numbers of stored lines does not match with size of lines to be read!";
									      return;
									      }  
									      */


									 // create groups of 14 (17) bits for 2S (PS) clusters, joining consecutive lines if needed
									 int nAvailableBits = N_BITS_PER_WORD;
									 int iLine = 0;
									 int bitsToRead = 0;
									 int nFullClustersStrip = 0;
									 int nFullClustersPix = 0;
									 //  make sure nLines <= MaxPayloadLines
									 for (auto k : independent_group_elements(acc, nLines)) {
										 int byteIdx = getLineIndex(idx, static_cast<unsigned int>(k));
										 lines[k] = readLine(dataPtr, byteIdx);
									 }
									 alpaka::syncBlockThreads(acc);
									 // Print out for intermediate check
									 // Print lines[0] with limit of 10
									 if (printoutCount < 10 && nLines > 0) {
										 printf("lines[0] = %u\n", lines[0]);
										 printoutCount++;
									 }
									 // Read payloads for 2S 
									 if (is2SModule) {

										 // Consider only one thread per thread block with once_per_block 
										 if (cms::alpakatools::once_per_block(acc)){
											 readPayload(stripClusterWords, lines, numStripClusters, nAvailableBits, iLine, bitsToRead, nFullClustersStrip, SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false );

										 }
									 } 
									 // For the PS models 
									 else {

										 if (cms::alpakatools::once_per_block(acc)){
											 readPayload(stripClusterWords, lines, numStripClusters, nAvailableBits, iLine, bitsToRead, nFullClustersStrip, SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false );
											 readPayload(pixelClusterWords, lines, numPixelClusters, nAvailableBits, iLine, bitsToRead, nFullClustersPix, PX_CLUSTER_BITS, PX_CLUSTER_WORD_MASK, true, nFullClustersStrip );
										 }
									 }

									 // TODOs: (after unpacked)
									 // 1. Need to define an output SoA with two columns for seed clus and corr clusters with the max size numclus * channels * dtcs * slinks = 128 *36 *~800 * ~4 
									 // 2. Need to define a global atomic which will decide where a block write out the results in the output soA 
									 // 3. Read payload is done and after reading just use these dynmem and send them to unpackers   
									 /*
									 // unpack the cluster words and create Phase2TrackerCluster1D objects
									 int count_clusters = 0;
									 if (is2SModule) {
									 // create the Phase2TrackerCluster1D objects for 2S modules
									 for (auto icluster : stripClustersWords){
									 std::pair<Phase2TrackerCluster1D, bool> thisCluster = unpack2S(icluster, iChannel);
									 if (thisCluster.second)
									 thisChannel1DSeedClusters.push_back(thisCluster.first);
									 else  
									 thisChannel1DCorrClusters.push_back(thisCluster.first);
									 count_clusters++;
									 } // end loop on cluster words 
									 } else {
									 // create the Phase2TrackerCluster1D objects for PS modules
									 // first loop on strip clusters
									 for (auto icluster : stripClustersWords){
									 Phase2TrackerCluster1D thisCluster = unpackStripOnPS(icluster, iChannel);
									 // for PS, strip is always correlated sensor
									 thisChannel1DCorrClusters.push_back(thisCluster);
									 count_clusters++; 
									 } 
									 // then loop on pixel clusters
									 for (auto icluster : pixelClustersWords){
									 Phase2TrackerCluster1D thisCluster = unpackPixelOnPS(icluster, iChannel);
									 // for PS, pixel is always seed sensor
									 thisChannel1DSeedClusters.push_back(thisCluster);
									 count_clusters++;
									 } 
									 }

									 // use FastFiller to fill the output DetSetVector output collection
									 // fill every time that 2 channels are read
									 if (iChannel%2 != 1) 
									 continue;

									 // Store clusters of this channel
									 std::vector<Phase2TrackerCluster1D>::iterator it;
									 {
									 // inner detid is defined as module detid + 1. First int in the pair from the map
									 edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller spcs(*outputClusterCollection, stackMap_[thisDetId].first);
									 for (it = thisChannel1DSeedClusters.begin(); it != thisChannel1DSeedClusters.end(); it++) {
									 spcs.push_back(*it);
									 }
									 }
									 {
									 // outer detid is defined as inner detid + 1 or module detid + 2. Second int in the pair from the map
									 edmNew::DetSetVector<Phase2TrackerCluster1D>::FastFiller spcc(*outputClusterCollection, stackMap_[thisDetId].second);
									 for (it = thisChannel1DCorrClusters.begin(); it != thisChannel1DCorrClusters.end(); it++) {
									 spcc.push_back(*it);
									 }
									 }
									 */ 
								 } // end loop on channels for this dtc

							 } // end fed data size > 0

						 } // independatn group elements 
						 alpaka::syncBlockThreads(acc);
					 } // call operator  
	};


	// Kernel for 2S modules: unpack only stripClustersWords into (x,y,width) TODO :: chabnge the unpackers from kernal to functions (remove the operator and like the functions on top )
	// 2. move uniform elements outside the fucniton now 
	// 3. in the previous todo we added the reserved places before the unpacking and this needs to write the output to the reserved the places : 
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

	// Launch the generic Unpacker kernel (header + payload) on device
	void launchUnpacker(
			Queue& queue,
			cms::alpakatools::device_buffer<Device, unsigned char[]> rawdatabuff,
			cms::alpakatools::device_buffer<Device, size_t[]> sizedatabuff,
			cms::alpakatools::device_buffer<Device, size_t[]> offsetdatabuff,
			cms::alpakatools::device_buffer<Device, int[]> inmap) {
		const uint32_t threadsPerBlock = 128;
		// +1 added for normalization of the 3D indexing
		const uint32_t blocks = (MAX_DTC_ID - MIN_DTC_ID + 1) * SLINKS_PER_DTC ;
		auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);

		alpaka::exec<Acc1D>(
				queue,
				workDiv,
				Unpacker{},
				rawdatabuff.data(),
				sizedatabuff.data(),
				offsetdatabuff.data(),
				inmap.data()//,
				);
	}
	/*
	// Launch the 2S unpacker kernel on device
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
	*/
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
					 typename InMapView
						 >
						 ALPAKA_FN_HOST_ACC static std::size_t
						 getBlockSharedMemDynSizeBytes(
								 Unpacker const & /*kernel*/,
								 Vec1D threads, 
								 Vec1D elements,
								 RawBufView, //const* /*in*/,
								 SizeBufView, //const* /*sizes*/,
								 OffBufView, //const* /*offsets*/,
								 InMapView //const* /*detIdxMap*/
								 ) {
							 return static_cast<std::size_t>(MaxTotalSharedWords) * sizeof(uint32_t);
						 }
		};
} // namespace alpaka::trait


