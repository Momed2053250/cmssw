// alpaka-related imports
#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/FEDRawData/interface/alpaka/StripPixelDeviceCollection.h"
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropDeviceCollection.h"
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

// writing the kernals 
  struct Unpacker {
          template <typename Acc, typename InMap, typename RawDataBuff>
          ALPAKA_FN_ACC void operator()(Acc const& acc,
                                        RawDataBuff in, InMap const& detIdxModuleTypeMap
                                        ) const {
		  for (auto frdId : cms::alpakatools::independent_groups(acc, (MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC)){
				  
				  


		  			 if (in[frdId].size() > 0) 
        {/*
          const unsigned char* dataPtr = fedData.data();
          
          // read the header
          std::vector<uint32_t> headerWords;
          for (size_t i = 0; i < HEADER_N_LINES*N_BYTES_PER_WORD; i += N_BYTES_PER_WORD) // Read 4 bytes (32 bits) at a time
          {
            // Extract 4 bytes (32 bits) and pack them into a uint32_t word
            headerWords.push_back(readLine(dataPtr, i));
          }
          theHeader.setValue(headerWords);
    
          // read the offsets: each 32 bit word contains two offset words of 16 bit each
          std::vector<uint32_t> offsetWords;
          size_t nOffsetsLines = OFFSET_BITS * CICs_PER_SLINK / N_BITS_PER_WORD;
          size_t initByte = HEADER_N_LINES*N_BYTES_PER_WORD;
          size_t endByte = (nOffsetsLines-1)*N_BYTES_PER_WORD + initByte; // -1 because we only need the starting i of the line
  
          for (size_t i = initByte; i <= endByte; i += N_BYTES_PER_WORD) // Read 4 bytes (32 bits) at a time
            offsetWords.push_back(readLine(dataPtr, i));
          theOffsets.setValue(offsetWords);
          
          // now read the payload (channel header + clusters)
          // all channel headers should be there, even if 0 clusters are found
          // the loop is not on the actual channel number, as in the ClusterToRaw conversion each channel is split by CIC0_CIC1
          // NOTE: we need to save into the Phase2TrackerCluster1D collection two "channels" at the time 
          // in order to get all the clusters from the same lpGBT and fill them once at the end
          std::vector<Phase2TrackerCluster1D> thisChannel1DSeedClusters, thisChannel1DCorrClusters;
          for (unsigned int iChannel = 0; iChannel < CICs_PER_SLINK; iChannel++)
          {
            // clear the collection if iChannel is even
            if (iChannel%2==0){
              thisChannel1DSeedClusters.clear();
              thisChannel1DCorrClusters.clear();
            }  
  
            // retrieve the module type: 
            // first we need to construct the DTCElinkId object ## dtc_id, gbtlink_id, elink_id
            // to get the gbt_id we should reverse what is done in the packer function,
            // where clusters from channel X are split into 2*i and 2*i+1 based on being from CIC0 or CIC1
            unsigned int gbt_id = iSlink * MODULES_PER_SLINK + std::div(iChannel, 2).quot;
            DTCELinkId thisDTCElinkId(dtcID, gbt_id, 0);
  
            int thisDetId = -1;
            bool is2SModule = false;
            // then pass it to the map to get the detid
            if (cablingMap_->knowsDTCELinkId(thisDTCElinkId))
            {
              auto possibleDetIds = cablingMap_->dtcELinkIdToDetId(thisDTCElinkId); // this returns a pair, detid will be an uint32_t (not a DetId)
              thisDetId = possibleDetIds->second;
              LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID) << "\tiGBT: " << unsigned(gbt_id) 
                                               << "\tielink: " << unsigned(0) << "\t -> detId:" << thisDetId; 
              // check is 2S or PS 
              is2SModule = trackerGeometry_->getDetectorType( stackMap_[thisDetId].first) == TrackerGeometry::ModuleType::Ph2SS;
            }
            else {
              LogTrace("RawToClusterProducer") << "slink: " << iSlink << "\tiDTC: " << unsigned(dtcID) << "\tiGBT: " << unsigned(gbt_id) 
                                               << " -> not connected? " ; 
              continue;
            }
  
            // find the channel offset
            int initial_offset = (HEADER_N_LINES + MODULES_PER_SLINK) * N_BYTES_PER_WORD;
            int idx = initial_offset + theOffsets.getOffsetForChannel(iChannel) * N_BYTES_PER_WORD;
            
            // get the channel header and unpack it
            uint32_t headerWord = readLine(dataPtr, idx);
            // unsigned long eventID = (headerWord >> (N_BITS_PER_WORD - L1ID_BITS)) & L1ID_MAX_VALUE; // 9-bit field
            // int channelErrors = (headerWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS)) & CIC_ERROR_MASK; // 9-bit field
            unsigned int numStripClusters = (headerWord >> (N_BITS_PER_WORD - L1ID_BITS - CIC_ERROR_BITS - N_STRIP_CLUSTER_BITS)) & N_CLUSTER_MASK; // 7-bit field
            unsigned int numPixelClusters = (headerWord) & N_CLUSTER_MASK; // 7-bit field
            
            // define the number of lines of the payload
            unsigned int nLines = (numStripClusters + numPixelClusters > 0) ? 
                                   int((numStripClusters * SS_CLUSTER_BITS + numPixelClusters * PX_CLUSTER_BITS)/ N_BITS_PER_WORD) + 1 : 0;
  
            if (numStripClusters + numPixelClusters > 0 ){
              LogTrace("RawToClusterProducer") << "\tchannel " << iChannel << "\theader: " << std::bitset<N_BITS_PER_WORD>(headerWord) 
                                               << "\tn strip clusters = " << numStripClusters 
                                               << "\tn pixel clusters = " << numPixelClusters 
                                               << " (n lines = " << nLines << ")";
            }  
            
            // first retrieve all lines filled with clusters
            std::vector<uint32_t> lines;
            for (unsigned int iline = 0; iline < nLines; iline++) {
              lines.push_back(readLine(dataPtr, getLineIndex(idx, iline)));
            }        
            if ( lines.size() != nLines) {
              edm::LogError("RawtoClusterProducer") << "Numbers of stored lines does not match with size of lines to be read!";
              return;
            }  
    
            // first retrieve the cluster words
            // this was uint16, check if can be changed back 
            std::vector<uint32_t> stripClustersWords;
            stripClustersWords.resize(numStripClusters);
  
            std::vector<uint32_t> pixelClustersWords;
            pixelClustersWords.resize(numPixelClusters);
      
            // create groups of 14 (17) bits for 2S (PS) clusters, joining consecutive lines if needed
            int nAvailableBits = N_BITS_PER_WORD;
            int iLine = 0;
            int bitsToRead = 0;
            int nFullClustersStrip = 0;
            int nFullClustersPix = 0;
            
            readPayload(stripClustersWords, lines, numStripClusters, nAvailableBits, iLine, bitsToRead, nFullClustersStrip, SS_CLUSTER_BITS, SS_CLUSTER_WORD_MASK, false );
            readPayload(pixelClustersWords, lines, numPixelClusters, nAvailableBits, iLine, bitsToRead, nFullClustersPix, PX_CLUSTER_BITS, PX_CLUSTER_WORD_MASK, true, nFullClustersStrip );
            
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
    
          } // end loop on channels for this dtc
        */
	} // end fed data size > 0

	  } // independatn group elements 
	  } // call operator  
  };





  // Kernel for 2S modules: unpack only stripClustersWords into (x,y,width)
	struct S2UnpackerKernel {
	  template <typename Acc, typename InView, typename OutView>
	  ALPAKA_FN_ACC void operator()(Acc const& acc,
	                                InView const& in,
	                                OutView out) const {
			// what is this metadata class check if the same for OT or this way is not correct 
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

      // unpacker kernal launch
      template <typename InMap, typename RawDataBuff> 
       void launchUnpacker(
	Queue& queue,
	RawDataBuff rawdatabuff,
	InMap inmap) {
	       const uint32_t threadsPerBlock = 128;
	       const uint32_t blocks = (MAX_DTC_ID - MIN_DTC_ID) * SLINKS_PER_DTC ;
	       auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, threadsPerBlock);
          alpaka::exec<Acc1D>(queue, workDiv, Unpacker{}, rawdatabuff, inmap);
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


