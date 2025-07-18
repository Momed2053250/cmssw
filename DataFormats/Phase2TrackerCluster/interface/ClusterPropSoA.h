#ifndef DataFormats_Phase2TrackerCluster_interface_ClusterPropSoA_h
#define DataFormats_Phase2TrackerCluster_interface_ClusterPropSoA_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace Phase2RawToCluster {

  // Generate structure of arrays (SoA) layout
  GENERATE_SOA_LAYOUT(ClusterPropSoALayout,
                      SOA_COLUMN(uint32_t, size),        
                      SOA_COLUMN(uint32_t, strip),       
                      SOA_COLUMN(uint32_t, row),         
                      SOA_COLUMN(uint32_t, column),      
                      SOA_COLUMN(uint32_t, edge),        
                      SOA_COLUMN(bool, threshold),       
                      SOA_COLUMN(uint32_t, mipBit)       
  )
  using ClusterPropSoA = ClusterPropSoALayout<>;

}  // namespace Phase2RawToCluster

#endif  // DataFormats_Phase2TrackerCluster_interface_ClusterPropSoA_h