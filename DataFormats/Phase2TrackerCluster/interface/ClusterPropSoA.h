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
                      // columns: one value per element
                      SOA_COLUMN(uint32_t, width),
                      SOA_COLUMN(unsigned int, x),
		      SOA_COLUMN(unsigned int, y),
		      SOA_COLUMN(bool , seedFlag)
  )
  using ClusterPropSoA = ClusterPropSoALayout<>;
 
}  // namespace Phase2RawToCluster

#endif  // DataFormats_Phase2TrackerCluster_interface_ClusterPropSoA_h
