// Host-only converter: SoA -> edmNew::DetSetVector<Phase2TrackerCluster1D>
// added: host-only version that consumes the Host SoA and writes legacy AoS

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"

// added: host SoA only, no alpaka device types here
#include "DataFormats/Phase2TrackerCluster/interface/ClusterPropHostCollection.h"

#include <map>
#include <vector>
#include <algorithm>
#include <cstdint>

class ClusterPropSoAToLegacyED : public edm::stream::EDProducer<> {
public:
  explicit ClusterPropSoAToLegacyED(edm::ParameterSet const& iConfig)
      : soaToken_{consumes<Phase2RawToCluster::ClusterPropHostCollection>(
            iConfig.getParameter<edm::InputTag>("clusterSoASource"))},
        legacyOutToken_{produces<Phase2TrackerCluster1DCollectionNew>()} {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("clusterSoASource", edm::InputTag("Unpacker")); // added
    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::Event& iEvent, edm::EventSetup const&) override {
    auto const& host = iEvent.get(soaToken_);
    auto const view = host.view();

    // Gather by detId, then fill DetSetVector via FastFiller
    std::map<uint32_t, std::vector<Phase2TrackerCluster1D>> perDet;

    const int n = view.metadata().size();

    for (int i = 0; i < n; ++i) {
      const uint32_t detId  = view[i].detId();
      const uint32_t row    = view[i].row();
      const uint32_t column = view[i].column();
      const uint32_t size   = view[i].size();
      const bool     thr    = view[i].threshold();

      if (thr) {
        perDet[detId].emplace_back(row, column, size, 1u);
      } else {
        perDet[detId].emplace_back(row, column, size);
      }
    }

    // Sort clusters within each det by first strip (legacy expectation)
    for (auto& kv : perDet) {
      auto& v = kv.second;
      std::sort(v.begin(), v.end());
    }

    // Build the output with FastFiller
    auto out = std::make_unique<Phase2TrackerCluster1DCollectionNew>();
    out->reserve(perDet.size(), n);

    for (auto& [rawDetId, vec] : perDet) {
      DetId det(rawDetId);
      Phase2TrackerCluster1DCollectionNew::FastFiller filler(*out, det);
      filler.reserve(vec.size()); // optional
      for (auto const& c : vec) {
        filler.push_back(c);
      }
    }

    iEvent.put(legacyOutToken_, std::move(out));
  }

  const edm::EDGetTokenT<Phase2RawToCluster::ClusterPropHostCollection> soaToken_;
  const edm::EDPutTokenT<Phase2TrackerCluster1DCollectionNew> legacyOutToken_;
};

DEFINE_FWK_MODULE(ClusterPropSoAToLegacyED);
