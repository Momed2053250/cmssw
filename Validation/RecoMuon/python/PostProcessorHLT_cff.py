# post-processors for HLT Muon track validation in FullSim and FastSim
#
import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

from Validation.RecoMuon.PostProcessor_cff import postProcessorMuonTrack, postProcessorMuonTrackSummary

postProcessorMuonTrackHLT = postProcessorMuonTrack.clone(
    subDirs = ["HLT/Muon/MuonTrack/*"]
)

postProcessorMuonTrackHLTSummary = postProcessorMuonTrackSummary.clone(
    subDirs=cms.untracked.vstring("HLT/Muon/MuonTrack/")
)

postProcessorMuonTrackHLTComp = DQMEDHarvester("DQMGenericClient",
    subDirs=cms.untracked.vstring("HLT/Muon/MuonTrack/"),
    efficiency=cms.vstring(
        "Eff_L3Tk_Eta_mabh 'Eff_{L3,TK} vs #eta' hltL3Muons/effic_vs_eta hltL3TkFromL2/effic_vs_eta",
        "Eff_L3Tk_Pt_mabh 'Eff_{L3,TK} vs p_{T}' hltL3Muons/effic_vs_pt hltL3TkFromL2/effic_vs_pt",
        "Eff_L3Tk_Hit_mabh 'Eff_{L3,TK} vs n Hits' hltL3Muons/effic_vs_hit hltL3TkFromL2/effic_vs_hit",
        "Eff_L3L2_Eta_mabh 'Eff_{L3,L2} vs #eta' hltL3Muons/effic_vs_eta hltL2Muons_UpdAtVtx/effic_vs_eta",
        "Eff_L3L2_Pt_mabh 'Eff_{L3,L2} vs p_{T}' hltL3Muons/effic_vs_pt hltL2Muons_UpdAtVtx/effic_vs_pt",
        "Eff_L3L2_Hit_mabh 'Eff_{L3,L2} vs n Hits' hltL3Muons/effic_vs_hit hltL2Muons_UpdAtVtx/effic_vs_hit",
    ),
    resolution=cms.vstring(""),
    outputFileName=cms.untracked.string(""),
)

recoMuonPostProcessorsHLT = cms.Sequence(
    postProcessorMuonTrackHLT + postProcessorMuonTrackHLTSummary + postProcessorMuonTrackHLTComp 
)

from Configuration.Eras.Modifier_phase2_muon_cff import phase2_muon

phase2_muon.toReplaceWith(
    recoMuonPostProcessorsHLT,
    recoMuonPostProcessorsHLT.copyAndExclude(["postProcessorMuonTrackHLTComp"]),
)
