import FWCore.ParameterSet.Config as cms
# ---------- Add assigned jet-track association

from RecoJets.JetAssociationProducers.ak4JTA_cff import *
ak4JetTracksAssociatorAtVertexJPT = ak4JetTracksAssociatorAtVertex.clone(
    useAssigned = True,
    pvSrc       = "offlinePrimaryVertices"
)

# ---------- Tight Electron ID

from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff import egmGsfElectronIDs
JPTegmGsfElectronIDs = egmGsfElectronIDs.clone(
    physicsObjectsIDs = cms.VPSet(),
    physicsObjectSrc = 'gedGsfElectrons'
)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupVIDSelection
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff']
for id_module_name in my_id_modules:
    idmod= __import__(id_module_name, globals(), locals(), ['idName','cutFlow'])
    for name in dir(idmod):
        item = getattr(idmod,name)
        if hasattr(item,'idName') and hasattr(item,'cutFlow'):
            setupVIDSelection(JPTegmGsfElectronIDs,item)

# ---------- Seeds from TrackJets

from RecoJets.JetPlusTracks.jetPlusTrackAddonSeedProducer_cfi import *

JetPlusTrackAddonSeedReco = jetPlusTrackAddonSeedProducer.clone()

# ---------- Module definition
from RecoJets.JetPlusTracks.JetPlusTrackCorrections_cfi import *

JetPlusTrackZSPCorJetAntiKt4 = cms.EDProducer(
    "JetPlusTrackProducer",
    cms.PSet(JPTZSPCorrectorAntiKt4),
    src = cms.InputTag("ak4CaloJets"),
    srcTrackJets = cms.InputTag("ak4TrackJets"),
    srcAddCaloJets = cms.InputTag('JetPlusTrackAddonSeedReco'),
    extrapolations = cms.InputTag("trackExtrapolator"),
    tagName = cms.vstring('ZSP_CMSSW390_Akt_05_PU0'),
    tagNameOffset = cms.vstring(),
    PU = cms.int32(-1),
    FixedPU = cms.int32(0),
    UseZSP = cms.bool(False),
    srcPVs = cms.InputTag('offlinePrimaryVertices'),    
    alias = cms.untracked.string('JetPlusTrackZSPCorJetAntiKt4'),
    ptCUT = cms.double(15.),
    dRcone = cms.double(0.4)
    )

JetPlusTrackZSPCorJetAntiKt4.JetTracksAssociationAtVertex   = "ak4JetTracksAssociatorAtVertexJPT"
JetPlusTrackZSPCorJetAntiKt4.JetTracksAssociationAtCaloFace = "ak4JetTracksAssociatorAtCaloFace"
JetPlusTrackZSPCorJetAntiKt4.JetSplitMerge = 2

### ---------- Sequences

# Anti-Kt

JetPlusTrackCorrectionsAntiKt4Task = cms.Task(
    JetPlusTrackAddonSeedReco,
    ak4JetTracksAssociatorAtVertexJPT,
    ak4JetTracksAssociatorAtCaloFace,
    JetPlusTrackZSPCorJetAntiKt4
    )
JetPlusTrackCorrectionsAntiKt4 = cms.Sequence(JetPlusTrackCorrectionsAntiKt4Task)

# For backward-compatiblity (but to be deprecated!)
JetPlusTrackCorrections = cms.Sequence(JetPlusTrackCorrectionsAntiKt4)
