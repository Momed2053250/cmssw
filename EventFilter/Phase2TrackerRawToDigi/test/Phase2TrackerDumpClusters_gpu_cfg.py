import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("DumpClustersSoA")

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

options = VarParsing.VarParsing('analysis')
options.register("onUnpacked",
                False,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.bool,
                "define if running on unpacked clusters or on original ones")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.Phase2TrackerDumpClustersSoA = cms.EDAnalyzer(
    'Phase2TrackerDumpClustersSoA',
    ProductLabel=cms.InputTag("phase2RawToCluster", "ClusterPropSoA")
)

process.TFileService = cms.Service('TFileService',
    fileName=cms.string(options.outputFile),
    closeFileFast=cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles)
)

if options.onUnpacked:
    process.Phase2TrackerDumpClustersSoA.ProductLabel = cms.InputTag("Unpacker", "", "PACKANDUNPACK")

process.load('Configuration.Geometry.GeometryExtendedRun4D98Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun4_realistic_v1', '')

process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'frontier://FrontierProd/CMS_CONDITIONS'
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    DumpStat=cms.untracked.bool(True),
    toGet=cms.VPSet(cms.PSet(
        record=cms.string('TrackerDetToDTCELinkCablingMapRcd'),
        tag=cms.string("TrackerDetToDTCELinkCablingMap__OT800_IT711__T33__OTOnly"),
    )),
)
process.es_prefer_local_cabling = cms.ESPrefer("PoolDBESSource", "")

process.p = cms.EndPath(process.Phase2TrackerDumpClustersSoA)