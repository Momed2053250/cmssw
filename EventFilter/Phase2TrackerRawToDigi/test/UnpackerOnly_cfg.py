## cfg file to run ONLY the unpacking and converting steps for Phase2 OT clusters
## identical services/accelerators/tracers/conditions; path runs just Unpacker and Converter

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.Utilities.FileUtils as FileUtils
import os

# flag for using the conversion
Legacy_Format = True

process = cms.Process("PACKANDUNPACK")
process.options.numberOfThreads = 8
process.options.numberOfStreams = 8

def get_input_mc_line(dataset_database, line_number):
    with open(dataset_database, 'r') as file:
        lines = file.readlines()
        if line_number < 0 or line_number >= len(lines):
            raise IndexError("Line number out of range")
        return lines[line_number].strip()

options = VarParsing.VarParsing('analysis')
options.register('cluster',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Cluster ID from HTCondor")
options.register('process',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Process ID from HTCondor")
options.parseArguments()

GEOMETRY = "D98"

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Accelerators_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

if GEOMETRY == "D88" or GEOMETRY == 'D98':
    process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY + 'Reco_cff')
    process.load('Configuration.Geometry.GeometryExtendedRun4' + GEOMETRY + '_cff')
else:
    print("this is not a valid geometry!!!")

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun4_realistic_v1', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# ---- IMPORTANT: read the RAW produced by the pack step ----
# If you prefer to pass via CLI, you can: cmsRun Phase2_unpack_convert_only_cfg.py inputFiles=file:raw2clusters.root
# Otherwise this default points to the pack step’s output:
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:raw2clusters.root"
    )
)

process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.connect = 'frontier://FrontierProd/CMS_CONDITIONS'

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDB,
    DumpStat = cms.untracked.bool(True),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('TrackerDetToDTCELinkCablingMapRcd'),
        tag = cms.string("TrackerDetToDTCELinkCablingMap__OT800_IT711__T33__OTOnly"),
    )),
)

process.es_prefer_local_cabling = cms.ESPrefer("PoolDBESSource", "")

# No Clusterizer and no Packer in this file.
# Only consume the FEDRawData from the previous step and run Unpacker (+ optional converter).

process.Unpacker = cms.EDProducer("Phase2RawToClusterProducer@alpaka",
#process.Unpacker = cms.EDProducer("alpaka_serial_sync::Phase2RawToClusterProducer",
#process.Unpacker = cms.EDProducer("alpaka_cuda_async::Phase2RawToClusterProducer",
#process.Unpacker = cms.EDProducer("alpaka_rocm_async::Phase2RawToClusterProducer",
    fedRawDataCollection = cms.InputTag("Packer"),
)

process.ClusterConverter = cms.EDProducer("ClusterPropSoAToLegacyED",
    clusterSoASource = cms.InputTag("Unpacker")
)

process.out = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep FEDRawDataCollection_*_*_*',
        'keep *_ClustersFromPhase2TrackerDigis_*_*',
        'keep *_Packer_*_*',
        'keep *_Unpacker_*_*',
        'keep *_mix_Tracker_*',
        'keep *_ClusterConverter_*_*'
    ),
    fileName = cms.untracked.string('unpack_convert_out.root')
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
# (no effect here since we don't run the Clusterizer in this file, but kept identical)
# premix_stage2.toModify(process.ClustersFromPhase2TrackerDigis, rawHits = ["mixData:Tracker"])

process.Timing = cms.Service("Timing",
    summaryOnly = cms.untracked.bool(True),
    useJobReport = cms.untracked.bool(True)
)
# mark framework transitions in the NVIDIA profiler
process.NVProfilerService = cms.Service("NVProfilerService",
    showModulePrefetching = cms.untracked.bool(False)
)

# ------------------ SPLIT HERE: ONLY run Unpacker (+ Converter) ------------------
if Legacy_Format:
    process.dtc = cms.Path(
        process.Unpacker *
        process.ClusterConverter
    )
else:
    process.dtc = cms.Path(
        process.Unpacker
    )

process.output = cms.EndPath(process.out)

# mark framework transitions in the NVIDIA profiler
process.NVProfilerService = cms.Service("NVProfilerService",
    showModulePrefetching = cms.untracked.bool(False)
)
