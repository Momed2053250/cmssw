#!/bin/bash
#
# wrapper to run the CMSSW cfg under HTCondor
#

# 1) bring in CMS environment
#source /cvmfs/cms.cern.ch/cmsset_default.sh

# 2) go to your CMSSW area
cd /afs/cern.ch/user/m/mmomed/unpacker-13-05-25/CMSSW_15_0_4/src/
cmsenv

# 3) run the cfg (assumes it lives in the condor scratch dir)
cmsRun /afs/cern.ch/user/m/mmomed/unpacker-13-05-25/CMSSW_15_0_4/src/EventFilter/Phase2TrackerRawToDigi/test/SLinkProducerAndUnpacker_alpaka_cfg.py

