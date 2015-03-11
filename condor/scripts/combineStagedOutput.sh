#!/bin/bash

WORK_DIR=`pwd`

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

hadd signal_contamination_DATASETNAME.root signal_contamination_DATASETNAME_*.root
rm signal_contamination_DATASETNAME_*.root
