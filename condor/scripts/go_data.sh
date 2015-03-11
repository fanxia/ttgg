#!/bin/bash

JOB_NUMBER=$1
WORK_DIR=`pwd`

tar -xzf fileLists.tgz

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3/src/
mv $WORK_DIR/src.tgz .
tar -xzf src.tgz
eval `scramv1 runtime -sh`

mv $WORK_DIR/ANALYZER .
mv $WORK_DIR/filelist_$JOB_NUMBER .

while read file
do
  	sed -i '9i chain.Add("'$file'");' ANALYZER
done < filelist_$JOB_NUMBER

root -b -q -l ANALYZER

mv hist_analysis_CSVM.root $WORK_DIR/hist_analysis_CSVM_$JOB_NUMBER.root
cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
rm fileLists.tgz
rm filelist_*
