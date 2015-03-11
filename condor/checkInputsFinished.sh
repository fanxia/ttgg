#!/bin/bash

folder=$1

for x in `ls -d *_$folder | grep -v "SingleElectron" | grep -v "SingleMu"`
do

	[ "$(ls $x/pileupReweighting_*.root)" ] || continue
	[ "$(ls $x/btagEfficiency_*.root)" ] || continue

	[ -f $x/btagEfficiency_${x%_$folder}.root ] && continue
	[ -f $x/pileupReweighting_${x%_$folder}.root ] && continue


	if [ `ls $x/pileupReweighting_*.root | wc -l` -eq `ls $x/filelist_* | wc -l` ] && [ `ls $x/btagEfficiency_*.root | wc -l` -eq `ls $x/filelist_* | wc -l` ]
	then
		cd $x
		./combineStagedInput.sh
		cp pileupReweighting*.root ../saved_weights/
		cp btagEfficiency_*.root ../saved_weights/
		condor_submit acceptance_go_mc.jdl
		cd ..
	fi
done
