#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Usage: ./go_plots.sh Jan15(the date)"
	echo
	exit 0
fi

folder=$1

for x in `ls -d *_$folder | grep -v "SingleElectron" | grep -v "SingleMu"`
do

	[ "$(ls $x/signal_contamination_*.root)" ] || continue

	[ -f $x/signal_contamination_${x%_$folder}.root ] && continue

	if [ `ls $x/signal_contamination_*.root | wc -l` -eq `ls $x/filelist_* | wc -l` ]
	then
		echo would move signal_contamination_${x%_$folder}.root somewhere
		cd $x
		./combineStagedOutput.sh
		cp signal_contamination_${x%_$folder}.root ../../makePlots/inputs/
		cd ..
		rm -r $x
	fi

done
