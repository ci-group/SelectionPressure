#!/bin/bash
#
# $Id: run.sh 164 2015-10-12 08:12:34Z ehaasdi $
#
# Runs independent repeats of SelectionPressure experiment for all problems and analyses the individual runs.
#

OUTPUT_RAW=/Volumes/Slartibartfast/experiments/selection_pressure/results
OUTPUT_ANALYSIS=/Volumes/Slartibartfast/experiments/selection_pressure/analysis

for run in {1..50}
do
	#Yasga runs here
	for problem in 1
	do
		# Tournament based schemes
		for scheme in 0 1
		do
			for tournament in 2 5 10
			do
			(
				log=${OUTPUT_RAW}/$problem.$scheme.$tournament.$run
				../Debug/SelectionPressure --s $scheme --t $tournament --p $problem --l $log
				bash selection_pressure.sh $log --output $OUTPUT_ANALYSIS  --method both  --optimisation min
			)&
			done
		done

		# Fitness-proportionate selection
		(
			log=${OUTPUT_RAW}/$problem.2.0.$run
			../Debug/SelectionPressure --s 2 --p $problem --l $log
			bash selection_pressure.sh $log --output $OUTPUT_ANALYSIS  --method both --optimisation min
		)&

#		wait
	done

    #nsga runs here
    (
	log=${OUTPUT_RAW}/2.$run
    ../Debug/SelectionPressure --p 2 --l $log
	bash selection_pressure.sh $log --output $OUTPUT_ANALYSIS --method both --optimisation min
	)&

	(
	log=${OUTPUT_RAW}/3.$run
    ../Debug/SelectionPressure --p 3 --l $log
	bash selection_pressure.sh $log --output $OUTPUT_ANALYSIS --method both --optimisation min
    )&

	wait
done
