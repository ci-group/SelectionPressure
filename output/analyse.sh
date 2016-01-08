#!/bin/bash
#
# $Id: analyse.sh 165 2016-01-08 11:35:17Z ehaasdi $
#
# For defined experiments (see definition of EXPERIMENTS below), generate summary statistics over multiple runs and plot them for experiments processed with selection_pressure.sh
#
#BASEDIR=~/Documents/Projecten/QuantifyingSP/
BASEDIR=~/projects/workspace/SelectionPressure/output/

AWKDIR=~/projects/scripts/awk

TEMP=`mktemp XXXX`

# $1: experiment name
# $2: objective postfix. May be empty
function analyse_FET {

	if [ -z "$2" ]
	then
  		POSTFIX_IN=.1
	else
		POSTFIX_IN=$2
	fi

	RAW=$1.FET$2

    paste -d ' ' ${SRC_DIR}/$1.*${POSTFIX_IN}.FET.txt | awk '{printf("%d",$1); for (i=2; i<=NF; i=i+2) { printf("%s %f", FS, -(log($i)/log(10))) } printf("%s", RS) } ' > $RAW

	gawk -v skip=1 -f ${AWKDIR}/moments-per-line.awk $RAW > $TEMP
	awk 'BEGIN{print "#generation"}{print $1}'  ${SRC_DIR}/$1.1${POSTFIX_IN}.FET.txt | paste - $TEMP > $RAW.moments
	gnuplot -e "set term pngcairo; set output '$RAW.png'; set title '$1$2'; set ylabel 'FET'; raw_input='$RAW'; moments_input='$RAW.moments'; nr_observations=50; skip=100" $BASEDIR/raw_and_average.gpl
}

# $1: experiment name
# $2: objective postfix. May be empty
function analyse_rankbased {

	if [ -z "$2" ]
	then
  		POSTFIX_IN=.1
	else
		POSTFIX_IN=$2
	fi

	RAW=$1.rankbased$2

    paste -d ' ' ${SRC_DIR}/$1.*${POSTFIX_IN}.rankbased.txt | awk '{printf("%d",$1); for (i=2; i<=NF; i=i+2) { printf("%s %f", FS, $i) } printf("%s", RS) } ' > $RAW

	gawk -v skip=1 -f ${AWKDIR}/moments-per-line.awk $RAW > $TEMP
	awk 'BEGIN{print "#generation"}{print $1}'  ${SRC_DIR}/$1.1${POSTFIX_IN}.rankbased.txt | paste - $TEMP > $RAW.moments
	gnuplot -e "set term pngcairo; set output '$RAW.png'; set title '$1$2'; set ylabel 'rankbased'; set yrange [0:0.6];raw_input='$RAW'; moments_input='$RAW.moments'; nr_observations=50; skip=25" $BASEDIR/raw_and_average.gpl
}

#SRC_DIR=analysis
SRC_DIR=./

EXPERIMENTS="1.0.2 1.0.5 1.0.10 1.1.2 1.1.5 1.1.10 1.2.0"
#0.0.10 0.0.2 0.0.5 0.1.10 0.1.2 0.1.5 0.2.0 1.0.10 1.0.2 1.0.5 1.1.10 1.1.2 1.1.5 1.2.0"
#"bugs"

# for EXPERIMENT in $EXPERIMENTS
# do
#  	echo -n $EXPERIMENT ": FET..."
#  	analyse_FET $EXPERIMENT
#  	echo "rank-based"
#  	analyse_rankbased $EXPERIMENT
# done

#MULTI_OBJ="2.short 3.short monee"
MULTI_OBJ="monee"
for EXPERIMENT in $MULTI_OBJ
do
	echo -n $EXPERIMENT "FET..."
 	analyse_FET $EXPERIMENT .1
 	analyse_FET $EXPERIMENT .2

  	echo "rank-based"
 	analyse_rankbased $EXPERIMENT .1
 	analyse_rankbased $EXPERIMENT .2
done

rm $TEMP
