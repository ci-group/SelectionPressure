#!/bin/bash
#
# $Id: selection_pressure.sh 164 2015-10-12 08:12:34Z ehaasdi $
#
# Quantify selection pressure using either Fisher's exact test (FET) or rank-based
# correlation (Kendall's tau-b).
#
# Call as: selection_pressure.sh <input file(s)> <options>
#
# Expected input files layout:
# Time(generation) Id Offspring Fitness(es)
#
#BASEDIR=~/Documents/Projecten/QuantifyingSP/
BASEDIR=~/projects/workspace/SelectionPressure/output/

. ${HOME}/opt/lib/shflags

#define the command line options
DEFINE_boolean 'plotOnly' false "Don't regenerate data, just generate a plot"
DEFINE_boolean 'skip10' true "Only analyse every 10th generation (useful for long runs)"
DEFINE_string 'output' '.' "Path of output directory. Defaults to current working directory."
DEFINE_string 'method' 'both' "Method to calculate selection pressure. One of 'fisher', 'rank', both"
DEFINE_string 'optimisation' 'max' "Is it a maximisation or minimisation problem ('max', 'min')"

# Parse the flags
FLAGS "$@" || exit 1
eval set -- "${FLAGS_ARGV}"

do_fisher=${FLAGS_FALSE}
do_rankbased=${FLAGS_FALSE}
if [ ${FLAGS_method} == 'fisher' ]
then
	do_fisher=${FLAGS_TRUE}
fi
if [ ${FLAGS_method} == 'rank' ]
then
	do_rankbased=${FLAGS_TRUE}
fi
if [ ${FLAGS_method} == 'both' ]
then
	do_rankbased=${FLAGS_TRUE}
	do_fisher=${FLAGS_TRUE}
fi

do_max=0
if [ ${FLAGS_optimisation} == 'max' ]
then
	do_max=1
fi
if [ ${FLAGS_optimisation} == 'min' ]
then
	do_max=0
fi

if [ ${FLAGS_skip10} -eq ${FLAGS_TRUE} ]
then
	RECORD_FILTER='^[0-9]*0 '
else
	RECORD_FILTER='^[0-9]* '
fi


for experiment in $@
do
(
	if [[ $file =~ \.bzip2$ ]] #check for bzipped results
	then
		experiment_name=`basename $experiment .bzip2`
		GREP=bzgrep
		CAT=bzcat
	else
		experiment_name=`basename $experiment`
		GREP=grep
		CAT=cat
	fi
	#
	# Count number of objectives
	#
	NR_COLS=`$CAT $experiment | head -n 20 | grep '^[0-9]' | head -n 1 | wc -w`
	let "NR_OBJECTIVES = $NR_COLS - 3"

	for objective in `seq 1 $NR_OBJECTIVES`
	do
		if [ ${FLAGS_plotOnly} -eq ${FLAGS_FALSE} ]
		then

			if [ ${do_fisher} -eq ${FLAGS_TRUE} ]
			then
				echo "Fisher-based analysis of" $experiment ", objective nr. " $objective

				#
				# FET requires 2x2 matrix. Classes we use: above or below median
				# offspring counts, above or below median fitness
				# The fishers_exact_test.awk script generates these counts.
				#
				# grep filters out comments
				#
				# Awk script expected pressure-stats columns:
				# Time (e.g. Generation), Id, offspring, fitness(es)
				#
				tempfile=`mktemp ${TMPDIR}/FETXXXX`
                ${GREP} "${RECORD_FILTER}" $experiment | gawk -v max=$do_max -v objective=$objective -f $BASEDIR/fishers_exact_test.awk  > $tempfile

				# Generate script to call R's Fisher's exact test for calculation
				rscript=`mktemp ${TMPDIR}/rscriptXXXX`

				echo 'theTest <- function(x) {' > $rscript
				echo '	result = fisher.test(matrix(c(x["A"],x["C"],x["B"],x["D"]),nrow=2), alternative = "greater")' >> $rscript
				echo '	cat(paste(x["time"], result["p.value"],"\n", sep=" "))' >> $rscript
				echo '}' >> $rscript
				echo "myData <- read.table(\"$tempfile\", sep=\" \", col.names=c(\"time\", \"A\", \"B\", \"C\", \"D\"))" >> $rscript
				echo 'apply(myData, 1, theTest)' >> $rscript

				# somehow, the R code adds 'NULL' line in some cases, use sed to filter
				# that out. sort command makes sure that output is sorted by
				# timestep/generation
				R --no-save --slave < $rscript | sed '/NULL/d' | sort -n > ${FLAGS_output}/$experiment_name.$objective.FET.txt

                rm $tempfile
				rm $rscript
			fi

			if [ ${do_rankbased} -eq ${FLAGS_TRUE} ]
			then
				echo "rank-based analysis of" $experiment ", objective nr. " $objective

				#
				# use awk to strip comments and select objective
				#
				tempfile=`mktemp ${TMPDIR}/rankbasedXXXX`
				${CAT} $experiment | gawk -v objective=$objective "/${RECORD_FILTER}/{print \$1, \$2, \$3, \$(3+objective)}" > $tempfile

				#
				# Generate script to call R's kendall test for calculation
				#
				rscript=`mktemp ${TMPDIR}/rscriptXXXX`
                #echo 'install.packages("Kendall")' >>$rscript
				echo 'names <- c("generation", "Id", "offspring", "fitness")' > $rscript
				echo "myData <- read.table(\"$tempfile\", sep=\" \", col.names=names)" >> $rscript
                echo 'test <- function(x) { Kendall(x$offspring, x$fitness)[1] }' >> $rscript
                #echo 'test <- function(x) { cor(x$offspring, x$fitness, method="kendall") }' >> $rscript
				echo 'coefs <- by(myData, myData$generation, test)' >> $rscript
				echo 'cat(paste(names(coefs), coefs), sep = "\n")' >> $rscript

				# Added sed to filter out spurious warnings by Kendall package
				R --no-save --slave < $rscript | sed '/^WARNING/d' > ${FLAGS_output}/$experiment_name.$objective.rankbased.txt

				rm $tempfile
				rm $rscript
			fi
		fi

# 		if [ ${do_fisher} -eq ${FLAGS_TRUE} ]
# 		then
# 			gnuplot -e "set term pngcairo; set output '${FLAGS_output}/$experiment_name.$objective.selective-pressure.fisher.png'; input='$experiment.$objective.FET.txt';set title '$experiment'" $BASEDIR/fisher.gpl
# 		fi
#
# 		if [ ${do_rankbased} -eq ${FLAGS_TRUE} ]
# 		then
# 			gnuplot -e "set term pngcairo; set output '$experiment.$objective.selective-pressure.rankbased.png'; input='${FLAGS_output}/$experiment_name.$objective.rankbased.txt';set title '$experiment'" $BASEDIR/spearman.gpl
# 		fi
	done
	) &
done

wait

echo "Analysis complete"