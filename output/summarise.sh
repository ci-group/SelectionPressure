#!/bin/bash

# Command line flags
. ${HOME}/opt/lib/shflags

DEFINE_string 'method' 'fisher' "Method to calculate selection pressure. Either 'fisher' or  'spearman'"
DEFINE_string 'output' 'out' "Base name for output"

FLAGS "$@" || exit 1
eval set -- "${FLAGS_ARGV}"

let firsttime=1

paste_buffer=`mktemp statsXXXX`

#
# Collate all runs
#
for run in $@
do
	if [ $firsttime -eq 1 ]
	then
		firsttime=0

		if [ ${FLAGS_method} == 'fisher' ]
		then
			awk '/^[0-9]*00 /{print $1, -log($2)}' $run > $paste_buffer
		fi
		if [ ${FLAGS_method} == 'spearman' ]
		then
			awk '/^[0-9]*00 /{print $1, -($2)}' $run > $paste_buffer
		fi

	else
		paste_buffer_2=`mktemp statsXXXX`

		if [ ${FLAGS_method} == 'fisher' ]
		then
			awk '/^[0-9]*00 /{print -log($2)}' $run | paste $paste_buffer - > $paste_buffer_2
		fi
		if [ ${FLAGS_method} == 'spearman' ]
		then
			awk '/^[0-9]*00 /{print -($2)}' $run  | paste $paste_buffer - > $paste_buffer_2
		fi

		mv $paste_buffer_2 $paste_buffer
	fi
done

output=${FLAGS_output}
#
# Calculate statistics per timestep
#
paste_buffer2=`mktemp statsXXXX`
echo "#generation" | cat - $paste_buffer | cut -f 1 -d ' ' > $paste_buffer2
cut -f 2- -d ' ' $paste_buffer | gawk -f ~/projects/awk/moments-per-line.awk | paste -d ' ' $paste_buffer2 - > $output.stats

rm -f $paste_buffer2 $paste_buffer

gnuplot -e "set terminal pngcairo enhanced font 'Helvetica,11' size 1024, 756; set output '$output.png'; set style fill transparent solid 0.2 border; plot '$output.stats' u 1:4:5 every 10 w filledcu lc rgb('#2B6088') title '$output', '' u 1:3 ps 0.5 lc rgb('#2B6088') notitle; quit"