tournament=`head $1 | grep 'Tournament size:'  | awk '{print $4}'`

tmpfile=`mktemp XXXX`

BASEDIR=~/projects/workspace/SelectionPressure/output/

for generation in 1 2 5 10 20 50 100 200 500 1000 5000 10000 20000 49999
do
	grep "^$generation " $1  | sort -n -k 4 > $tmpfile
	gnuplot -e "set label 1 'generation $generation' at 1,9.7; \
	set title 'tournament size $tournament';\
	set output '$1.$generation.png';\
	input='$tmpfile'" \
	$BASEDIR/sarmadejong.gpl
done

rm -f $tmpfile