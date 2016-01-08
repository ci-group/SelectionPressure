BASEDIR=~/projects/workspace/SelectionPressure/output/

for problem in 0 1
do
	for scheme in 0 1
	do
		for tournament in 2 5 10
		do
		(
			log=$problem.$scheme.$tournament
			bash $BASEDIR/summarise.sh --method=fisher --output $log.FET $log.*.FET.txt
			bash $BASEDIR/summarise.sh --method=spearman --output $log.spearman $log.*.spearman.txt
		)&
		done

		log=$problem.2.0
		(
		bash $BASEDIR/summarise.sh --method=fisher --output $log.FET $log.*.FET.txt
		bash $BASEDIR/summarise.sh --method=spearman --output $log.spearman $log.*.spearman.txt
		)&
	done

	wait
done