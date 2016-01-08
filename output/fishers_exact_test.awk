#!/opt/local/bin/gawk
#
# $Id: $
#
# Output for R's Fischer Exact Test as follows:
#
# 				Low/Equal Fitness	High Fitness
# Low/Equal Off			A 				B
# High Off				C 				D
#
# Matrix per time slice
#
# Input: generation Id offspring fitness(es)
# Output: generation A B C D
#
#
# The objective variable defines which of the fitnesses to use;
# specify through "-v objective=1". If unspecified, defaults to 1,
# specify whether the problem is a min or max problem. when max=1 it is a max problem
# so that works for single-objective experiments.
#

BEGIN{
	# If objective nr is not specified, default to 1
	if (objective == 0)
		objective = 1


	if (max){
		do_max_term = 1 # max problem
		} else {
		do_max_term = -1 # min problem, so multiply fitness with -1
		}

	lastGeneration = -1
}

!/^\#/{
	i = $1

	if (lastGeneration != i) {
		lastGeneration = i
		generations[maxGen++] = i
	}

	sizes[i]++
#	Ids[i][sizes[i]] = $2
	offspring[i][sizes[i]] = $3;
	fitnesses[i][sizes[i]] = $(3+objective);
}
# Note: 0-based basic arrays, but subarrays 1-based (for asort compatibility)

END {
	# calculate median values per generation and fill matrices for FET
	for (g in generations) {
		i = generations[g]

		if (!(i in fitnesses)) continue

		n = sizes[i];

		#
		# Calculate threshold values: sort fitness and offspring arrays
		# and set threshold values as specified by percentile.
		#
		percentile = 2  # 2: median, 4: quartile, etc.
		thresholdIndex = int(n/percentile) + 1;

    	asort(fitnesses[i], fit);
 		thresholdFitness = (n % percentile) ?
 			fit[thresholdIndex] :
 			(fit[thresholdIndex] + fit[thresholdIndex-1])/2;

    	asort(offspring[i], children);
    	thresholdOffspring = (n % percentile) ?
    		children[thresholdIndex] :
    		(children[thresholdIndex] + children[thresholdIndex-1])/2;

		#
		# If median offspring count is 0, set it to 1 (or all individuals
		# are on the same side of the cutoff)
		#
		if (thresholdOffspring == 0)
			thresholdOffspring = 1;

		a[i] = 0;
		b[i] = 0;
		c[i] = 0;
		d[i] = 0;

		for (j=1; j<=n; j++)
		{
#			print Ids[i][j], offspring[i][j], fitnesses[i][j] > "/dev/stderr"
			if (offspring[i][j] <= thresholdOffspring)  # Low offspring
			{
				if (fitnesses[i][j]*do_max_term > thresholdFitness*do_max_term)
					b[i]++;	# high fitness, low offspring
				else
					a[i]++;	# low fitness, low offspring

			} else {	# High offspring
				if (fitnesses[i][j]*do_max_term > thresholdFitness*do_max_term)
					d[i]++;	# high fitness, high offspring
				else
					c[i]++;	# low fitness, high offspring

			}
		}


		print i, a[i], b[i], c[i], d[i]
	}
}