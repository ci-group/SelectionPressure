/*****************************************************************************
 *   Copyright (C) 2014 VU University Amsterdam                              *
 *   Based on PaGMO sga implementation                                       *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *
 *   $Id: yasga.cpp 164 2015-10-12 08:12:34Z ehaasdi $
 *
 *****************************************************************************/

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/round.hpp>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include <pagmo/exceptions.h>
#include <pagmo/population.h>
#include <pagmo/types.h>
#include <pagmo/algorithm/base.h>
#include "yasga.h"
#include <pagmo/problem/base_stochastic.h>

#define _VERSION_H_AS_HEADER_
#include "version.h"

using namespace pagmo;
using namespace pagmo::algorithm;

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability (of each allele if binomial crossover)
 * @param[in] m Mutation probability (of each allele)
 * @param[in] elitism The best individual is reinserted in the population each elitism generations
 * @param[in] mut Mutation type. One of sga::mutation::GAUSSIAN, sga::mutation::RANDOM
 * @param[in] width Mutation width. When gaussian mutation is selected is the width of the mutation
 * @param[in] tournament_size Size of the tournament for tournament selection. A value of 1 boils down to random selection.
 * @param[in] cro Crossover type. One of sga::crossover::BINOMIAL, sga::crossover::EXPONENTIAL
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1]\f$, mutation probability is not \f$ \in [0,1]\f$,
 * elitism is <=0, tournament_size < 1
 *
 */

yasga::yasga(int gen, const double &cr, const double &m, int elitism, mutation::type mut, double width, pagmo::population::size_type tournament_size,
	    selection_operator selection, double initial_temperature, double temperature_decrement, crossover::type cro)
	:base(),
	 m_gen(gen),
	 m_cr(cr),
	 m_m(m),
	 m_elitism(elitism),
	 m_mut(mut,width),
	 m_tournament_size(tournament_size),
	 m_selection_operator(selection),
	 m_initial_temperature(initial_temperature),
	 m_temperature_decrement(temperature_decrement),
	 m_cro(cro),
	 m_logfile(std::string("offspring.log"))
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (cr > 1 || cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1] range");
	}
	if (m < 0 || m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (elitism < 1) {
			pagmo_throw(value_error,"elitisim must be greater than zero");
	}
	if (tournament_size < 1) {
		pagmo_throw(value_error,"tournament_size must be greater than zero");
	}
	if (width <0 || width >1) {
		pagmo_throw(value_error,"mutation width must be in the [0,1] range");
	}

	if (m_initial_temperature == 0) {
		m_initial_temperature = gen;
	}
}

/// Clone method.
base_ptr yasga::clone() const
{
	return base_ptr(new yasga(*this));
}

/// Evolve implementation.
/**
 * Run the simple genetic algorithm for the number of generations specified in the constructors.
 * At each improvment velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void yasga::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), Di = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - Di;


	//We perform some checks to determine wether the problem/population are suitable for SGA
	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and SGA is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and SGA is not suitable to solve it");
	}

	if (NP < 5) {
		pagmo_throw(value_error,"for SGA at least 5 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy), Xnew(NP,dummy);

	std::vector<fitness_vector > fit(NP);		//fitness

	fitness_vector bestfit;
	decision_vector bestX(D,0);

	std::vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);
	std::vector <int> selection(NP);

	std::vector<int> fitnessID(NP);

	// Initialise the chromosomes and their fitness to that of the initial deme
	for (pagmo::population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Find the best member and store in bestX and bestfit
	double bestidx = pop.get_best_idx();
	bestX = pop.get_individual(bestidx).cur_x;
	bestfit = pop.get_individual(bestidx).cur_f;

	std::ofstream offspring_log(m_logfile.c_str());

	offspring_log
		<< "# Generated w. version " << version.v_long
		<< "\n# Selection scheme: " << m_selection_operator
		<< "\n# Tournament size: " << m_tournament_size
		<< "\n# Gen Id Offspring Fitness\n# Problem: "
		<< pop.problem().get_name() << '\n';

	// Main SGA loop
	for (int j = 0; j<m_gen; j++) {
//		if (j %100 ==0) std::cout << '.' << std::flush;

		// For logging offspring counts
		std::vector<unsigned > offspring_counts(NP,0);

		double temperature(m_initial_temperature - (m_temperature_decrement * j));

		if (m_selection_operator == ROULETTE) {
			// Preparations for fitness proportionate selection
			// We scale all fitness values from 0 (worst) to absolute value of the best fitness
			fitness_vector worstfit=fit[0];
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if (prob.compare_fitness(worstfit,fit[i])) worstfit=fit[i];
			}

			for (pagmo::population::size_type i = 0; i < NP; i++) {
				selectionfitness[i] = fabs(worstfit[0] - fit[i][0]);
			}

			// We build and normalise the cumulative sum
			cumsumTemp[0] = selectionfitness[0];
			for (pagmo::population::size_type i = 1; i< NP; i++) {
				cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
			}
			for (pagmo::population::size_type i = 0; i < NP; i++) {
				cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
			}
		}

		for (pagmo::population::size_type i=0; i<NP; ++i) {
			// Generate m_tournament_size random indices into population, select the one with the best fitness
			// Note that an individual may be entered into a tournament more than once
			pagmo::population::size_type selIndex = boost::uniform_int<int>(0,NP - 1)(m_urng);
			fitness_vector tmpFit = fit[selIndex];
			selection[i] = selIndex;
			for (pagmo::population::size_type k = 0; k < m_tournament_size - 1; ++k) {
				selIndex = boost::uniform_int<int>(0,NP - 1)(m_urng);

				switch (m_selection_operator) {
				case TOURNAMENT: {
					if (prob.compare_fitness(fit[selIndex], tmpFit)) {
						tmpFit = fit[selIndex];
						selection[i] = selIndex;
					}
					break;
				}
				case BOLTZMANN: {

					//
					// First, make sure that current selection is better of the two;
					// we may select the worse one after all through boltzmann comparison
					//
					if (prob.compare_fitness(fit[selIndex], tmpFit)) {
						tmpFit = fit[selIndex];
						selection[i] = selIndex;
					}

					// FIXME: should have difference calculated by the problem!! Now assume single objective, minimisation
					if (boltzmann_compare(tmpFit[0], fit[selIndex][0], temperature)) {
						tmpFit = fit[selIndex];
						selection[i] = selIndex;
					}
					break;
				}

				case ROULETTE: {
					//we throw a die and pick up the corresponding index
					double r2 = m_drng();

					for (pagmo::population::size_type selIndex = 0; selIndex < NP; selIndex++) {
						if (cumsum[selIndex] > r2) {
							selection[i]=selIndex;
							break;
						}
					}
					break;
				}


				default:
					pagmo_throw(value_error,"unknown selection operator");
				}
			}
		}

		//Xnew stores the new selected generation of chromosomes
		for (pagmo::population::size_type i = 0; i < NP; ++i) {
			Xnew[i]=X[selection[i]];

			// Count # offspring
			offspring_counts[selection[i]]++;
		}

		//2 - Crossover
		{
			int r1,L;
			decision_vector  member1,member2;

			for (pagmo::population::size_type i=0; i< NP; i++) {
				//for each chromosome selected i.e. in Xnew
				member1 = Xnew[i];
				//we select a mating patner different from the self (i.e. no masturbation)
				do {
					r1 = boost::uniform_int<int>(0,NP - 1)(m_urng);
				} while ( r1 == boost::numeric_cast<int>(i) );
				member2 = Xnew[r1];

				// Count # offspring for 2nd parent
				offspring_counts[r1]++;

				//and we operate crossover
				switch (m_cro) {
					//0 - binomial crossover
				case crossover::BINOMIAL: {
					size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
					for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
						if ((m_drng() < m_cr) || L + 1 == D) { /* change at least one parameter */
							member1[n] = member2[n];
						}
						n = (n+1)%D;
					}
					break; }
					//1 - exponential crossover
				case crossover::EXPONENTIAL: {
					size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
					L = 0;
					do {
						member1[n] = member2[n];
						n = (n+1) % D;
						L++;
					}  while ( (m_drng() < m_cr) && (L < boost::numeric_cast<int>(D)) );
					break; }
				}
				Xnew[i] = member1;

			} }

		//3 - Mutation
		switch (m_mut.m_type) {
		case mutation::GAUSSIAN: {
			boost::normal_distribution<double> dist;
			boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > delta(m_drng,dist);
			for (pagmo::problem::base::size_type k = 0; k < Dc;k++) { //for each continuous variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = (delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
					}
				}
			}
			for (pagmo::problem::base::size_type k = Dc; k < D;k++) { //for each integer variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = boost::math::iround(delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
					}
				}
			}
			break;
			}
		case mutation::RANDOM: {
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				for (pagmo::problem::base::size_type j = 0; j < Dc;j++) { //for each continuous variable
					if (m_drng() < m_m) {
						Xnew[i][j] = boost::uniform_real<double>(lb[j],ub[j])(m_drng);
					}
				}
				for (pagmo::problem::base::size_type j = Dc; j < D;j++) {//for each integer variable
					if (m_drng() < m_m) {
						Xnew[i][j] = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
					}
				}
			}
			break;
			}
		}

		// Log offspring
		for (pagmo::population::size_type i = 0; i < NP; ++i) {
			offspring_log << j << ' ' << i << ' ' << offspring_counts[i] << ' ' << fit[i][0] << '\n';
		}

		// If the problem is a stochastic optimization chage the seed and re-evaluate taking care to update also best and local bests
		try
		{
			//4 - Evaluate the new population (stochastic problem)
			dynamic_cast<const pagmo::problem::base_stochastic &>(prob).set_seed(m_urng());
			pop.clear(); // Removes memory based on different seeds (champion and best_x, best_f, best_c)
			
			// We re-evaluate the best individual (for elitism)
			prob.objfun(bestfit,bestX);
			// Re-evaluate wrt new seed the particle position and memory
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				// We evaluate here the new individual fitness
				prob.objfun(fit[i],Xnew[i]);
				// We update the velocity (in case coupling with PSO via archipelago)
				//dummy = Xnew[i];
				//std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				///We now set the cleared pop. cur_x is the best_x, re-evaluated with new seed.
				pop.push_back(Xnew[i]);
				//pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		catch (const std::bad_cast& e)
		{
			//4 - Evaluate the new population (deterministic problem)
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				prob.objfun(fit[i],Xnew[i]);
				dummy = Xnew[i];
				std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				//updates x and v (cache avoids to recompute the objective function and constraints)
				pop.set_x(i,Xnew[i]);
				pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		
		//5 - Reinsert best individual every m_elitism generations
		if (j % m_elitism == 0) {
			int worst=0;
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if ( prob.compare_fitness(fit[worst],fit[i]) ) worst=i;
			}
			Xnew[worst] = bestX;
			fit[worst] = bestfit;
			dummy = Xnew[worst];
			std::transform(dummy.begin(), dummy.end(), pop.get_individual(worst).cur_x.begin(), dummy.begin(),std::minus<double>());
			//updates x and v (cache avoids to recompute the objective function)
			pop.set_x(worst,Xnew[worst]);
			pop.set_v(worst,dummy);
		}
		X = Xnew;

	} // end of main SGA loop
}

/// Algorithm name
std::string yasga::get_name() const
{
	return "A Simple Genetic Algorithm";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string yasga::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "CR:" << m_cr << ' ';
	s << "M:" << m_m << ' ';
	s << "elitism:" << m_elitism << ' ';
	s << "mutation:";
	switch (m_mut.m_type) {
		case mutation::RANDOM: {
		      s << "RANDOM "; 
		      break;
		      }
		case mutation::GAUSSIAN: {
		      s << "GAUSSIAN (" << m_mut.m_width << ") "; 
		      break;
		      }
	}
	s << "tournament_size:" << m_tournament_size << ' ';
	s << "crossover:";
	switch (m_cro) {
		case crossover::EXPONENTIAL: {
		      s << "EXPONENTIAL "; 
		      break;
		      }
		case crossover::BINOMIAL: {
		      s << "BINOMIAL "; 
		      break;
		      }
	}

	return s.str();
}


/**
 *  Comparison function for Boltzmann selection
 *
 *  Assumes minimisation problem
 */
bool yasga::boltzmann_compare(double champion_fitness, double challenger_fitness, double temperature) const {

	double diff(challenger_fitness - champion_fitness);

	if (diff <= 0.0) // Challenger beats champion
	 return true;
	//else

	if (temperature <= 0.0) // Annealing has ended: only select improvements
	 return false;
	// else

	double rnd(m_drng());
//	std::cout << rnd << ' ' << exp( - diff/temperature ) <<' '<< diff <<' '<< temperature << std::endl;
	return (rnd < exp( - diff/temperature ));
}

