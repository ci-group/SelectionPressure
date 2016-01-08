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
 *   $Id: yasga.h 158 2014-12-03 12:43:55Z ehaasdi $
 *
 *****************************************************************************/

#ifndef _YASGA_H
#define _YASGA_H

#include <pagmo/config.h>
#include <pagmo/problem/base.h>
#include <pagmo/serialization.h>
#include <pagmo/algorithm/base.h>

/// Yet Another Simple Genetic Algorithm (YASGA)
/**
 * Evert's modified version of Dario Izzo's SGA (part of base PaGMO).
 *
 * Modified to have tournament selection.
 *
 * From the oOriginal documentation:
 *
 * Genetic algorithms are very popular algorithms used widely by people of very different backgrounds.
 * As a consequence there are a large number of different implementations and toolboxes that are available
 * and can be used to construct a genetic algorithm. We decided not to choose one of these and, instead, to
 * provide only a basic implementation of the algorithm implementing a floating point encoding (not binary)
 * and some common mutation and crossover strategies, hence the name Simple Genetic Algorithm.
 *
 * Mutation is gaussian or random, crossover exponential or binomial.
 *
 * The algorithm works on single objective, box constrained problems. The mutation operator acts
 * differently on continuous and discrete variables.
 *
 * @author Dario Izzo (dario.izzo@googlemail.com)
 * @author Evert Haasdijk (e.haasdijk@vu.nl)
 *
 */

class yasga: public pagmo::algorithm::base
{
public:

	/// Mutation operator info
	struct mutation {
			/// Mutation type, gaussian or random
			enum type {GAUSSIAN = 0, RANDOM = 1};
			/// Constructor
			/**
			* \param[in] t the mutation type
			* \param[in] width the width of the gaussian bell in case of a gaussian mutation. The
			*		parameter is otherwise ignored. width is a percentage with respect to the
			*		ub[i]-lb[i] width.
			*/
			mutation(mutation::type t, double width) : m_type(t),m_width(width) {}
			/// Mutation type
			type m_type;
			/// Mutation width
			double m_width;
		private:
			friend class boost::serialization::access;
			template <class Archive>
			void serialize(Archive &ar, const unsigned int)
			{
				ar & m_type;
				ar & m_width;
			}
	};

	/// Crossover operator info
	struct crossover {
		/// Crossover type, binomial or exponential
		enum type {BINOMIAL = 0, EXPONENTIAL = 1};
	};

	/// Selection operator
	enum selection_operator {
		/// Regular rank-based tournament selection
		TOURNAMENT = 0,
		/// Tournament with Boltzmann selection
		BOLTZMANN = 1,
		/// Fitness-proportionate selection
		ROULETTE = 2
	};

	yasga(int gen  = 1, const double &cr = .95, const double &m = .02, int elitism = 1,
	    mutation::type mut  = mutation::GAUSSIAN, double width = 0.1,
	    pagmo::population::size_type tournament_size = 2,
	    selection_operator selection = TOURNAMENT,
	    double initial_temperature = 0,
	    double temperature_decrement = 1,
	    crossover::type cro = crossover::EXPONENTIAL);

	pagmo::algorithm::base_ptr clone() const;
	void evolve(pagmo::population &) const;
	std::string get_name() const;

	void set_logfile(std::string logfile) {
		m_logfile = logfile;
	}
protected:
	std::string human_readable_extra() const;
private:
	bool boltzmann_compare(double champion_fitness, double challenger_fitness, double temperature) const;
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & const_cast<int &>(m_gen);
		ar & const_cast<double &>(m_cr);
		ar & const_cast<double &>(m_m);
		ar & const_cast<int &>(m_elitism);
		ar & const_cast<mutation &>(m_mut);
		ar & const_cast<size_t &>(m_tournament_size);
		ar & const_cast<crossover::type &>(m_cro);
	}  
	//Number of generations
	const int m_gen;
	//Crossover rate
	const double m_cr;
	//Mutation rate
	const double m_m;

	//Elitism (number of generations after which to reinsert the best)
	const int m_elitism;
	//Mutation
	const mutation m_mut;
	//Selection_type
	const pagmo::population::size_type m_tournament_size;

	selection_operator m_selection_operator;
	double m_initial_temperature;
	const double m_temperature_decrement;

	//Crossover_type
	const crossover::type m_cro;

	std::string m_logfile;
};

//BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::yasga);

#endif // _YASGA_H
