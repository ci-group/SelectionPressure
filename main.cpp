/*****************************************************************************
 *   Copyright (C) 2014 VU University Amsterdam
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the
 *   Free Software Foundation, Inc.,
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *   $Id: main.cpp 161 2015-01-19 07:56:05Z ehaasdi $
 *
 *****************************************************************************/
#include <iostream>
#include <iomanip>
#include <pagmo/pagmo.h>
#include <pagmo/algorithm/de.h>
#include "yasga.h"
#include "my_nsga2.h"

#include "version.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int ac, char** av)
{
	// Declare the supported options.
	po::options_description desc;
	desc.add_options()
	    ("help", "produce help message")
	    ("s", po::value<unsigned>(), "selection scheme (0,1,2) for (tournament, boltzmann, roulettewheel)")
	    ("t", po::value<unsigned>(), "tournament size")
	    ("p", po::value<unsigned>(), "problem (0,1,2,3) for (Himmelblau, Schwefel, Fonseca and Fleming, Schaffer's study)")
	    ("l", po::value<std::string>(), "logfile")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	std::string opt = "help";
	if (vm.count(opt)) {
	    std::cout << desc << "\n";
	    return 1;
	}

	unsigned tournament(2);
	opt = "t";
	if (vm.count(opt)) {
	    tournament = vm[opt].as<unsigned>();
	}

	unsigned p(0);
	opt = "p";
	if (vm.count(opt)) {
		p = vm[opt].as<unsigned>();
	}

	pagmo::problem::base_ptr problem;
	double xp(0.7), mp(1.0), mw(0.6); // settings from rapid tuning run for himmelblau
	switch (p) {
		case 0 : {
			problem.reset(new pagmo::problem::himmelblau);
			break;
		}
		case 1: {
			problem.reset(new pagmo::problem::schwefel(10));
			xp = 1; mp = .25;  mw = 1; 	// settings from rapid tuning run for schwefel
			break;
		}
		case 2: {
			problem.reset(new pagmo::problem::fon);
		}
		case 3: {
			problem.reset(new pagmo::problem::sch);
		}
	}

	pagmo::population pop(*problem, 200);

	if (problem->get_f_dimension() == 1) {

		unsigned s(0);
		opt = "s";
		if (vm.count(opt)) s = vm[opt].as<unsigned>();

		yasga::selection_operator sel( static_cast<yasga::selection_operator>(s) );

		yasga yasga_optimiser(50000, xp, mp, 1, yasga::mutation::GAUSSIAN,  mw, tournament, sel, 1000, 0.02);

		opt = "l";
		if (vm.count(opt)) {
			yasga_optimiser.set_logfile(vm[opt].as<std::string>());
		}

		yasga_optimiser.evolve(pop);
	} else {
		double m_prob = 1.0 / problem->get_ub().size();
		pagmo::algorithm::my_nsga2 nsga(50000, 0.9, 20, m_prob, 20);
		// settings from Deb, K. and Pratap, A. and Agarwal, S. and Meyarivan, T., "A fast and elitist multiobjective genetic algorithm: NSGA-II"
		opt = "l";
		if (vm.count(opt)) {
			nsga.set_logfile(vm[opt].as<std::string>());
		}

		nsga.evolve(pop);
	}

	std::cout << pop.champion().f << std::endl;

	return 0;
}
