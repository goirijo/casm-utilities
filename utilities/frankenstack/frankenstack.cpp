#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>
#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"

namespace Utilities {

void frankenstack_initializer(po::options_description& frankenstack_desc) {
	UtilityProgramOptions::add_help_suboption(frankenstack_desc);
	UtilityProgramOptions::add_desc_suboption(frankenstack_desc);
	UtilityProgramOptions::add_output_suboption(frankenstack_desc);
	frankenstack_desc.add_options()(
	    "superstructures,s",
	    po::value<std::vector<fs::path>>()->multitoken()->required(),
	    "POS.vasp like files you want to "
	    "stack on top of each other.");
	frankenstack_desc.add_options()(
	    "vector,v",
	    po::value<std::vector<double>>()->multitoken()->required(),
	    "a vector of ab shift values in fractional lengths");

	frankenstack_desc.add_options()("number,n", po::value<int>(),
					"number of times to repeat the "
					"stacking if only one unit and shift "
					"value is given");
	return;
}
}

using namespace Utilities;

int main(int argc, char* argv[]) {
	Handler frankenstack_launch(argc, argv, frankenstack_initializer);

	if (frankenstack_launch.count("help")) {
		std::cout << frankenstack_launch.desc() << std::endl;
		return 1;
	}

	try {
		frankenstack_launch.notify();
	}

	catch (po::required_option& e) {
		std::cerr << e.what() << std::endl;
		return 2;
	}

	auto super_paths =
	    frankenstack_launch.fetch<std::vector<fs::path>>("superstructures");
	auto unrolled_shifts =
	    frankenstack_launch.fetch<std::vector<double>>("vector");
	if (unrolled_shifts.size() % 2 == 1) {
		std::cerr
		    << "Odd amount of numbers need pairs of a and b shifts"
		    << std::endl;
		return 3;
	}
	std::vector<std::pair<double, double>> shift_values;
	for (int i = 0; i < unrolled_shifts.size() / 2; i++) {
		shift_values.push_back(std::make_pair(
		    unrolled_shifts[2 * i], unrolled_shifts[2 * i + 1]));
	}
	if (!super_paths.size() || !shift_values.size()) {
		std::cerr << "Missing superstructure path or shift vector"
			  << std::endl;
	}
	if (super_paths.size() == 1) {
		auto unit = Rewrap::Structure(super_paths[0]);
		if (shift_values.size() == 1) {
			auto shift_value = shift_values[0];
			int num_stacks = 1;
			if (frankenstack_launch.vm().count("number")) {
				num_stacks =
				    frankenstack_launch.fetch<int>("number");
			}
			std::cout << num_stacks << std::endl;
			auto big_struc =
			    structure_stacker(unit, shift_value, num_stacks);
			if (frankenstack_launch.vm().count("output")) {
				Simplicity::write_poscar(
				    big_struc,
				    frankenstack_launch.fetch<fs::path>(
					"output"));
				return 0;
			}
			Simplicity::print_poscar(big_struc, std::cout);
			return 0;
		}
		auto big_struc = structure_stacker(unit, shift_values);
		if (frankenstack_launch.vm().count("output")) {
			Simplicity::write_poscar(
			    big_struc,
			    frankenstack_launch.fetch<fs::path>("output"));
			return 0;
		}
		Simplicity::print_poscar(big_struc, std::cout);
		return 0;

	} else {
		std::vector<Rewrap::Structure> units;
		for (auto& item : super_paths) {
			units.push_back(Rewrap::Structure(item));
		}
		if (shift_values.size() == 1) {
			auto shift_value = shift_values[0];
			auto big_struc = structure_stacker(units, shift_value);
			if (frankenstack_launch.vm().count("output")) {
				Simplicity::write_poscar(
				    big_struc,
				    frankenstack_launch.fetch<fs::path>(
					"output"));
				return 0;
			}
			Simplicity::print_poscar(big_struc, std::cout);
			return 0;
		} else {
			auto big_struc = structure_stacker(units, shift_values);
			if (frankenstack_launch.vm().count("output")) {
				Simplicity::write_poscar(
				    big_struc,
				    frankenstack_launch.fetch<fs::path>(
					"output"));
				return 0;
			}
			Simplicity::print_poscar(big_struc, std::cout);
			return 0;
		}
	}
}
