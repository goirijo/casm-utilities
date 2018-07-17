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
	    "substructures,s",
	    po::value<std::vector<fs::path>>()->multitoken()->required(),
	    "POS.vasp like files you want to "
	    "stack on top of each other.");
	frankenstack_desc.add_options()("number,n", po::value<int>(),
					"number of times to repeat the "
					"stacking if only one unit is given");
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

	auto sub_paths =
	    frankenstack_launch.fetch<std::vector<fs::path>>("substructures");

	if (!sub_paths.size()) {
		std::cerr << "Missing substructure path " << std::endl;
	}
	auto big_struc = Rewrap::Structure(sub_paths[0]);
	if (sub_paths.size() == 1) {
		auto unit = Rewrap::Structure(sub_paths[0]);
		int num_stacks = 1;
		if (frankenstack_launch.vm().count("number")) {
			num_stacks = frankenstack_launch.fetch<int>("number");
		}
		std::vector<Rewrap::Structure> struc_vec;
		struc_vec.insert(struc_vec.end(), num_stacks, unit);
		big_struc = Frankenstein::structure_stacker(struc_vec);

	} else {
		std::vector<Rewrap::Structure> units;
		for (auto& item : sub_paths) {
			std::cout << item << std::endl;
			units.push_back(Rewrap::Structure(item));
		}
		big_struc = Frankenstein::structure_stacker(units);
	}
	if (frankenstack_launch.vm().count("output")) {
		Simplicity::write_poscar(
		    big_struc, frankenstack_launch.fetch<fs::path>("output"));
		return 0;
	}
	Simplicity::print_poscar(big_struc, std::cout);
	return 0;
}
