#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>
#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"

namespace Utilities {

void originshift_initializer(po::options_description& originshift_desc) {
	UtilityProgramOptions::add_help_suboption(originshift_desc);
	UtilityProgramOptions::add_desc_suboption(originshift_desc);
	UtilityProgramOptions::add_output_suboption(originshift_desc);
	originshift_desc.add_options()("structure",
				       po::value<fs::path>()->required(),
				       "POS.vasp like file you want to "
				       "shift all the coordinates of");
	originshift_desc.add_options()(
	    "shift", po::value<std::vector<double>>()->multitoken()->required(),
	    "shift value that will be added to all coordinates (negative "
	    "origin shift) units are fractional");
	return;
}
}

using namespace Utilities;

int main(int argc, char* argv[]) {
	Handler originshift_launch(argc, argv, originshift_initializer);

	if (originshift_launch.count("help")) {
		std::cout << originshift_launch.desc() << std::endl;
		return 1;
	}

	try {
		originshift_launch.notify();
	}

	catch (po::required_option& e) {
		std::cerr << e.what() << std::endl;
		return 2;
	}

	auto struc_path = originshift_launch.fetch<fs::path>("structure");
	auto struc = Rewrap::Structure(struc_path);
	auto out_struc = struc;
	auto vec = originshift_launch.fetch<std::vector<double>>("shift");
	origin_shift(&out_struc, Eigen::Map<Eigen::Vector3d>(&vec[0]));
	if (originshift_launch.vm().count("output")) {
		Simplicity::write_poscar(
		    out_struc, originshift_launch.fetch<fs::path>("output"));
		return 0;
	}
	Simplicity::print_poscar(out_struc, std::cout);
	return 0;
}
