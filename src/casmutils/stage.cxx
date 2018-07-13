#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <boost/filesystem.hpp>
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <fstream>



///
/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction)
/// and splits it into n+1 structures where n is the size of slice_loc and the
/// entries in slice_loc
/// dictate where the slicing happens in increasing order. Entries in slice_loc
/// are fractional lengths of the c vector
std::vector<Rewrap::Structure> structure_slicer(const Rewrap::Structure &big_struc,
					std::vector<double> &slice_loc) {
	double lower_bound = 0.0;
	double upper_bound = 0.0;
	auto it = slice_loc.begin();
	std::vector<Rewrap::Structure> struc_array;
	for (; it != slice_loc.end(); it++) {
		upper_bound = *it;
		Eigen::Matrix3d lat_mat = big_struc.lattice().lat_column_mat();
		lat_mat.col(2) = lat_mat.col(2) * (upper_bound - lower_bound);
		CASM::Lattice sliced_lat(lat_mat);
		CASM::Structure sliced_struc(sliced_lat);
		for (const auto &item : big_struc.basis) {
			/// only move basis sites below slice pivot
			if (item.const_frac()(2) > lower_bound &&
			    item.const_frac()(2) < upper_bound) {
				CASM::Site new_site = item;
				/// adjust c coord by lower bound
				Eigen::Vector3d altered = new_site.const_frac();
				altered(2) =
				    new_site.const_frac()(2) - lower_bound;
				new_site.frac() = altered;
				sliced_struc.basis.push_back(
				    CASM::Site(CASM::Coordinate(new_site.const_cart(),
						    sliced_lat, CASM::CART),
					 item.occ_name()));
			}
		}
		struc_array.push_back(sliced_struc);
		lower_bound = upper_bound;
		if (std::distance(it, slice_loc.end()) == 1) {
			upper_bound = 1;
			Eigen::Matrix3d lat_mat =
			    big_struc.lattice().lat_column_mat();
			lat_mat.col(2) =
			    lat_mat.col(2) * (upper_bound - lower_bound);
			CASM::Lattice sliced_lat(lat_mat);
			CASM::Structure sliced_struc(sliced_lat);
			for (auto &item : big_struc.basis) {
				/// only move basis sites below slice pivot
				if (item.const_frac()(2) >= lower_bound &&
				    item.const_frac()(2) < upper_bound) {
					CASM::Site new_site = item;
					/// adjust c coord by lower bound
					Eigen::Vector3d altered =
					    new_site.const_frac();
					altered(2) = new_site.const_frac()(2) -
						     lower_bound;
					new_site.frac() = altered;
					sliced_struc.basis.push_back(CASM::Site(
					    CASM::Coordinate(new_site.const_cart(),
						       sliced_lat, CASM::CART),
					    item.occ_name()));
				}
			}
			struc_array.push_back(Rewrap::Structure(sliced_struc));
		}
	}
	return struc_array;
}
/// This function splits a structure in equally sized slices along the c -axis
/// (layering)
/// which is perpendicular to the ab plane.
std::vector<Rewrap::Structure> structure_slicer(const Rewrap::Structure &big_struc,
					int n_pieces) {
	std::vector<double> slice_locs;
	for (int i = 1; i < n_pieces; i++) {
		slice_locs.push_back((double)i / (double)n_pieces - 0.0001);
	}
	return structure_slicer(big_struc, slice_locs);
}

/// This function takes a vector of structures with the same ab lattice vectors
/// and stacks them along the c direction. A corresponding vector of ab in plane
/// shift values relative to the first structure is given
Rewrap::Structure structure_stacker(
    std::vector<Rewrap::Structure> &sub_strucs,
    std::vector<std::pair<double, double>> &shift_vals) {
	Eigen::Matrix3d stacked_lat_mat =
	    sub_strucs[0].lattice().lat_column_mat();
	// extend c axis accordingly
	for (int i = 1; i < sub_strucs.size(); i++) {
		Eigen::Matrix3d lat_mat =
		    sub_strucs[i].lattice().lat_column_mat();
		stacked_lat_mat.col(2) =
		    stacked_lat_mat.col(2) + lat_mat.col(2);
	}
	CASM::Lattice stacked_lat(stacked_lat_mat);
	CASM::Structure stacked_struc(stacked_lat);
	Eigen::Vector3d c_shift = Eigen::Vector3d::Zero();
	for (int i = 0; i < sub_strucs.size(); i++) {
		// determine appropriate c-axis shift for position in stacking
		c_shift =
		    c_shift + sub_strucs[i].lattice().lat_column_mat().col(2);
		Eigen::Vector3d a_shift =
		    shift_vals[i].first *
		    stacked_struc.lattice().lat_column_mat().col(0);
		Eigen::Vector3d b_shift =
		    shift_vals[i].second *
		    stacked_struc.lattice().lat_column_mat().col(1);
		for (const auto &item : sub_strucs[i].basis) {
			CASM::Site new_site = item;
			/// adjust site by a, b, and c shifts
			Eigen::Vector3d altered = new_site.const_cart();
			altered =
			    new_site.const_cart() + a_shift + b_shift + c_shift;
			new_site.cart() = altered;
			stacked_struc.basis.push_back(
			    CASM::Site(CASM::Coordinate(new_site.const_cart(), stacked_lat,
					    CASM::CART),
				 item.occ_name()));
		}
	}
	return stacked_struc;
}

/// This version of structure sets the same horizontal shift between each
/// consecutive layer
/// Layers are possibly different
Rewrap::Structure structure_stacker(std::vector<Rewrap::Structure> &sub_strucs,
			    const std::pair<double, double> &shift_value) {
	std::vector<std::pair<double, double>> shift_vals;
	for (int i = 0; i < sub_strucs.size(); i++) {
		auto pair =std::make_pair(i * shift_value.first,
						    i * shift_value.second);
		shift_vals.push_back(pair);
	}
	return structure_stacker(sub_strucs, shift_vals);
}

/// This version of structure sets the same layer as the stacking unit
/// Shifts are possibly different
Rewrap::Structure structure_stacker(
    Rewrap::Structure &unit, std::vector<std::pair<double, double>> &shift_vals) {
	std::vector<Rewrap::Structure> sub_strucs;
	sub_strucs.insert(sub_strucs.end(), shift_vals.size(), unit);
	return structure_stacker(sub_strucs, shift_vals);
}

/// This version of structure sets the same layer as the stacking unit and the
/// same shift between each layer n times
Rewrap::Structure structure_stacker(Rewrap::Structure &unit,
			    const std::pair<double, double> &shift_value,
			    int n_times) {
	std::vector<Rewrap::Structure> sub_strucs;
	sub_strucs.insert(sub_strucs.end(), n_times, unit);
	return structure_stacker(sub_strucs, shift_value);
}
