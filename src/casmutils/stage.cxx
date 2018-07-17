#include "casmutils/stage.hpp"
#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include "casmutils/structure.hpp"

namespace Frankenstein {
/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction) and splits it into two structures
/// slice_loc is in fractional length of the c vector
std::pair<Rewrap::Structure, Rewrap::Structure> structure_slicer(
    const Rewrap::Structure &big_struc, double slice_loc, double tol) {
	Rewrap::Structure cpy_big = big_struc;
	Simplicity::mod_coordinates(&cpy_big);
	std::vector<Rewrap::Structure> struc_array;
	Eigen::Matrix3d lat_mat = big_struc.lattice().lat_column_mat();
	lat_mat.col(2) = lat_mat.col(2) * (slice_loc);
	CASM::Lattice bottom_lat(lat_mat);
	CASM::Structure bottom_struc(bottom_lat);
	auto b_struc = Rewrap::Structure(bottom_struc);
	lat_mat = big_struc.lattice().lat_column_mat();
	lat_mat.col(2) = lat_mat.col(2) * (1 - slice_loc);
	CASM::Lattice top_lat(lat_mat);
	CASM::Structure top_struc(top_lat);
	auto t_struc=Rewrap::Structure(top_struc);
	for (const auto &item : cpy_big.basis) {
		/// only move basis sites below slice pivot
		if (item.const_frac()(2) >= 0 &&
		    item.const_frac()(2) < slice_loc + tol) {
			auto coord = CASM::Coordinate(item.const_cart(),
						      bottom_lat, CASM::CART);
			auto site = CASM::Site(coord, item.occ_name());
			bottom_struc.basis.push_back(site);
			b_struc=Rewrap::Structure(bottom_struc);
			Simplicity::mod_coordinates(&b_struc);
		}
		else {
			CASM::Site new_site = item;
			/// adjust c coord by slice location
			Eigen::Vector3d altered = new_site.const_frac();
			altered(2) = new_site.const_frac()(2) - slice_loc;
			new_site.frac() = altered;
			auto coord = CASM::Coordinate(new_site.const_cart(),
						      top_lat, CASM::CART);
			auto site = CASM::Site(coord, item.occ_name());
			top_struc.basis.push_back(site);
			t_struc=Rewrap::Structure(top_struc);
			Simplicity::mod_coordinates(&t_struc);

		}
	}
	return std::make_pair(b_struc, t_struc);
}

///
/// This function takes a structure ( where c is perpendicular to ab plane and
/// is layering direction)
/// and splits it into n+1 structures where n is the size of slice_loc and the
/// entries in slice_loc
/// dictate where the slicing happens in increasing order. Entries in slice_loc
/// are fractional lengths of the c vector
std::vector<Rewrap::Structure> multi_slice(const Rewrap::Structure &big_struc,
					   const Eigen::VectorXd &slice_loc,
					   double tol) {
	auto struc_pair = structure_slicer(big_struc, slice_loc(0), tol);
	std::vector<Rewrap::Structure> slices;
	slices.push_back(struc_pair.first);
	if (slice_loc.size() == 1) {
		slices.push_back(struc_pair.second);
		return slices;
	}
	Eigen::VectorXd tmp =
	    Eigen::VectorXd::Constant(slice_loc.size() - 1, slice_loc(0));
	tmp = slice_loc.tail(slice_loc.size() - 1) - tmp;
	auto after_piv =
	    multi_slice(struc_pair.second, tmp / (1.0 - slice_loc(0)), tol);
	slices.insert(slices.end(), after_piv.begin(), after_piv.end());
	return slices;
}
/// This function splits a structure in equally sized slices along the c -axis
/// (layering)
/// which is perpendicular to the ab plane.
std::vector<Rewrap::Structure> equal_slice(const Rewrap::Structure &big_struc,
					   int n_pieces) {
	std::vector<double> slice_locs;
	for (int i = 1; i < n_pieces; i++) {
		slice_locs.push_back((double)i / (double)n_pieces);
	}
	Eigen::Map<Eigen::VectorXd> eigen_locs(&slice_locs[0],
					       slice_locs.size());
	return multi_slice(big_struc, eigen_locs, CASM::TOL);
}

/// This function takes a vector of structures with the same ab lattice vectors
/// and stacks them along the c direction.
Rewrap::Structure structure_stacker(
    const std::vector<Rewrap::Structure> &sub_strucs) {
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
		if (i > 0 ){
		c_shift =
		    c_shift + sub_strucs[i].lattice().lat_column_mat().col(2);
		}
		Rewrap::Structure cpy_i=sub_strucs[i];
		Simplicity::mod_coordinates(&cpy_i);
		for (const auto &item : cpy_i.basis) {
			CASM::Site new_site = item;
			/// adjust site by a, b, and c shifts
			Eigen::Vector3d altered = new_site.const_cart();
			altered = new_site.const_cart() + c_shift;
			new_site.cart() = altered;
			stacked_struc.basis.push_back(CASM::Site(
			    CASM::Coordinate(new_site.const_cart(), stacked_lat,
					     CASM::CART),
			    item.occ_name()));
		}
	}
	auto rw_struc=Rewrap::Structure(stacked_struc);
	Simplicity::mod_coordinates(&rw_struc);
	return rw_struc;
}
}
/// This function takes a structures and shifts the origin by shift val
/// shift val is in fractional coordinates of the lattice
Rewrap::Structure *origin_shift(Rewrap::Structure *struc,
			       const Eigen::Vector3d &shift_val) {
	for (auto &item : struc->basis) {
		item +=
		    CASM::Coordinate(shift_val, struc->lattice(), CASM::FRAC);
	}
	Simplicity::mod_coordinates(struc);
	return struc;
}

/// This function alters the coordinates of the given struc to have fractional
/// coordinates between 0 and 1 only
Rewrap::Structure *Simplicity::mod_coordinates(Rewrap::Structure *struc) {
	for (auto &site : struc->basis) {
		site.within();
	}
	return struc;
}
