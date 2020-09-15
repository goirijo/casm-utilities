#include "casm/crystallography/io/VaspIO.hh"

#include <iomanip>
#include "casm/container/algorithm.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
  namespace VaspIO {

    // --- Definitions -------------------------------------------------------- //


    /// \brief Construct PrintPOSCAR object
    ///
    /// By default:
    /// - title = config.name()
    /// - scale = 1.0
    /// - coordinate mode = frac (Direct)
    /// - atom names line is printed
    /// - Atoms printed in order appearing in the Configuration. (No sorting by atom types)
    /// - No selective dynamics
    /// - {"Va", "va", "VA"} atoms not printed
    /// - also adds displacements and deformation if they are present as DoF
    ///
    /// Currently:
    /// - assumes molecules should be printed as their individual atoms.
    ///
    PrintPOSCAR::PrintPOSCAR(xtal::SimpleStructure _struc, std::string _title, PrintPOSCAR::SpeciesMode _mode) :
      m_title(std::move(_title)),
      m_species_mode(_mode),
      m_struc(std::move(_struc)),
      m_permute(sequence<Index>(0, m_struc.info(_mode).size() - 1)),
      m_scale(1.0),
      m_coord_mode(FRAC),
      m_atom_names(true),
      m_sel_dynamics(false),
      m_append_atom_names(true),
      m_ignore {"VA", "Va", "va"} {
      m_sel_dynamics = _struc.info(m_species_mode).properties.count(AnisoValTraits::selective_dynamics().name());
    }


    /// \brief Default sort is by species name
    void PrintPOSCAR::sort() {
      m_permute = m_struc.info(m_species_mode).sort_by_name();
    }

    /// \brief Print POSCAR, providing a range of std::tuple<AtomName, Coordinate, SelectiveDynamics>
    void PrintPOSCAR::print(std::ostream &sout) const {

      int tprec = sout.precision();
      std::ios::fmtflags tflags = sout.flags();
      int prec = 8;
      sout.precision(prec);

      auto const &info = m_struc.info(m_species_mode);

      if(m_permute.size() != info.size()) {
        throw std::runtime_error("Error in PrintPOSCAR::print: m_permute.size() != info.size()");
      }

      // first filter out all atoms we are going to ignore, indices of the remaining atoms get put in 'atom'
      std::vector<Index> atom;

      // Count size of each contiguous block of non-ignored species
      std::vector<std::pair<std::string, Index> > atom_count;

      for(Index i : m_permute) {
        // if Atom's name is not found in the ignore list, add it to 'atom'
        if(ignore().count(info.names[i]) == 0) {
          atom.push_back(i);
          if(atom_count.empty() || info.names[i] != atom_count.back().first)
            atom_count.emplace_back(info.names[i], 0);

          atom_count.back().second++;
        }
      }

      // Initialize coordinate transformation matrix
      Eigen::Matrix3d c2f;
      if(m_coord_mode == CART)
        c2f.setIdentity();
      else
        c2f = m_struc.lat_column_mat.inverse();

      // Determine printing width and construct Eigen::IOFormat object
      int width = 12;
      width = print_matrix_width(sout, m_struc.lat_column_mat.transpose(), width);
      for(Index i : atom) {
        Eigen::Vector3d vec = c2f * info.cart_coord(i);
        width = print_matrix_width(sout, vec.transpose(), width);
      }
      Eigen::IOFormat format(prec, width + 1);

      // Begin printing:
      // print title, scale, and lattice
      sout <<  m_title << "\n";
      sout <<  std::fixed << std::setprecision(prec) << m_scale << "\n";

      sout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

      for(Index i = 0; i < 3; ++i)
        sout <<  m_struc.lat_column_mat.col(i).transpose().format(format) << '\n';

      // if after filtering out ignored atoms none are left, return
      if(atom.size() == 0) {
        return;
      }

      // optionally print atom names line

      if(m_atom_names) {
        for(auto const &p : atom_count) {
          sout << p.first << " ";
        }
        sout << "\n";
      }

      for(auto const &p : atom_count) {
        sout << p.second << " ";
      }
      sout << "\n";

      // print 'Selective Dynamics' if using selective dynamics
      if(m_sel_dynamics) {
        sout << "Selective Dynamics\n";
      }

      // print coord mode
      sout <<  xtal::COORD_MODE::NAME(m_coord_mode) << "\n";

      // print all coordinates, and seletive dynamics settings, and atom names if applicable
      for(Index i : atom) {
        sout <<  (c2f * info.cart_coord(i)).transpose().format(format);


        if(m_sel_dynamics) {
          Eigen::Vector3i sd = iround(info.properties.at(AnisoValTraits::selective_dynamics().name()).col(i));
          sout << " ";
          for(Index j = 0; j < 3; ++j)
            sout << (sd[j] ? "T " : "F ");
        }

        if(m_append_atom_names) {
          sout << " " << info.names[i];
        }
        sout << "\n";
      }
      sout << "\n";

      sout.precision(tprec);
      sout.flags(tflags);
    }
  }
}
