#ifndef __MCW_HPP__
#define __MCW_HPP__

/////////////////////////////////////////////////////////////////////
// MCW model is defined by the Hamiltonian:                        //
// H = -1/N*\sum_{i<=j}N_iN_jJ_{ij}m_im_j - \sum_{i=1}^ph_iN_im_i, //
// so the intracomponent interactions are counted twice, and the   //
// couplings J correspond to \tilde J from the thesis.             //
/////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>

#include "matrices/col_vector.hpp"
#include "matrices/symm_matrix.hpp"
#include "randomisation/prng.hpp"
#include "MCW_node.hpp"

class MCW_model {
  symm_matrix<double> J;
  col_vector<double> h;
  col_vector<MCW::node> nodes;
  col_vector<col_vector<MCW::node*>> components; // Gives random access to individual components
public:
  MCW_model();
  MCW_model(const unsigned int &p);
  MCW_model(const symm_matrix<double> &J, const col_vector<double> &h, const col_vector<unsigned int> &Numbers_of_spins);
  MCW_model(const symm_matrix<double> &J, const col_vector<double> &h, const col_vector<unsigned int> &Numbers_of_spins, prng& rnd);
  // MCW_model(const symm_matrix<double> &J, const col_vector<double> &h, const col_vector<unsigned int> &Numbers_of_spins, const col_vector<double> &magnetisations, prng& rnd);

  MCW_model& set_fields(const col_vector<double> &fields);
  MCW_model& set_couplings(const symm_matrix<double> &couplings);
  MCW_model& set_magnetisations(const col_vector<double> &magnetisations, prng& rnd);
  MCW_model& randomise(prng& rnd);

  // MCW_model& Glauber_dynamics(prng& rnd, const double &Temperature, const unsigned int &N_iters);
  MCW_model& Metropolis_dynamics(const unsigned int &N_iters, prng& rnd, const unsigned int &N_dynamic_nodes_max=1, const bool &initialise_randomly=false, const double &temp=1);
  MCW_model& single_spin_Metropolis_dynamics(const unsigned int &N_iters, prng& rnd, const bool &initialise_randomly=false, const double &temp=1);

  col_vector<double> magnetisations() const;
  col_vector<int> Magnetisations() const;

  // void save(const std::string &file_name) const;
  
  friend std::ostream& operator<<(std::ostream&, const MCW_model&);
};

#include "matrices/matrix.tpp"
#include "matrices/symm_matrix.tpp"
#include "matrices/col_vector.tpp"

#endif
