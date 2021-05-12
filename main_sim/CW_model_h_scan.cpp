//////////////////////////////////////////////////////////////////////
// This code performs repetitive sampling for various parameters    //
// of MCW model defined by the Hamiltonian:                         //
// H = -1/N*\sum_{i<=j}N_iN_jJ_{ij}m_im_j - \sum_{i=1}^{p}h_iN_im_i //
// Note that the couplings J in the Hamiltonian correspond to       //
// \tilde J from the thesis, so they are rescaled.                  //
//////////////////////////////////////////////////////////////////////

#include "../MCW_lib/MCW.hpp"
using namespace std;
//J_ are J from the thesis (not \tilde J)
#define J_ 2
#define N_NODES 50
#define H_STEP .001
#define H_START .08
#define H_FIN .1
#define N_AVRGING 10000
#define N_ITERS_PER_NODE 10
#define N_ITERS (N_NODES*N_ITERS_PER_NODE)
#define N_DYNAMIC_NODES N_NODES
#define INITIALISE_RANDOMLY false

int main() {
  // Rescale the couplings to compensate for double
  // counting of intracomponent interactions as follows:
  // J[i][j]=.5*J_[i][j],  i==j
  // J[i][j]=.5*J_[i][j],  i!=j
  symm_matrix<double> J(1);
  J[0][0]=J_*.5;
  col_vector<double> h(1);
  h[0] = 1;
  col_vector<unsigned int> N(1);
  N[0] = N_NODES;

  prng rnd(RAND_MAX);

  MCW_model MCW(J, h, N, rnd);
  MCW.randomise(rnd);

  cout << "../data/CW_model/h_scanning/CW__J_" + to_string(J_) + "__H_START_" + to_string(H_START) + "__H_FIN_" + to_string(H_FIN) + "__H_STEP_"  + to_string(H_STEP) + "__N_NODES_" + to_string(N_NODES) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv\n";
  ofstream ofs("../data/CW_model/h_scanning/CW__J_" + to_string(J_) + "__H_START_" + to_string(H_START) + "__H_FIN_" + to_string(H_FIN) + "__H_STEP_"  + to_string(H_STEP) + "__N_NODES_" + to_string(N_NODES) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv");
  if(!ofs.is_open()) {
    cerr << "ERROR openning file to save the data!\n" << endl;
    exit(1);
  }

  for(h[0]=H_START; h[0]<H_FIN; h[0]+=H_STEP) {
    cerr << h[0] << '\t';
    MCW.set_fields(h);
    for(unsigned int i=0; i<N_AVRGING; ++i) {
      MCW.Metropolis_dynamics(N_ITERS, rnd, N_DYNAMIC_NODES);
      //!!! SET_DYNAMIC_NODES to 1!!! MCW.single_spin_Metropolis_dynamics(N_ITERS, rnd);
      ofs << h[0] << ',' << MCW.magnetisations()[0] << '\n';
      ofs.flush();
    }
  }  
  return 0;
}
