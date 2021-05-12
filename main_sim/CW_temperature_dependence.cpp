//////////////////////////////////////////////////////////////////////
// This code performs repetitive sampling for various parameters    //
// of MCW model defined by the Hamiltonian:                         //
// H = -1/N*\sum_{i<=j}N_iN_jJ_{ij}m_im_j - \sum_{i=1}^{p}h_iN_im_i //
// Note that the couplings J in the Hamiltonian correspond to       //
// \tilde J from the thesis, so they are rescaled.                  //
//////////////////////////////////////////////////////////////////////

#include "../MCW_lib/MCW.hpp"
using namespace std;
#define J_ 1
#define N_NODES 10000
#define T_STEP .001
#define T_START T_STEP
#define T_FIN 2
#define H .01
#define N_AVRGING 10
#define N_ITERS_PER_NODE 10
#define N_ITERS (N_NODES*N_ITERS_PER_NODE)
#define N_DYNAMIC_NODES 1
#define INITIALISE_RANDOMLY false
#define INIT_m -1

int main() {
  // Rescale the couplings to compensate for double
  // counting of intracomponent interactions as follows:
  // J[i][j]=.5*J_[i][j],  i==j
  // J[i][j]=   J_[i][j],  i!=j
  symm_matrix<double> J(1);
  J[0][0]=J_*.5;
  col_vector<double> h(1);
  h[0] = H;
  col_vector<unsigned int> N(1);
  N[0] = N_NODES;
  col_vector<double> init_m(1);
  init_m[0]=INIT_m;

  prng rnd(RAND_MAX);

  MCW_model MCW(J, h, N, rnd);
  MCW.set_magnetisations(init_m, rnd);
  cout << MCW << '\n';
  cout << MCW.magnetisations() << '\n';

  cout << "../data/CW_model/T_scanning/CW__J_" + to_string(J_) + "__H_" + to_string(H) + "__INIT_m_" + to_string(INIT_m) + "__T_START_" + to_string(T_START) + "__T_FIN_" + to_string(T_FIN) + "__T_STEP_"  + to_string(T_STEP) + "__N_NODES_" + to_string(N_NODES) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv\n";
  ofstream ofs("../data/CW_model/T_scanning/CW__J_" + to_string(J_) + "__H_" + to_string(H) + "__INIT_m_" + to_string(INIT_m) + "__T_START_" + to_string(T_START) + "__T_FIN_" + to_string(T_FIN) + "__T_STEP_"  + to_string(T_STEP) + "__N_NODES_" + to_string(N_NODES) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv");
  if(!ofs.is_open()) {
    cerr << "ERROR openning file to save the data!\n" << endl;
    exit(1);
  }
  
  if(!ofs.is_open()) {
    cerr << "ERROR openning file to save the data!\n" << endl;
    exit(1);
  }

  for(double T=T_START; T<T_FIN; T+=T_STEP) {
    cerr << T << '\t';
    MCW.set_fields(h);
    for(unsigned int i=0; i<N_AVRGING; ++i) {
      // MCW.Metropolis_dynamics(N_ITERS, rnd, N_DYNAMIC_NODES);
      MCW.single_spin_Metropolis_dynamics(N_ITERS, rnd, false, T);
      ofs << T << ',' << MCW.magnetisations()[0] << '\n';
      ofs.flush();
    }
  }  
  return 0;
}
