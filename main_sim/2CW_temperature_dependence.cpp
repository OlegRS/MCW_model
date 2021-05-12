//////////////////////////////////////////////////////////////////////
// This code performs repetitive sampling for various parameters    //
// of MCW model defined by the Hamiltonian:                         //
// H = -1/N*\sum_{i<=j}N_iN_jJ_{ij}m_im_j - \sum_{i=1}^{p}h_iN_im_i //
// Note that the couplings J in the Hamiltonian correspond to       //
// \tilde J from the thesis, so they are rescaled.                  //
//////////////////////////////////////////////////////////////////////

#include "../MCW_lib/MCW.hpp"
using namespace std;
#define J11 3
#define J22 3
#define J12 -.42
#define N1 8000
#define N2 10000
#define N_NODES (N1+N2)
#define T_STEP .001
#define T_START T_STEP
#define T_FIN 2.5
#define H2 0
#define H1 0
#define N_AVRGING 10
#define N_ITERS_PER_NODE 500
#define N_ITERS (N_NODES*N_ITERS_PER_NODE)
#define N_DYNAMIC_NODES 1
#define INITIALISE_RANDOMLY false
#define INIT_m1 -1
#define INIT_m2 1

int main() {
  // Rescale the couplings to compensate for double
  // counting of intracomponent interactions as follows:
  // J[i][j]=.5*J_[i][j],  i==j
  // J[i][j]=   J_[i][j],  i!=j
  symm_matrix<double> J(2);
  J[0][0]=J11*.5; J[0][1]=J12;
  J[1][0]=J12;    J[1][1]=J22*.5;
  col_vector<double> h(2);
  h[0] = H1; h[1] = H2;
  col_vector<unsigned int> N(2);
  N[0] = N1;   N[1] = N2;
  col_vector<double> init_m(2);
  init_m[0]=INIT_m1; init_m[1]=INIT_m2;

  prng rnd(RAND_MAX);

  MCW_model MCW(J, h, N, rnd);
  MCW.set_magnetisations(init_m, rnd);
  cout << MCW << '\n';
  cout << MCW.magnetisations() << '\n';

  cout << "../data/2CW_model/T_scanning/2CW__J11" + to_string(J11) + "__J22_" + to_string(J22) + "__J12_" + to_string(J12) + "__T_START_" + to_string(T_START) + "__T_FIN_" + to_string(T_FIN) + "__T_STEP_"  + to_string(T_STEP) + "__H1_" + to_string(H1) + "__H2_" + to_string(H2) + "__N1_" + to_string(N1)+ "__N2_" + to_string(N2)+ "__INIT_m1_" + to_string(INIT_m1) + "__INIT_m2_" + to_string(INIT_m2) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv\n";

  ofstream ofs("../data/2CW_model/T_scanning/2CW__J11" + to_string(J11) + "__J22_" + to_string(J22) + "__J12_" + to_string(J12) + "__T_START_" + to_string(T_START) + "__T_FIN_" + to_string(T_FIN) + "__T_STEP_"  + to_string(T_STEP) + "__H1_" + to_string(H1) + "__H2_" + to_string(H2) + "__N1_" + to_string(N1)+ "__N2_" + to_string(N2)+ "__INIT_m1_" + to_string(INIT_m1) + "__INIT_m2_" + to_string(INIT_m2) + "__ITERS_PER_NODE_" + to_string(N_ITERS_PER_NODE) + "__N_AVRGING_" + to_string(N_AVRGING) + "__RANDOM_INITIALISATION_" + to_string(INITIALISE_RANDOMLY) + "__N_DYNAMIC_NODES_" + to_string(N_DYNAMIC_NODES) +".csv");
  
  if(!ofs.is_open()) {
    cerr << "ERROR openning file to save the data!\n" << endl;
    exit(1);
  }

  col_vector<double> magnetisations(2);
  for(double T=T_START; T<T_FIN; T+=T_STEP) {
    cerr << T << '\t';
    MCW.set_fields(h);
    for(unsigned int i=0; i<N_AVRGING; ++i) {
      // MCW.Metropolis_dynamics(N_ITERS, rnd, N_DYNAMIC_NODES);
      MCW.single_spin_Metropolis_dynamics(N_ITERS, rnd, false, T);
      magnetisations = MCW.magnetisations();
      ofs << T << ',' << magnetisations[0] << ',' << magnetisations[1] << '\n';
      ofs.flush();
    }
  }  
  return 0;
}
