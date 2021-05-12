#include "MCW.hpp"

MCW_model::MCW_model(const symm_matrix<double> &J_, const col_vector<double> &h_, const col_vector<unsigned int> &Numbers_of_spins) : J(J_), h(h_), nodes(col_vector<MCW::node>(Numbers_of_spins.sum())), components(col_vector<col_vector<MCW::node*>>(Numbers_of_spins.size())) {
  unsigned int offset = 0;
  for(unsigned int i=0; i<Numbers_of_spins.size(); ++i) {
    components[i] = col_vector<MCW::node*>(Numbers_of_spins[i]);
    for(unsigned int j=0; j<components[i].size(); ++j) {
       components[i][j] = &nodes[offset+j];
       components[i][j]->component = i;
    }
    offset += Numbers_of_spins[i];
  }
}

MCW_model::MCW_model(const symm_matrix<double> &J_, const col_vector<double> &h_, const col_vector<unsigned int> &Numbers_of_spins, prng& rnd) : J(J_), h(h_), nodes(col_vector<MCW::node>(Numbers_of_spins.sum())), components(col_vector<col_vector<MCW::node*>>(Numbers_of_spins.size())) {
  unsigned int offset = 0;
  for(unsigned int i=0; i<Numbers_of_spins.size(); ++i) {
    components[i] = col_vector<MCW::node*>(Numbers_of_spins[i]);
    for(unsigned int j=0; j<components[i].size(); ++j) {
      components[i][j] = &nodes[offset+j];
      components[i][j]->spin = 2*(rnd()%2)-1;
      components[i][j]->component = i;
    }
    offset += Numbers_of_spins[i];
  }
}

MCW_model& MCW_model::set_fields(const col_vector<double> &fields) {
  h = fields;
  return *this;
}
MCW_model& MCW_model::set_couplings(const symm_matrix<double> &couplings) {
  J = couplings;
  return *this;
}

MCW_model& MCW_model::set_magnetisations(const col_vector<double> &magnetisations, prng& rnd) {
  if(magnetisations.size() != components.size()) {
    std::cerr << "-----------------------------------------------------------------------------------------\n"
              << "set_magnetisations ERROR: Size of magnetisations vector should match number of components\n"
              << "-----------------------------------------------------------------------------------------\n";
    exit(1);
  }
  
  col_vector<double> mags = this->magnetisations();
  for(unsigned int j=0; j<components.size(); ++j) {
    double m = mags[j];
    if(magnetisations[j]>1 || magnetisations[j]<-1) {
      std::cerr << "--------------------------------------------------------------\n"
                << "set_magnetisations ERROR: Magnetisations out of [-1,1] range! \n"
                << "--------------------------------------------------------------\n";
      exit(1);
    }  

    unsigned int sz = components[j].size();
    double dm = 2./sz;
    if(rnd()%2) { // Preventing bias
      while(m < magnetisations[j]-.1*dm) {
        unsigned int i = rnd()%sz; // Selecting random node from component
        if(components[j][i]->spin==-1) {
          components[j][i]->spin = 1;
          m+=dm;
        }
      }
      while(m > magnetisations[j]+.1*dm) {
        unsigned int i = rnd()%sz; // Selecting random node from component
        if(components[j][i]->spin==1) {
          components[j][i]->spin = -1;
          m-=dm;
        }
      }
    }
    else {
      while(m > magnetisations[j]+.1*dm) {
        unsigned int i = rnd()%sz; // Selecting random node from component
        if(components[j][i]->spin==1) {
          components[j][i]->spin = -1;
          m-=dm;
        }
      }
      while(m < magnetisations[j]-.1*dm) { 
        unsigned int i = rnd()%sz; // Selecting random node from component
        if(components[j][i]->spin==-1) {
          components[j][i]->spin = 1;
          m+=dm;
        }
      }
    }
  }
  return *this;
}

MCW_model& MCW_model::randomise(prng& rnd) {
  for(unsigned int j=0; j<components.size(); ++j) {
    double P = .5*rnd.rand_max();
    for(unsigned int i=0; i<components[j].size(); ++i)
      if(rnd() < P)
        components[j][i]->spin = 1;
      else
        components[j][i]->spin = -1;
  }
  return *this;
}

col_vector<double> MCW_model::magnetisations() const {
  col_vector<double> M(components.size());
  for(unsigned int i=0; i<components.size(); ++i) {
    M[i] = 0;
    for(unsigned int j=0; j<components[i].size(); ++j)
      M[i] += components[i][j]->spin;
    M[i] /= components[i].size();
  }
  return M;
}

col_vector<int> MCW_model::Magnetisations() const {
  col_vector<int> M(components.size());
  for(unsigned int i=0; i<components.size(); ++i) {
    M[i] = 0;
    for(unsigned int j=0; j<components[i].size(); ++j)
      M[i] += components[i][j]->spin;
  }
  return M;
}

MCW_model& MCW_model::Metropolis_dynamics(const unsigned int &N_iters, prng& rnd, const unsigned int &N_dynamic_nodes_max, const bool &initialise_randomly, const double &temp) {
  unsigned int rand_max = rnd.rand_max();
  if(initialise_randomly)
    randomise(rnd);

  unsigned int N_nodes = nodes.size();

  col_vector<int> delta_M(components.size()); //Total difference in magnetisations of components
  col_vector<int> M = Magnetisations();
  // Metropolis dynamics main loop
  for(unsigned int n=0; n<N_iters; ++n) {
    unsigned int N_flips = rnd()%N_dynamic_nodes_max+1;
    col_vector<MCW::node*> nodes_to_flip(N_flips);
    for(unsigned int l=0; l<N_flips; ++l) {
      // Finding N_flips nodes proposed to be flipped
      nodes_to_flip[l] = &nodes[rnd()%N_nodes];
      for(unsigned int k=0; k<l; ++k)
        if(nodes_to_flip[k] == nodes_to_flip[l])
          --l;
    }
    // Computing delta_M
    delta_M.set_to(0);
    for(unsigned int l=0; l<N_flips; ++l)
      delta_M[nodes_to_flip[l]->component] -= 2*nodes_to_flip[l]->spin;
    // Computing delta_H
    double delta_H=0;
    double delta_H_=0;
    for(unsigned int l=0; l<components.size(); ++l) {
      delta_H -= h[l]*delta_M[l];
      for(unsigned int s=0; s<=l; ++s)
        delta_H_ -= J[l][s]*(M[l]*delta_M[s]+M[s]*delta_M[l]+delta_M[s]*delta_M[l]);
    }
    delta_H += 1/(double)N_nodes*delta_H_;
    // Accepting or rejecting the proposal
    if(rnd() < exp(-1/temp*delta_H)*rand_max) {//Accept the proposal
      for(unsigned int l=0; l<N_flips; ++l)
        nodes_to_flip[l]->flip_spin();
      for(unsigned int l=0; l<components.size(); ++l)
        M[l]+=delta_M[l];
    }
  }
  return *this;
}

MCW_model& MCW_model::single_spin_Metropolis_dynamics(const unsigned int &N_iters, prng& rnd, const bool &initialise_randomly, const double &temp) {
  unsigned int rand_max = rnd.rand_max();
  if(initialise_randomly)
    randomise(rnd);
  // Rescaling J for convenience of computing delta_H
  symm_matrix<double> Jr(J.size());
  for(unsigned int i=0; i<components.size(); ++i)
    for(unsigned int j=0; j<components.size(); ++j)
      i!=j ? Jr[i][j] = J[i][j] : Jr[i][j] = 2*J[i][j];

  unsigned int N_nodes = nodes.size();
  col_vector<int> M = Magnetisations();
  // Running Metropolis dynamics
  for(unsigned int n=0; n<N_iters; ++n) {
    // Choosing a node at random
    MCW::node &node_to_flip = nodes[rnd()%N_nodes];
    // Computing delta_H
    double delta_H=0;
    int delta_M = -2*node_to_flip.spin;
    for(unsigned int s=0; s<components.size(); ++s)
      delta_H -= Jr[node_to_flip.component][s]*M[s];
    delta_H = delta_M/(double)N_nodes*(delta_H - .5*Jr[node_to_flip.component][node_to_flip.component]*delta_M) - h[node_to_flip.component]*delta_M;
    // Accepting or rejecting the proposal
    if(rnd() < exp(-1/temp*delta_H)*rand_max) {//Accept the proposal
      node_to_flip.flip_spin();
      M[node_to_flip.component] += delta_M;
    }
  }
  return *this;
}

std::ostream& operator<<(std::ostream &os, const MCW_model &MCWM) {
  os << "N:\n" << MCWM.nodes.size() << std::endl
     << "J:\n" << MCWM.J << std::endl
     << "h:\n" << MCWM.h << std::endl;
  return os;
}
