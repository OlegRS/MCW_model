#include "MCW_node.hpp"

MCW::node& MCW::node::flip_spin() {
  spin==1 ? spin=-1 : spin=1;
  return *this;
}

std::ostream& operator<<(std::ostream &os, const MCW::node &nd) {
  os << "spin=" << nd.spin << '\t' << "component=" << nd.component;
  return os;
}
