#ifndef __MCW_NODE_HPP__
#define __MCW_NODE_HPP__
#include<iostream>

namespace MCW {
  struct node {
    short spin;
    unsigned int component;

    node& flip_spin();
    
    friend std::ostream& operator<<(std::ostream&, const node&);
  };
}

#endif
