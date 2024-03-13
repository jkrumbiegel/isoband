#include <iostream>
using namespace std;

#include "polygon.h"

template <typename P> 
ostream & operator<<(ostream &out, const P &p) {
  out << "(" << p.x << ", " << p.y << ")";
  return out;
}

template <typename P> 
bool operator==(const P &p1, const P &p2) {
  return (p1.x == p2.x) && (p1.y == p2.y);
}

ostream & operator<<(ostream &out, const in_polygon_type &t) {
  switch(t) {
  case inside:
    out << "inside";
    break;
  case outside:
    out << "outside";
    break;
  default:
    out << "undetermined";
  }
  return out;
}
