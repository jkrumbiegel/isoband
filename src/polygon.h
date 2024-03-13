#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <unordered_map>
#include <ostream>

using namespace std;

// point in x-y space
template <typename T> 
struct point {
  T x, y; // x and y coordinates

  point(T x_in = 0, T y_in = 0) : x(x_in), y(y_in) {}
};

using pointd = point<double>;
using pointf = point<float>;

bool operator==(const pointd &p1, const pointd &p2);
ostream & operator<<(ostream &out, const pointd &p);

bool operator==(const pointf &p1, const pointf &p2);
ostream & operator<<(ostream &out, const pointf &p);

typedef vector<pointd> polygon;

enum in_polygon_type {
  inside,       // point is inside a polygon
  outside,      // point is outside a polygon
  undetermined // point lies right on the boundary
};

#endif // POLYGON_H
