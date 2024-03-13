// This file implements the 2D isoline and isoband algorithms described
// here: https://en.wikipedia.org/wiki/Marching_squares
// Includes merging of line segments and polygons.
// Written by Claus O. Wilke

#define R_NO_REMAP

// #include <R.h>
// #include <Rinternals.h>

#include <iostream>
#include <vector>
#include <unordered_map>
#include <math.h>       /* isfinite */

using namespace std;

#include "polygon.h" // for point

// point in abstract grid space
enum point_type {
  grid,  // point on the original data grid
  hintersect_lo, // intersection with horizontal edge, low value
  hintersect_hi, // intersection with horizontal edge, high value
  vintersect_lo, // intersection with vertical edge, low value
  vintersect_hi  // intersection with vertical edge, high value
};

// return type for extern C functions
template <typename T>  
struct resultStruct {
  T *x;
  T *y;
  int *id;
  int len;
};

struct grid_point {
  int r, c; // row and column
  point_type type; // point type

  // default constructor; negative values indicate non-existing point off grid
  grid_point(int r_in = -1, int c_in = -1, point_type type_in = grid) : r(r_in), c(c_in), type(type_in) {}
  // copy constructor
  grid_point(const grid_point &p) : r(p.r), c(p.c), type(p.type) {}
};

// hash function for grid_point
struct grid_point_hasher {
  size_t operator()(const grid_point& p) const
  {
    // this should work up to about 100,000,000 rows/columns
    return hash<long long>()(
      (static_cast<long long>(p.r) << 30) ^
        (static_cast<long long>(p.c) << 3) ^
          static_cast<long long>(p.type));
  }
};

bool operator==(const grid_point &p1, const grid_point &p2) {
  return (p1.r == p2.r) && (p1.c == p2.c) && (p1.type == p2.type);
}

ostream & operator<<(ostream &out, const grid_point &p) {
  out << "(" << p.c << ", " << p.r << ", " << p.type << ")";
  return out;
}

// connection between points in grid space
struct point_connect {
  grid_point prev, next; // previous and next points in polygon
  grid_point prev2, next2; // alternative previous and next, when two separate polygons have vertices on the same grid point

  bool altpoint;  // does this connection hold an alternative point?
  bool collected, collected2; // has this connection been collected into a final polygon?

  point_connect() : altpoint(false), collected(false), collected2(false) {};
};

ostream & operator<<(ostream &out, const point_connect &pc) {
  out << "prev: " << pc.prev << "; next: " << pc.next << " ";
  if (pc.altpoint) {
    out << "AP prev: " << pc.prev2 << "; next2: " << pc.next2 << " ";
  }
  return out;
}

template <typename T>  
class isobander {
protected:
  int nrow, ncol; // numbers of rows and columns
  // SEXP grid_x, grid_y, grid_z;
  T *grid_x_p, *grid_y_p, *grid_z_p;
  T vlo, vhi; // low and high cutoff values
  grid_point tmp_poly[8]; // temp storage for elementary polygons; none has more than 8 vertices
  point_connect tmp_point_connect[8];
  int tmp_poly_size; // current number of elements in tmp_poly

  typedef unordered_map<grid_point, point_connect, grid_point_hasher> gridmap;
  gridmap polygon_grid;

  bool interrupted;

  void reset_grid() {
    polygon_grid.clear();

    for (auto i=0; i<8; i++) {
      tmp_point_connect[i] = point_connect();
    }
  }

  // internal member functions

  T central_value(int r, int c) {// calculates the central value of a given cell
    return (grid_z_p[r + c * nrow] + grid_z_p[r + (c + 1) * nrow] + grid_z_p[r + 1 + c * nrow] + grid_z_p[r + 1 + (c + 1) * nrow])/4;
  }

  void poly_start(int r, int c, point_type type) { // start a new elementary polygon
    tmp_poly[0].r = r;
    tmp_poly[0].c = c;
    tmp_poly[0].type = type;

    tmp_poly_size = 1;
  }

  void poly_add(int r, int c, point_type type) { // add point to elementary polygon
    tmp_poly[tmp_poly_size].r = r;
    tmp_poly[tmp_poly_size].c = c;
    tmp_poly[tmp_poly_size].type = type;

    tmp_poly_size++;
  }

  void poly_merge() { // merge current elementary polygon to prior polygons
    bool to_delete[] = {false, false, false, false, false, false, false, false};

    // first, we figure out the right connections for current polygon
    for (auto i = 0; i < tmp_poly_size; i++) {
      // create defined state in tmp_point_connect[]
      // for each point, find previous and next point in polygon
      tmp_point_connect[i].altpoint = false;
      tmp_point_connect[i].next = tmp_poly[(i+1<tmp_poly_size) ? i+1 : 0];
      tmp_point_connect[i].prev = tmp_poly[(i-1>=0) ? i-1 : tmp_poly_size-1];

      // now merge with existing polygons if needed
      const grid_point &p = tmp_poly[i];
      if (polygon_grid.count(p) > 0) { // point has been used before, need to merge polygons
        if (!polygon_grid[p].altpoint) {
          // basic scenario, no alternative point at this location
          int score = 2 * (tmp_point_connect[i].next == polygon_grid[p].prev) + (tmp_point_connect[i].prev == polygon_grid[p].next);
          switch (score) {
          case 3: // 11
            // both prev and next cancel, point can be deleted
            to_delete[i] = true;
            break;
          case 2: // 10
            // merge in "next" direction
            tmp_point_connect[i].next = polygon_grid[p].next;
            break;
          case 1: // 01
            // merge in "prev" direction
            tmp_point_connect[i].prev = polygon_grid[p].prev;
            break;
          default: // 00
            // if we get here, we have two polygon vertices sharing the same grid location
            // in an unmergable configuration; need to store both
            tmp_point_connect[i].prev2 = polygon_grid[p].prev;
            tmp_point_connect[i].next2 = polygon_grid[p].next;
            tmp_point_connect[i].altpoint = true;
          }
        } else {
          // case with alternative point at this location
          int score =
            8 * (tmp_point_connect[i].next == polygon_grid[p].prev2) + 4 * (tmp_point_connect[i].prev == polygon_grid[p].next2) +
            2 * (tmp_point_connect[i].next == polygon_grid[p].prev) + (tmp_point_connect[i].prev == polygon_grid[p].next);
          switch (score) {
          case 9: // 1001
            // three-way merge
            tmp_point_connect[i].next = polygon_grid[p].next2;
            tmp_point_connect[i].prev = polygon_grid[p].prev;
            break;
          case 6: // 0110
            // three-way merge
            tmp_point_connect[i].next = polygon_grid[p].next;
            tmp_point_connect[i].prev = polygon_grid[p].prev2;
            break;
          case 8: // 1000
            // two-way merge with alt point only
            // set up merged alt point
            tmp_point_connect[i].next2 = polygon_grid[p].next2;
            tmp_point_connect[i].prev2 = tmp_point_connect[i].prev;
            // copy over existing point as is
            tmp_point_connect[i].prev = polygon_grid[p].prev;
            tmp_point_connect[i].next = polygon_grid[p].next;
            tmp_point_connect[i].altpoint = true;
            break;
          case 4: // 0100
            // two-way merge with alt point only
            // set up merged alt point
            tmp_point_connect[i].prev2 = polygon_grid[p].prev2;
            tmp_point_connect[i].next2 = tmp_point_connect[i].next;
            // copy over existing point as is
            tmp_point_connect[i].prev = polygon_grid[p].prev;
            tmp_point_connect[i].next = polygon_grid[p].next;
            tmp_point_connect[i].altpoint = true;
            break;
          case 2: // 0010
            // two-way merge with original point only
            // merge point
            tmp_point_connect[i].next = polygon_grid[p].next;
            // copy over existing alt point as is
            tmp_point_connect[i].prev2 = polygon_grid[p].prev2;
            tmp_point_connect[i].next2 = polygon_grid[p].next2;
            tmp_point_connect[i].altpoint = true;
            break;
          case 1: // 0100
            // two-way merge with original point only
            // merge point
            tmp_point_connect[i].prev = polygon_grid[p].prev;
            // copy over existing alt point as is
            tmp_point_connect[i].prev2 = polygon_grid[p].prev2;
            tmp_point_connect[i].next2 = polygon_grid[p].next2;
            tmp_point_connect[i].altpoint = true;
            break;
          default:
            throw std::runtime_error("undefined merging configuration");
          }
        }
      }
    }

    // then we copy the connections into the polygon matrix
    for (auto i = 0; i < tmp_poly_size; i++) {
      const grid_point &p = tmp_poly[i];

      if (to_delete[i]) { // delete point if needed
        polygon_grid.erase(p);
      } else {            // otherwise, copy
        polygon_grid[p] = tmp_point_connect[i];
      }
    }

    // print_polygons_state();
  }

  void print_polygons_state() {
    for (auto it = polygon_grid.begin(); it != polygon_grid.end(); it++) {
      cout << it->first << ": " << it->second << endl;
    }
    cout << endl;
  }

  // linear interpolation of boundary intersections
  T interpolate(T x0, T x1, T z0, T z1, T value) {
    T d = (value - z0) / (z1 - z0);
    T x = x0 + d * (x1 - x0);
    return x;
  }

  // calculate output coordinates for a given grid point
  point<T> calc_point_coords(const grid_point &p) {
    switch(p.type) {
    case grid:
      return point<T>(grid_x_p[p.c], grid_y_p[p.r]);
    case hintersect_lo: // intersection with horizontal edge, low value
      return point<T>(interpolate(grid_x_p[p.c], grid_x_p[p.c+1], grid_z_p[p.r + p.c * nrow], grid_z_p[p.r + (p.c + 1) * nrow], vlo), grid_y_p[p.r]);
    case hintersect_hi: // intersection with horizontal edge, high value
      return point<T>(interpolate(grid_x_p[p.c], grid_x_p[p.c+1], grid_z_p[p.r + p.c * nrow], grid_z_p[p.r + (p.c + 1) * nrow], vhi), grid_y_p[p.r]);
    case vintersect_lo: // intersection with vertical edge, low value
      return point<T>(grid_x_p[p.c], interpolate(grid_y_p[p.r], grid_y_p[p.r+1], grid_z_p[p.r + p.c * nrow], grid_z_p[p.r + 1 + p.c * nrow], vlo));
    case vintersect_hi: // intersection with vertical edge, high value
      return point<T>(grid_x_p[p.c], interpolate(grid_y_p[p.r], grid_y_p[p.r+1], grid_z_p[p.r + p.c * nrow], grid_z_p[p.r + 1 + p.c * nrow], vhi));
    default:
      return point<T>(0, 0); // should never get here
    }
  }

public:
  isobander(T *x, int lenx, T *y, int leny, T *z, int nrow, int ncol, T value_low = 0, T value_high = 0) :
    grid_x_p(x), grid_y_p(y), grid_z_p(z), nrow(nrow), ncol(ncol),
    vlo(value_low), vhi(value_high), interrupted(false)
  {

    if (lenx != ncol) {throw std::invalid_argument("Number of x coordinates must match number of columns in density matrix.");}
    if (leny != nrow) {throw std::invalid_argument("Number of y coordinates must match number of rows in density matrix.");}
  }

  virtual ~isobander() {}

  bool was_interrupted() {return interrupted;}

  void set_value(T value_low, T value_high) {
    vlo = value_low;
    vhi = value_high;
  }

  virtual void calculate_contour() {
    // clear polygon grid and associated internal variables
    reset_grid();

    // setup matrix of ternarized cell representations
    vector<int> ternarized(nrow*ncol);
    vector<int>::iterator iv = ternarized.begin();

    for (auto i = 0; i < nrow * ncol; ++i) {
      *iv = (grid_z_p[i] >= vlo && grid_z_p[i] < vhi) + 2*(grid_z_p[i] >= vhi);
      iv++;
    }


    vector<int> cells((nrow - 1) * (ncol - 1));

    for (auto r = 0; r < nrow-1; r++) {
      for (auto c = 0; c < ncol-1; c++) {
        int index;
        if (!isfinite(grid_z_p[r + c * nrow]) || !isfinite(grid_z_p[r + (c + 1) * nrow]) ||
            !isfinite(grid_z_p[r + 1 + c * nrow]) || !isfinite(grid_z_p[r + 1 + (c + 1) * nrow])) {
          // we don't draw any contours if at least one of the corners is NA
          index = 0;
        } else {
          index = 27*ternarized[r + c * nrow] + 9*ternarized[r + (c + 1) * nrow] + 3*ternarized[r + 1 + (c + 1) * nrow] + ternarized[r + 1 + c * nrow];
        }
        cells[r + c * (nrow - 1)] = index;
        //cout << index << " ";
      }
      //cout << endl;
    }
    // if (checkInterrupt()) {
    //   interrupted = true;
    //   return;
    // }

    // all polygons must be drawn clockwise for proper merging
    for (auto r = 0; r < nrow-1; r++) {
      for (auto c = 0; c < ncol-1; c++) {
        //cout << r << " " << c << " " << cells(r, c) << endl;
        switch(cells[r + c * (nrow - 1)]) {
        // doing cases out of order, sorted by type, is easier to keep track of

        // no contour
        case 0: break;
        case 80: break;

        // single triangle
        case 1: // 0001
          poly_start(r, c, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 3: // 0010
          poly_start(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 9: // 0100
          poly_start(r, c, hintersect_lo);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_merge();
          break;
        case 27: // 1000
          poly_start(r, c, vintersect_lo);
          poly_add(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 79: // 2221
          poly_start(r, c, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 77: // 2212
          poly_start(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 71: // 2122
          poly_start(r, c, hintersect_hi);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_merge();
          break;
        case 53: // 1222
          poly_start(r, c, vintersect_hi);
          poly_add(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_merge();
          break;

          // single trapezoid
        case 78: // 2220
          poly_start(r, c, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 74: // 2202
          poly_start(r+1, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 62: // 2022
          poly_start(r, c+1, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_merge();
          break;
        case 26: // 0222
          poly_start(r, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 2: // 0002
          poly_start(r, c, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 6: // 0020
          poly_start(r+1, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 18: // 0200
          poly_start(r, c+1, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_merge();
          break;
        case 54: // 2000
          poly_start(r, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_merge();
          break;

          // single rectangle
        case 4: // 0011
          poly_start(r, c, vintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 12: // 0110
          poly_start(r, c, hintersect_lo);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 36: // 1100
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 28: // 1001
          poly_start(r, c, hintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, grid);
          poly_add(r, c, grid);
          poly_merge();
          break;
        case 76: // 2211
          poly_start(r, c, vintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 68: // 2112
          poly_start(r, c, hintersect_hi);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 44: // 1122
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 52: // 1221
          poly_start(r, c, hintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, grid);
          poly_add(r, c, grid);
          poly_merge();
          break;
        case 72: // 2200
          poly_start(r, c, vintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 56: // 2002
          poly_start(r, c, hintersect_hi);
          poly_add(r, c, hintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 8: // 0022
          poly_start(r, c, vintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 24: // 0220
          poly_start(r, c, hintersect_lo);
          poly_add(r, c, hintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;

        // single square
        case 40: // 1111
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;

        // single pentagon
        case 49: // 1211
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 67: // 2111
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_merge();
          break;
        case 41: // 1112
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 43: // 1121
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 31: // 1011
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 13: // 0111
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_merge();
          break;
        case 39: // 1110
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 37: // 1101
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 45: // 1200
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 15: // 0120
          poly_start(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 5: // 0012
          poly_start(r, c, vintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 55: // 2001
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c, hintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 35: // 1022
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 65: // 2102
          poly_start(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_merge();
          break;
        case 75: // 2210
          poly_start(r, c, vintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 25: // 0221
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c, hintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 29: // 1002
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 63: // 2100
          poly_start(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_merge();
          break;
        case 21: // 0210
          poly_start(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_merge();
          break;
        case 7: // 0021
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_merge();
          break;
        case 51: // 1220
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 17: // 0122
          poly_start(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 59: // 2012
          poly_start(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_merge();
          break;
        case 73: // 2201
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;

          // single hexagon
        case 22: // 0211
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_merge();
          break;
        case 66: // 2110
          poly_start(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_merge();
          break;
        case 38: // 1102
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 34: // 1021
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 58: // 2011
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_merge();
          break;
        case 14: // 0112
          poly_start(r, c+1, grid);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 42: // 1120
          poly_start(r, c, grid);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;
        case 46: // 1201
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r+1, c, grid);
          poly_merge();
          break;
        case 64: // 2101
          poly_start(r+1, c, grid);
          poly_add(r, c, vintersect_hi);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, grid);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c, hintersect_lo);
          poly_merge();
          break;
        case 16: // 0121
          poly_start(r, c+1, grid);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r+1, c, grid);
          poly_add(r, c, vintersect_lo);
          poly_add(r, c, hintersect_lo);
          poly_merge();
          break;
        case 32: // 1012
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_lo);
          poly_add(r, c+1, vintersect_lo);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_hi);
          poly_add(r, c, vintersect_hi);
          poly_merge();
          break;
        case 48: // 1210
          poly_start(r, c, grid);
          poly_add(r, c, hintersect_hi);
          poly_add(r, c+1, vintersect_hi);
          poly_add(r+1, c+1, grid);
          poly_add(r+1, c, hintersect_lo);
          poly_add(r, c, vintersect_lo);
          poly_merge();
          break;

        // 6-sided saddle
        case 10: // 0101
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_merge();
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_merge();
            } else {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c+1, grid);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_merge();
            }
          }
          break;
        case 30: // 1010
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_merge();
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_merge();
            } else {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_merge();
            }
          }
          break;
        case 70: // 2121
          {
            T vc = central_value(r, c);
            if (vc >= vhi) {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_merge();
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_merge();
            } else {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, grid);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_merge();
            }
          }
          break;
        case 50: // 1212
          {
            T vc = central_value(r, c);
            if (vc >= vhi) {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
            } else {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
            }
          }
          break;

        // 7-sided saddle
        case 69: // 2120
          {
            T vc = central_value(r, c);
            if (vc >= vhi) {
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_merge();
              poly_start(r, c, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_merge();
            } else {
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_merge();
            }
          }
          break;
        case 61: // 2021
          {
            T vc = central_value(r, c);
              if (vc >= vhi) {
                poly_start(r+1, c, grid);
                poly_add(r, c, vintersect_hi);
                poly_add(r+1, c, hintersect_hi);
                poly_merge();
                poly_start(r, c+1, vintersect_hi);
                poly_add(r, c, hintersect_hi);
                poly_add(r, c, hintersect_lo);
                poly_add(r, c+1, vintersect_lo);
                poly_merge();
              } else {
                poly_start(r+1, c, grid);
                poly_add(r, c, vintersect_hi);
                poly_add(r, c, hintersect_hi);
                poly_add(r, c, hintersect_lo);
                poly_add(r, c+1, vintersect_lo);
                poly_add(r, c+1, vintersect_hi);
                poly_add(r+1, c, hintersect_hi);
                poly_merge();
              }
            }
          break;
        case 47: // 1202
          {
            T vc = central_value(r, c);
            if (vc >= vhi) {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
              poly_start(r+1, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_merge();
            } else {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
            }
          }
          break;
        case 23: // 0212
          {
            T vc = central_value(r, c);
            if (vc >= vhi) {
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
              poly_start(r, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_merge();
            } else {
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
            }
          }
          break;
        case 11: // 0102
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_merge();
              poly_start(r, c, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
            } else {
              poly_start(r, c+1, grid);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_merge();
            }
          }
          break;
        case 19: // 0201
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_merge();
              poly_start(r, c+1, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
            } else {
              poly_start(r+1, c, grid);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_merge();
            }
          }
          break;
        case 33: // 1020
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_merge();
              poly_start(r+1, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_merge();
            } else {
              poly_start(r, c, grid);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_merge();
            }
          }
          break;
        case 57: // 2010
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_merge();
              poly_start(r, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_merge();
            } else {
              poly_start(r+1, c+1, grid);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r, c, vintersect_lo);
              poly_add(r, c, vintersect_hi);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c+1, vintersect_lo);
              poly_merge();
            }
          }
          break;

        // 8-sided saddle
      case 60: // 2020
        {
          T vc = central_value(r, c);
          if (vc < vlo) {
            poly_start(r, c, vintersect_hi);
            poly_add(r, c, hintersect_hi);
            poly_add(r, c, hintersect_lo);
            poly_add(r, c, vintersect_lo);
            poly_merge();
            poly_start(r, c+1, vintersect_hi);
            poly_add(r+1, c, hintersect_hi);
            poly_add(r+1, c, hintersect_lo);
            poly_add(r, c+1, vintersect_lo);
            poly_merge();
          } else if (vc >= vhi) {
            poly_start(r, c, vintersect_hi);
            poly_add(r+1, c, hintersect_hi);
            poly_add(r+1, c, hintersect_lo);
            poly_add(r, c, vintersect_lo);
            poly_merge();
            poly_start(r, c+1, vintersect_hi);
            poly_add(r, c, hintersect_hi);
            poly_add(r, c, hintersect_lo);
            poly_add(r, c+1, vintersect_lo);
            poly_merge();
          } else {
            poly_start(r, c, vintersect_hi);
            poly_add(r, c, hintersect_hi);
            poly_add(r, c, hintersect_lo);
            poly_add(r, c+1, vintersect_lo);
            poly_add(r, c+1, vintersect_hi);
            poly_add(r+1, c, hintersect_hi);
            poly_add(r+1, c, hintersect_lo);
            poly_add(r, c, vintersect_lo);
            poly_merge();
          }
        }
        break;
        case 20: // 0202
          {
            T vc = central_value(r, c);
            if (vc < vlo) {
              poly_start(r, c, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
              poly_start(r, c+1, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
            } else if (vc >= vhi) {
              poly_start(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
              poly_start(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_merge();
            } else {
              poly_start(r, c, vintersect_lo);
              poly_add(r, c, hintersect_lo);
              poly_add(r, c, hintersect_hi);
              poly_add(r, c+1, vintersect_hi);
              poly_add(r, c+1, vintersect_lo);
              poly_add(r+1, c, hintersect_lo);
              poly_add(r+1, c, hintersect_hi);
              poly_add(r, c, vintersect_hi);
              poly_merge();
            }
          }
          break;
        }
      }
    }
  }

  virtual resultStruct<T> collect() {
    // Early exit if calculate_contour was interrupted
    // if (was_interrupted()) {
    //   return R_NilValue;
    // }

    // make polygons
    vector<T> x_out, y_out; vector<int> id;  // vectors holding resulting polygon paths
    int cur_id = 0;           // id counter for the polygon lines

    // iterate over all locations in the polygon grid
    for (auto it = polygon_grid.begin(); it != polygon_grid.end(); it++) {
      if (((it->second).collected && !(it->second).altpoint) ||
          ((it->second).collected && (it->second).collected2 && (it->second).altpoint)) {
        continue; // skip any grid points that are already fully collected
      }

      // we have found a new polygon line; process it
      cur_id++;

      grid_point start = it->first;
      grid_point cur = start;
      grid_point prev = (it->second).prev;
      // if this point has an alternative and it hasn't been collected yet then we start there
      if ((it->second).altpoint && !(it->second).collected2) prev = (it->second).prev2;

      int i = 0;
      do {
        point<T> p = calc_point_coords(cur);
        x_out.push_back(p.x);
        y_out.push_back(p.y);
        id.push_back(cur_id);

        // record that we have processed this point and proceed to next
        if (polygon_grid[cur].altpoint && polygon_grid[cur].prev2 == prev) {
          // if an alternative point exists and its previous point in the polygon
          // corresponds to the recorded previous point, then that's the point
          // we're working with here

          // mark current point as collected and advance
          polygon_grid[cur].collected2 = true;
          grid_point newcur = polygon_grid[cur].next2;
          prev = cur;
          cur = newcur;
        } else {
          // mark current point as collected and advance
          polygon_grid[cur].collected = true;
          grid_point newcur = polygon_grid[cur].next;
          prev = cur;
          cur = newcur;
        }
        i++;
      } while (!(cur == start)); // keep going until we reach the start point again
    }

    int len = x_out.size();

    T* xs = new T[len];
    T* ys = new T[len];
    int* ids = new int[len];

    copy(x_out.begin(), x_out.end(), xs);
    copy(y_out.begin(), y_out.end(), ys);
    copy(id.begin(), id.end(), ids);

    return resultStruct<T>{xs, ys, ids, len};
  }
};


template <typename T>  
class isoliner : public isobander<T> {
  using isobander<T>::tmp_poly_size;
  using isobander<T>::polygon_grid;
  using isobander<T>::tmp_poly;
  using isobander<T>::vlo;
  using isobander<T>::vhi;
  using isobander<T>::nrow;
  using isobander<T>::ncol;
  using isobander<T>::grid_x_p;
  using isobander<T>::grid_y_p;
  using isobander<T>::grid_z_p;
  using isobander<T>::reset_grid;
  using isobander<T>::central_value;
  using isobander<T>::calc_point_coords;
  using isobander<T>::poly_start;
  using isobander<T>::poly_add;

protected:
  void line_start(int r, int c, point_type type) { // start a new line segment
    tmp_poly[0].r = r;
    tmp_poly[0].c = c;
    tmp_poly[0].type = type;

    tmp_poly_size = 1;
  }

  void line_add(int r, int c, point_type type) { // add point to line
    tmp_poly[tmp_poly_size].r = r;
    tmp_poly[tmp_poly_size].c = c;
    tmp_poly[tmp_poly_size].type = type;

    tmp_poly_size++;
  }

  void line_merge() { // merge current elementary polygon to prior polygons
    int score = 2*polygon_grid.count(tmp_poly[1]) + polygon_grid.count(tmp_poly[0]);

    switch(score) {
    case 0: // completely unconnected line segment
      polygon_grid[tmp_poly[0]].next = tmp_poly[1];
      polygon_grid[tmp_poly[1]].prev = tmp_poly[0];
      break;
    case 1: // only first point connects
      if (polygon_grid[tmp_poly[0]].next == grid_point()) {
        polygon_grid[tmp_poly[0]].next = tmp_poly[1];
        polygon_grid[tmp_poly[1]].prev = tmp_poly[0];
      } else if (polygon_grid[tmp_poly[0]].prev == grid_point()) {
        polygon_grid[tmp_poly[0]].prev = tmp_poly[1];
        polygon_grid[tmp_poly[1]].next = tmp_poly[0];
      } else {
        // should never go here
        throw std::runtime_error("cannot merge line segment at interior of existing line segment");
      }
      break;
    case 2: // only second point connects
      if (polygon_grid[tmp_poly[1]].next == grid_point()) {
        polygon_grid[tmp_poly[1]].next = tmp_poly[0];
        polygon_grid[tmp_poly[0]].prev = tmp_poly[1];
      } else if (polygon_grid[tmp_poly[1]].prev == grid_point()) {
        polygon_grid[tmp_poly[1]].prev = tmp_poly[0];
        polygon_grid[tmp_poly[0]].next = tmp_poly[1];
      } else {
        // should never go here
        throw std::runtime_error("cannot merge line segment at interior of existing line segment");
      }
      break;
    case 3: // two-way merge
      //break; // two-way merge doesn't work yet
      {
        int score2 =
          8*(polygon_grid[tmp_poly[0]].next == grid_point()) +
          4*(polygon_grid[tmp_poly[0]].prev == grid_point()) +
          2*(polygon_grid[tmp_poly[1]].next == grid_point()) +
          (polygon_grid[tmp_poly[1]].prev == grid_point());

        switch(score2) {
        case 9: // 1001
          polygon_grid[tmp_poly[0]].next = tmp_poly[1];
          polygon_grid[tmp_poly[1]].prev = tmp_poly[0];
          break;
        case 6: // 0110
          polygon_grid[tmp_poly[0]].prev = tmp_poly[1];
          polygon_grid[tmp_poly[1]].next = tmp_poly[0];
          break;
        case 10: // 1010
          {
            polygon_grid[tmp_poly[0]].next = tmp_poly[1];
            polygon_grid[tmp_poly[1]].next = tmp_poly[0];

            // need to reverse connections
            grid_point cur = tmp_poly[1];
            int i = 0;
            do {
              grid_point tmp = polygon_grid[cur].prev;
              polygon_grid[cur].prev = polygon_grid[cur].next;
              polygon_grid[cur].next = tmp;
              cur = tmp;
              i++;
            } while (!(cur == grid_point()));
          }
          break;
        case 5: // 0101
          {
            polygon_grid[tmp_poly[0]].prev = tmp_poly[1];
            polygon_grid[tmp_poly[1]].prev = tmp_poly[0];

            // need to reverse connections
            grid_point cur = tmp_poly[0];
            int i = 0;
            do {
              grid_point tmp = polygon_grid[cur].next;
              polygon_grid[cur].next = polygon_grid[cur].prev;
              polygon_grid[cur].prev = tmp;
              cur = tmp;
              i++;
            } while (!(cur == grid_point()));
          }
          break;
        default:  // should never go here
          throw std::runtime_error("cannot merge line segment at interior of existing line segment");
        }
      }
    break;
    default:
      throw std::runtime_error("unknown merge state");
    }
  }

public:
  isoliner(T *x, int lenx, T *y, int leny, T *z, int nrow, int ncol, T value = 0) :
    isobander<T>(x, lenx, y, leny, z, nrow, ncol, value, 0) {}

  void set_value(T value) {
    vlo = value;
  }

  virtual void calculate_contour() {
    // clear polygon grid and associated internal variables
    reset_grid();

    // setup matrix of binarized cell representations
    vector<int> binarized(nrow*ncol);
    vector<int>::iterator iv = binarized.begin();
    for (auto i = 0; i < nrow * ncol; ++i) {
      *iv = (grid_z_p[i] >= vlo);
      iv++;
    }

    vector<int> cells((nrow - 1) * (ncol - 1));

    for (auto r = 0; r < nrow-1; r++) {
      for (auto c = 0; c < ncol-1; c++) {
        int index;
        if (!isfinite(grid_z_p[r + c * nrow]) || !isfinite(grid_z_p[r + (c + 1) * nrow]) ||
            !isfinite(grid_z_p[r + 1 + c * nrow]) || !isfinite(grid_z_p[r + 1 + (c + 1) * nrow])) {
          // we don't draw any contours if at least one of the corners is NA
          index = 0;
        } else {
          index = 8*binarized[r + c * nrow] + 4*binarized[r + (c + 1) * nrow] + 2*binarized[r + 1 + (c + 1) * nrow] + 1*binarized[r + 1 + c * nrow];
        }

        // two-segment saddles
        if (index == 5 && (central_value(r, c) < vlo)) {
          index = 10;
        } else if (index == 10 && (central_value(r, c) < vlo)) {
          index = 5;
        }

        cells[r + c * (nrow - 1)] = index;
      }
    }

    // if (checkInterrupt()) {
    //   interrupted = true;
    //   return;
    // }

    for (auto r = 0; r < nrow-1; r++) {
      for (auto c = 0; c < ncol-1; c++) {
        switch(cells[r + c * (nrow - 1)]) {
        case 0: break;
        case 1:
          line_start(r, c, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        case 2:
          line_start(r, c+1, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        case 3:
          line_start(r, c, vintersect_lo);
          line_add(r, c+1, vintersect_lo);
          line_merge();
          break;
        case 4:
          line_start(r, c, hintersect_lo);
          line_add(r, c+1, vintersect_lo);
          line_merge();
          break;
        case 5:
          // like case 2
          line_start(r, c+1, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          // like case 7
          line_start(r, c, hintersect_lo);
          line_add(r, c, vintersect_lo);
          line_merge();
          break;
        case 6:
          line_start(r, c, hintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        case 7:
          line_start(r, c, hintersect_lo);
          line_add(r, c, vintersect_lo);
          line_merge();
          break;
        case 8:
          line_start(r, c, hintersect_lo);
          line_add(r, c, vintersect_lo);
          line_merge();
          break;
        case 9:
          line_start(r, c, hintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        case 10:
          // like case 1
          line_start(r, c, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          // like case 4
          line_start(r, c, hintersect_lo);
          line_add(r, c+1, vintersect_lo);
          line_merge();
          break;
        case 11:
          line_start(r, c, hintersect_lo);
          line_add(r, c+1, vintersect_lo);
          line_merge();
          break;
        case 12:
          line_start(r, c, vintersect_lo);
          line_add(r, c+1, vintersect_lo);
          line_merge();
          break;
        case 13:
          line_start(r, c+1, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        case 14:
          line_start(r, c, vintersect_lo);
          line_add(r+1, c, hintersect_lo);
          line_merge();
          break;
        default: break; // catch everything, just in case
        }
      }
    }
  }

  virtual resultStruct<T> collect() {
    // // Early exit if calculate_contour was interrupted
    // if (was_interrupted()) {
    //   return R_NilValue;
    // }

    // make line segments
    vector<T> x_out, y_out; vector<int> id;  // vectors holding resulting polygon paths
    int cur_id = 0;           // id counter for individual line segments

    // iterate over all locations in the polygon grid
    for (auto it = polygon_grid.begin(); it != polygon_grid.end(); it++) {
      if ((it->second).collected) {
        continue; // skip any grid points that are already collected
      }

      // we have found a new polygon line; process it
      cur_id++;

      grid_point start = it->first;
      grid_point cur = start;

      int i = 0;
      if (!(polygon_grid[cur].prev == grid_point())) {
        // back-track until we find the beginning of the line or circle around once
        do {
          cur = polygon_grid[cur].prev;
          i++;
        } while (!(cur == start || polygon_grid[cur].prev == grid_point()));
      }

      start = cur; // reset starting point
      i = 0;
      do {
        point<T> p = calc_point_coords(cur);

        x_out.push_back(p.x);
        y_out.push_back(p.y);
        id.push_back(cur_id);

        // record that we have processed this point and proceed to next
        polygon_grid[cur].collected = true;
        cur = polygon_grid[cur].next;
        i++;
      } while (!(cur == start || cur == grid_point())); // keep going until we reach the start point again
      // if we're back to start, need to output that point one more time
      if (cur == start) {
        point<T> p = calc_point_coords(cur);
        x_out.push_back(p.x);
        y_out.push_back(p.y);
        id.push_back(cur_id);
      }
    }

    int len = x_out.size();

    T* xs = new T[len];
    T* ys = new T[len];
    int* ids = new int[len];

    copy(x_out.begin(), x_out.end(), xs);
    copy(y_out.begin(), y_out.end(), ys);
    copy(id.begin(), id.end(), ids);

    return resultStruct<T>{xs, ys, ids, len};
  }
};
extern "C" resultStruct<float>* isobands32_impl(float *x, int lenx, float *y, int leny, float *z, int nrow, int ncol, float *values_low, float *values_high, int n_bands) {

  isobander<float> ib(x, lenx, y, leny, z, nrow, ncol, 0.0, 0.0);

  resultStruct<float>* returnstructs = new resultStruct<float>[n_bands];

  for (auto i = 0; i < n_bands; ++i) {
    ib.set_value(values_low[i], values_high[i]);
    ib.calculate_contour();

    resultStruct<float> result = ib.collect();

    returnstructs[i] = result;
  }

  return returnstructs;
}

extern "C" resultStruct<float>* isolines32_impl(float *x, int lenx, float *y, int leny, float *z, int nrow, int ncol, float *values, int n_values) {

  isoliner<float> il(x, lenx, y, leny, z, nrow, ncol);

  resultStruct<float>* returnstructs = new resultStruct<float>[n_values];

  for (auto i = 0; i < n_values; ++i) {
    il.set_value(values[i]);
    il.calculate_contour();

    resultStruct<float> result = il.collect();

    returnstructs[i] = result;
  }

  return returnstructs;
}

extern "C" resultStruct<double>* isobands_impl(double *x, int lenx, double *y, int leny, double *z, int nrow, int ncol, double *values_low, double *values_high, int n_bands) {

  isobander<double> ib(x, lenx, y, leny, z, nrow, ncol, 0.0, 0.0);

  resultStruct<double>* returnstructs = new resultStruct<double>[n_bands];

  for (auto i = 0; i < n_bands; ++i) {
    ib.set_value(values_low[i], values_high[i]);
    ib.calculate_contour();

    resultStruct<double> result = ib.collect();

    returnstructs[i] = result;
  }

  return returnstructs;
}

extern "C" resultStruct<double>* isolines_impl(double *x, int lenx, double *y, int leny, double *z, int nrow, int ncol, double *values, int n_values) {

  isoliner<double> il(x, lenx, y, leny, z, nrow, ncol);

  resultStruct<double>* returnstructs = new resultStruct<double>[n_values];

  for (auto i = 0; i < n_values; ++i) {
    il.set_value(values[i]);
    il.calculate_contour();

    resultStruct<double> result = il.collect();

    returnstructs[i] = result;
  }

  return returnstructs;
}