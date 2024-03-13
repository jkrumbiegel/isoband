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
#include "isoband.h"

template <typename T> 
void isobander<T>::reset_grid() {
  polygon_grid.clear();

  for (int i=0; i<8; i++) {
    tmp_point_connect[i] = point_connect();
  }
}

template <typename T> 
T isobander<T>::central_value(int r, int c) {// calculates the central value of a given cell
  return (grid_z_p[r + c * nrow] + grid_z_p[r + (c + 1) * nrow] + grid_z_p[r + 1 + c * nrow] + grid_z_p[r + 1 + (c + 1) * nrow])/4;
}

template <typename T> 
void isobander<T>::poly_start(int r, int c, point_type type) { // start a new elementary polygon
  tmp_poly[0].r = r;
  tmp_poly[0].c = c;
  tmp_poly[0].type = type;

  tmp_poly_size = 1;
}

template <typename T> 
void isobander<T>::poly_add(int r, int c, point_type type) { // add point to elementary polygon
  tmp_poly[tmp_poly_size].r = r;
  tmp_poly[tmp_poly_size].c = c;
  tmp_poly[tmp_poly_size].type = type;

  tmp_poly_size++;
}

template <typename T> 
void isobander<T>::poly_merge() { // merge current elementary polygon to prior polygons
  bool to_delete[] = {false, false, false, false, false, false, false, false};

  // first, we figure out the right connections for current polygon
  for (int i = 0; i < tmp_poly_size; i++) {
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
  for (int i = 0; i < tmp_poly_size; i++) {
    const grid_point &p = tmp_poly[i];

    if (to_delete[i]) { // delete point if needed
      polygon_grid.erase(p);
    } else {            // otherwise, copy
      polygon_grid[p] = tmp_point_connect[i];
    }
  }
}

template <typename T>
void isobander<T>::print_polygons_state() {
  for (auto it = polygon_grid.begin(); it != polygon_grid.end(); it++) {
    cout << it->first << ": " << it->second << endl;
  }
  cout << endl;
}

// linear interpolation of boundary intersections
template <typename T>
T isobander<T>::interpolate(T x0, T x1, T z0, T z1, T value) {
  T d = (value - z0) / (z1 - z0);
  T x = x0 + d * (x1 - x0);
  return x;
}

// calculate output coordinates for a given grid point
template <typename T>
point<T> isobander<T>::calc_point_coords(const grid_point &p) {
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
    return point<T>(T(0), T(0)); // should never get here
    }
  }

template <typename T> 
void isobander<T>::set_value(T value_low, T value_high) {
  vlo = value_low;
  vhi = value_high;
}

template <typename T> 
void isobander<T>::calculate_contour() {
  // clear polygon grid and associated internal variables
  reset_grid();

  // setup matrix of ternarized cell representations
  vector<int> ternarized(nrow*ncol);
  vector<int>::iterator iv = ternarized.begin();

  for (int i = 0; i < nrow * ncol; ++i) {
    *iv = (grid_z_p[i] >= vlo && grid_z_p[i] < vhi) + 2*(grid_z_p[i] >= vhi);
    iv++;
  }

  vector<int> cells((nrow - 1) * (ncol - 1));

  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      int index;
      if (!isfinite(grid_z_p[r + c * nrow]) || !isfinite(grid_z_p[r + (c + 1) * nrow]) ||
          !isfinite(grid_z_p[r + 1 + c * nrow]) || !isfinite(grid_z_p[r + 1 + (c + 1) * nrow])) {
        // we don't draw any contours if at least one of the corners is NA
        index = 0;
      } else {
        index = 27*ternarized[r + c * nrow] + 9*ternarized[r + (c + 1) * nrow] + 3*ternarized[r + 1 + (c + 1) * nrow] + ternarized[r + 1 + c * nrow];
      }
      cells[r + c * (nrow - 1)] = index;
    }
  }

  // all polygons must be drawn clockwise for proper merging
  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
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

template <typename T> 
resultStruct<T> isobander<T>::collect() {
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

template <typename T> 
void isoliner<T>::line_merge() { // merge current elementary polygon to prior polygons
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

template <typename T> 
void isoliner<T>::set_value(T value) {
  vlo = value;
}

template <typename T> 
void isoliner<T>::calculate_contour() {
  // clear polygon grid and associated internal variables
  reset_grid();

  // setup matrix of binarized cell representations
  vector<int> binarized(nrow*ncol);
  vector<int>::iterator iv = binarized.begin();
  for (int i = 0; i < nrow * ncol; ++i) {
    *iv = (grid_z_p[i] >= vlo);
    iv++;
  }

  vector<int> cells((nrow - 1) * (ncol - 1));

  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
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

  for (int r = 0; r < nrow-1; r++) {
    for (int c = 0; c < ncol-1; c++) {
      switch(cells[r + c * (nrow - 1)]) {
      case 0: break;
      case 1:
        poly_start(r, c, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      case 2:
        poly_start(r, c+1, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      case 3:
        poly_start(r, c, vintersect_lo);
        poly_add(r, c+1, vintersect_lo);
        line_merge();
        break;
      case 4:
        poly_start(r, c, hintersect_lo);
        poly_add(r, c+1, vintersect_lo);
        line_merge();
        break;
      case 5:
        // like case 2
        poly_start(r, c+1, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        // like case 7
        poly_start(r, c, hintersect_lo);
        poly_add(r, c, vintersect_lo);
        line_merge();
        break;
      case 6:
        poly_start(r, c, hintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      case 7:
        poly_start(r, c, hintersect_lo);
        poly_add(r, c, vintersect_lo);
        line_merge();
        break;
      case 8:
        poly_start(r, c, hintersect_lo);
        poly_add(r, c, vintersect_lo);
        line_merge();
        break;
      case 9:
        poly_start(r, c, hintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      case 10:
        // like case 1
        poly_start(r, c, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        // like case 4
        poly_start(r, c, hintersect_lo);
        poly_add(r, c+1, vintersect_lo);
        line_merge();
        break;
      case 11:
        poly_start(r, c, hintersect_lo);
        poly_add(r, c+1, vintersect_lo);
        line_merge();
        break;
      case 12:
        poly_start(r, c, vintersect_lo);
        poly_add(r, c+1, vintersect_lo);
        line_merge();
        break;
      case 13:
        poly_start(r, c+1, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      case 14:
        poly_start(r, c, vintersect_lo);
        poly_add(r+1, c, hintersect_lo);
        line_merge();
        break;
      default: break; // catch everything, just in case
      }
    }
  }
}

template <typename T> 
resultStruct<T> isoliner<T>::collect() {
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

extern "C" resultStruct<float>* isobands32_impl(float *x, int lenx, float *y, int leny, float *z, int nrow, int ncol, float *values_low, float *values_high, int n_bands) {

  isobander<float> ib(x, lenx, y, leny, z, nrow, ncol, 0.0, 0.0);

  resultStruct<float>* returnstructs = new resultStruct<float>[n_bands];

  for (int i = 0; i < n_bands; ++i) {
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

  for (int i = 0; i < n_values; ++i) {
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

  for (int i = 0; i < n_bands; ++i) {
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

  for (int i = 0; i < n_values; ++i) {
    il.set_value(values[i]);
    il.calculate_contour();

    resultStruct<double> result = il.collect();

    returnstructs[i] = result;
  }

  return returnstructs;
}
