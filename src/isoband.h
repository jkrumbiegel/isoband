#ifndef ISOBAND_H
#define ISOBAND_H

// return type for extern C functions
template <typename T> 
struct resultStruct {
  T *x;
  T *y;
  int *id;
  int len;
};

// point in abstract grid space
enum point_type {
  grid,  // point on the original data grid
  hintersect_lo, // intersection with horizontal edge, low value
  hintersect_hi, // intersection with horizontal edge, high value
  vintersect_lo, // intersection with vertical edge, low value
  vintersect_hi  // intersection with vertical edge, high value
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
  T *grid_x_p, *grid_y_p, *grid_z_p;
  T vlo, vhi; // low and high cutoff values
  grid_point tmp_poly[8]; // temp storage for elementary polygons; none has more than 8 vertices
  point_connect tmp_point_connect[8];
  int tmp_poly_size; // current number of elements in tmp_poly

  typedef unordered_map<grid_point, point_connect, grid_point_hasher> gridmap;
  gridmap polygon_grid;

  bool interrupted;

  void reset_grid();
  // internal member functions

  T central_value(int r, int c);

  void poly_start(int r, int c, point_type type);
  void poly_add(int r, int c, point_type type);
  void poly_merge();

  void print_polygons_state();
  // linear interpolation of boundary intersections
  T interpolate(T x0, T x1, T z0, T z1, T value);
  point<T> calc_point_coords(const grid_point &p);

public:
  isobander(T *x, int lenx, T *y, int leny, T *z, int nrow, int ncol, T value_low = 0.0, T value_high = 0.0) :
    grid_x_p(x), grid_y_p(y), grid_z_p(z), nrow(nrow), ncol(ncol),
    vlo(value_low), vhi(value_high), interrupted(false)
  {
    if (lenx != ncol) {throw std::invalid_argument("Number of x coordinates must match number of columns in density matrix.");}
    if (leny != nrow) {throw std::invalid_argument("Number of y coordinates must match number of rows in density matrix.");}
  }

  virtual ~isobander() {};

  bool was_interrupted() {return interrupted;}

  void set_value(T value_low, T value_high);
  virtual void calculate_contour();
  virtual resultStruct<T> collect();
};

template <typename T> 
class isoliner : public isobander<T> {
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
  void line_merge();

public:
  isoliner(T *x, int lenx, T *y, int leny, T *z, int nrow, int ncol, T value = 0.0) :
    isobander<T>(x, lenx, y, leny, z, nrow, ncol, value, T(0)) {}

  void set_value(T value);
  virtual void calculate_contour();
  virtual resultStruct<T> collect();
};

#endif // ISOBAND_H
