/*
 *  Orthomesh -- Orthogonal Delaunay Mesh Generator
 *
 *  Copyright (C) 2013  Clifford Wolf <clifford@clifford.at>
 *  
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *  
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#ifndef OMESH2D
#define OMESH2D

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <set>
#include <map>

#define OMESH2D_LEFT   0
#define OMESH2D_RIGHT  1
#define OMESH2D_TOP    2
#define OMESH2D_BOTTOM 3

#define OMESH2D_BITS 30
#define OMESH2D_FULL (int32_t(1) << OMESH2D_BITS)
#define OMESH2D_HALF (int32_t(1) << (OMESH2D_BITS-1))
#define OMESH2D_STEP exp2(-30)

// The 2D orthomesher
struct Omesh2d
{
	// set this to a value > 0 for debug output
	int config_verbose;

	// set this to set the 'aggressiveness' of the split/slice operations
	// 0 = only split half width slices
	// 1 = also split at groups of 3 quarter width slices
	// 2 = perform a split/slice operation whenever possible
	int config_split_slice_aggressiveness;

	// main API
	Omesh2d(int32_t major_width, int32_t major_height);
	virtual ~Omesh2d() { };
	void run();

	// refinement interface for the user
	virtual void refine(int32_t major_x, int32_t major_y, int32_t minor_x, int32_t minor_y, int32_t size_x, int32_t size_y, bool &refine_x, bool &refine_y) { };

	// coloring interface (used by svg writers)
	virtual std::string getcolor(int32_t major_x, int32_t major_y, std::vector<int32_t> &minor_xy) { return "gray"; };

	// data structures
	struct Coordinate;
	struct CellType;
	struct Segment;
	struct Geometry;

	// the tessellation
	int32_t grid_w, grid_h;
	std::map<Omesh2d::Coordinate, Omesh2d::CellType> grid;
	std::map<Omesh2d::CellType, Omesh2d::Geometry> geometries;

	// geometry construction
	bool create_geometry(const Omesh2d::CellType &ctype);
	void create_simple_split_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype);
	void create_simple_slice_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype);
	void create_split_slice_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype);

	// various helper functions
	bool run_slice_refine(const Omesh2d::Coordinate &coord, char axis, std::vector<int32_t> &slice_at, int32_t start, int shift);
	void propagate_into(std::set<Omesh2d::Coordinate> &prop_queue, std::vector<Omesh2d::Coordinate> &split_queue, const Omesh2d::Coordinate &coord);
	void find_neighbours_to_split(std::vector<Omesh2d::Coordinate> &neigh, const Omesh2d::Coordinate &coord);
	void add_neighbourhood_to_queue(std::set<Omesh2d::Coordinate> &neigh, const Omesh2d::Coordinate &coord);

	// various static helper functions
	static int count_trailing_zeros(int32_t v);
	static void fixup_slice_vector(std::vector<int32_t> &slice_at);
	static void add_slice_hierarchy(std::set<int32_t> &slices, int32_t s);
	static void create_reduced_slices_copy(std::vector<int32_t> &to, const std::vector<int32_t> &from);
	static bool check_slice_compatibility(const std::vector<int32_t> &sliceA, const std::vector<int32_t> &sliceB);
	static int32_t next_indent(int32_t last_indent, int32_t last_gap, int32_t this_gap, int32_t next_gap);
	static int32_t prev_indent(int32_t this_indent, int32_t last_gap, int32_t this_gap, int32_t next_gap);
	static bool calc_indent_chain(std::vector<int32_t> &indents, const std::vector<int32_t> &slices, int aggressiveness);
	static bool calc_indent_chain(std::vector<int32_t> &indents_a, std::vector<int32_t> &indents_b, const std::vector<int32_t> &slices_a, const std::vector<int32_t> &slices_b, int aggressiveness);
	static bool boundary_check(const std::pair<int32_t, int32_t> &p1, const std::pair<int32_t, int32_t> &p2, const Omesh2d::Segment &seg);
	static void triangle_circumcenter(double &x, double &y, double x1, double y1, double x2, double y2, double x3, double y3);
	static double triangle_area(double x1, double y1, double x2, double y2, double x3, double y3);
	static double polygon_area(const std::vector<double> &xy_data);

	// SVG and HTML output functions
	static void svg_write_delaunay(FILE *f, const Omesh2d::Geometry &geom);
	static void svg_write_voronoi(FILE *f, const Omesh2d::Geometry &geom);
	void svg_write_delaunay(FILE *f, const Omesh2d::Coordinate &coord);
	void svg_write_grid(FILE *f, double scale = 100, bool container = true);
	void html_write_celltype(FILE *f, const Omesh2d::CellType &ctype);
};

// This struct describes the boundary conditions for a cell
// The cell position and absolute size is not part of this description
struct Omesh2d::CellType
{
	// true if there is a split-point on the edge
	bool split_at[4];

	// coordinates (strictly incrasing) of the slice points
	// at the edge. coordinates are in ]0 .. 2^30[. the upper
	// left corner is (0, 0). The step between slices must
	// be a power of two and the step size must only change
	// by a factor of two from one step to the other. a cell
	// must not have a split_at[] and a slive_at[] rule on
	// the same side. opposite sides must match each other
	// at least for every second step.
	std::vector<int32_t> slice_at[4];

	CellType();
	void check() const;
	bool operator <(const Omesh2d::CellType &other) const;
	bool operator !=(const Omesh2d::CellType &other) const;
};

// This struct describes a segment of the mesh. each segment
// is just the convex hull of a list of points.
struct Omesh2d::Segment
{
	// coordinates are in [0 .. 2^30]. this list must be sorted.
	// points[]->first is the x and points[]->second the y coordinate.
	std::vector<std::pair<int32_t, int32_t>> points;

	Segment() { }
	bool operator <(const Omesh2d::Segment &other) const { return points < other.points; }
	bool operator !=(const Omesh2d::Segment &other) const { return points < other.points; }

	void check() const;
};

// This struct describes the geometry of a cell type. It is just
// a list of segments.
struct Omesh2d::Geometry
{
	// this list must be sorted.
	std::vector<Omesh2d::Segment> segments;

	Geometry() { }
	bool operator <(const Omesh2d::Geometry &other) const { return segments < other.segments; }

	void check() const;
};

// This struct describes the size and position of a cell
struct Omesh2d::Coordinate
{
	// coordinate on the major grid
	int32_t major_x, major_y;

	// coordinate on the minor (refined) grid within the major cell
	// values are in the range [0 .. 2^30[
	int32_t minor_x, minor_y;

	// size as fraction of a major grid cell
	// the value 2^30 is the width of a major grid cell
	// only powers of two are allowed here
	int32_t size;

	void check() const;
	void normalize();
	bool operator <(const Omesh2d::Coordinate &other) const;

	bool is_above(const Omesh2d::Coordinate &other) const;
	bool is_below(const Omesh2d::Coordinate &other) const;
	bool is_left_of(const Omesh2d::Coordinate &other) const;
	bool is_right_of(const Omesh2d::Coordinate &other) const;

	int32_t scale_from(const Omesh2d::Coordinate &other, int32_t v) const;
	int32_t scale_x_from(const Omesh2d::Coordinate &other, int32_t v) const;
	int32_t scale_y_from(const Omesh2d::Coordinate &other, int32_t v) const;

	std::string to_string() const;
};

#endif
