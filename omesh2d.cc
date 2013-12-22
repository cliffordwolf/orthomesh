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

#include "omesh2d.h"
#include <assert.h>
#include <algorithm>
#include <list>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
USING_PART_OF_NAMESPACE_EIGEN


/*************************************************************************
 *                               main API                                *
 *************************************************************************/

Omesh2d::Omesh2d(int32_t major_width, int32_t major_height)
{
	config_verbose = 0;
	config_split_slice_aggressiveness = 0;

	grid_w = major_width;
	grid_h = major_height;

	for (int32_t x = 0; x < grid_w; x++)
	for (int32_t y = 0; y < grid_h; y++) {
		Omesh2d::Coordinate coord;
		coord.major_x = x;
		coord.major_y = y;
		coord.minor_x = 0;
		coord.minor_y = 0;
		coord.size = OMESH2D_FULL;
		coord.check();
		grid[coord] = Omesh2d::CellType();
	}
}

void Omesh2d::run()
{
	if (config_verbose > 0)
		printf("Entering grid refinement loop.\n");

	std::set<Omesh2d::Coordinate> check_queue;
	std::set<Omesh2d::Coordinate> prop_queue;

	for (auto &it : grid)
		check_queue.insert(it.first);

	while (!check_queue.empty())
	{
		std::vector<Omesh2d::Coordinate> split_queue;

		for (auto &coord : check_queue)
		{
			if (grid.count(coord) == 0)
				continue;

			bool refine_x = false, refine_y = false;
			refine(coord.major_x, coord.major_y, coord.minor_x, coord.minor_y, coord.size, coord.size, refine_x, refine_y);

			if (refine_x && !refine_y) {
				std::vector<int32_t> slice_at;
				if (!run_slice_refine(coord, 'x', slice_at, 0, 1))
					goto resolve_conflict_by_split;
				fixup_slice_vector(slice_at);
				grid[coord].slice_at[OMESH2D_TOP] = slice_at;
				grid[coord].slice_at[OMESH2D_BOTTOM] = slice_at;
				add_neighbourhood_to_queue(prop_queue, coord);
				prop_queue.insert(coord);
			} else
			if (!refine_x && refine_y) {
				std::vector<int32_t> slice_at;
				if (!run_slice_refine(coord, 'y', slice_at, 0, 1))
					goto resolve_conflict_by_split;
				fixup_slice_vector(slice_at);
				grid[coord].slice_at[OMESH2D_LEFT] = slice_at;
				grid[coord].slice_at[OMESH2D_RIGHT] = slice_at;
				add_neighbourhood_to_queue(prop_queue, coord);
				prop_queue.insert(coord);
			} else
			if (refine_x && refine_y) {
		resolve_conflict_by_split:
				if (config_verbose > 3)
					printf("      putting cell %s on split queue (refine).\n", coord.to_string().c_str());
				split_queue.push_back(coord);
			} else {
				add_neighbourhood_to_queue(prop_queue, coord);
				prop_queue.insert(coord);
			}
		}

		check_queue.clear();

		while (split_queue.empty() && !prop_queue.empty())
		{
			std::set<Omesh2d::Coordinate> next_prop_queue;

			if (config_verbose > 1)
				printf("   processing grid cells from prop queue: %zd\n", prop_queue.size());

			for (auto &coord : prop_queue)
				if (grid.count(coord) != 0)
					propagate_into(next_prop_queue, split_queue, coord);

			prop_queue.swap(next_prop_queue);
		}

		while (!split_queue.empty())
		{
			for (auto &coord : split_queue)
			{
				if (grid.count(coord) == 0)
					continue;

				if (config_verbose > 3)
					printf("      splitting cell %s.\n", coord.to_string().c_str());

				for (int x = 0; x < 2; x++)
				for (int y = 0; y < 2; y++) {
					Omesh2d::Coordinate c = coord;
					c.size = c.size >> 1;
					c.minor_x += x*c.size;
					c.minor_y += y*c.size;
					c.check();
					check_queue.insert(c);
					grid[c] = Omesh2d::CellType();
				}

				grid.erase(coord);
			}

			std::vector<Omesh2d::Coordinate> next_split_queue;
			for (auto &coord : split_queue)
				find_neighbours_to_split(next_split_queue, coord);
			split_queue.swap(next_split_queue);
		}

		if (config_verbose > 1)
			printf("   grid cells after refinement iteration: %zd\n", grid.size());
	}
}


/*************************************************************************
 *                         geometry construction                         *
 *************************************************************************/

bool Omesh2d::create_geometry(const Omesh2d::CellType &ctype)
{
	if (geometries.count(ctype) == 0)
	{
		Omesh2d::Geometry &geom = geometries[ctype];

		int split_points = 0, slice_points = 0;
		for (int i = 0; i < 4; i++) {
			if (ctype.split_at[i]) split_points++;
			slice_points += ctype.slice_at[i].size();
		}

		if (config_verbose > 2)
			printf("      creating geometry for a cell type with %d split points and %d slice points.\n", split_points, slice_points);

		if (slice_points == 0)
			create_simple_split_geometry(geom, ctype);
		else if (split_points == 0)
			create_simple_slice_geometry(geom, ctype);
		else if (split_points == 1)
			create_split_slice_geometry(geom, ctype);

		if (!geom.segments.empty())
			geom.check();
	}

	return !geometries.at(ctype).segments.empty();
}

void Omesh2d::create_simple_split_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype)
{
	uint16_t id = 0;
	if (ctype.split_at[OMESH2D_LEFT])   id |= 0x1000;
	if (ctype.split_at[OMESH2D_RIGHT])  id |= 0x0100;
	if (ctype.split_at[OMESH2D_TOP])    id |= 0x0010;
	if (ctype.split_at[OMESH2D_BOTTOM]) id |= 0x0001;

	// The P? macros are named after the layout of a numerical keypad.
	// for example P1 is the lower left corner and P8 is the center point of the top edge.
	#define S  geom.segments.push_back(Omesh2d::Segment());
	#define P1 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(           0, OMESH2D_FULL));
	#define P2 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_HALF, OMESH2D_FULL));
	#define P3 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, OMESH2D_FULL));
	#define P4 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(           0, OMESH2D_HALF));
	#define P5 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_HALF, OMESH2D_HALF));
	#define P6 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, OMESH2D_HALF));
	#define P7 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(           0,            0));
	#define P8 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_HALF,            0));
	#define P9 geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL,            0));

	// Type segments as 'S172' or 'S7496' and use the following regex to create P? commands:
	//	s/S\([0-9]\)\([0-9]\)\([0-9]\)\([0-9]\)/S P\1 P\2 P\3 P\4/g
	//	s/S\([0-9]\)\([0-9]\)\([0-9]\)/S P\1 P\2 P\3/g

	switch (id)
	{
	case 0x0000:
		S P1 P3 P7 P9
		break;

	/* geometries with one split point */

	case 0x0001:
		S P7 P1 P2
		S P7 P9 P2
		S P2 P9 P3
		break;

	case 0x0010:
		S P1 P7 P8
		S P1 P8 P3
		S P8 P9 P3
		break;

	case 0x0100:
		S P7 P9 P6
		S P7 P1 P6
		S P1 P6 P3
		break;

	case 0x1000:
		S P7 P4 P9
		S P4 P9 P3
		S P4 P1 P3
		break;
	
	/* geometries with two opposite split points */

	case 0x0011:
		S P7 P1 P8 P2
		S P8 P2 P9 P3
		break;

	case 0x1100:
		S P7 P9 P4 P6
		S P4 P6 P1 P3
		break;

	/* geometries with two split points next to each other */

	case 0x1010:
		S P7 P4 P8
		S P8 P9 P3
		S P4 P8 P3
		S P1 P4 P3
		break;

	case 0x0110:
		S P1 P7 P8
		S P8 P6 P9
		S P1 P8 P6
		S P1 P6 P3
		break;

	case 0x0101:
		S P7 P1 P2
		S P7 P2 P6
		S P7 P6 P9
		S P2 P6 P3
		break;

	case 0x1001:
		S P9 P7 P4
		S P9 P4 P2
		S P9 P2 P3
		S P4 P2 P1
		break;

	/* geometries with three split points */

	case 0x1110:
		S P8 P7 P4
		S P8 P4 P6
		S P8 P6 P9
		S P1 P4 P6 P3
		break;

	case 0x0111:
		S P6 P9 P8
		S P6 P8 P2
		S P6 P2 P3
		S P7 P8 P1 P2
		break;

	case 0x1101:
		S P2 P1 P4
		S P2 P4 P6
		S P2 P6 P3
		S P7 P4 P9 P6
		break;

	case 0x1011:
		S P4 P7 P8
		S P4 P8 P2
		S P4 P2 P1
		S P8 P2 P9 P3
		break;

	default:
		abort();
	}

	#undef S
	#undef P1
	#undef P2
	#undef P3
	#undef P4
	#undef P5
	#undef P6
	#undef P7
	#undef P8
	#undef P9

	for (auto &seg : geom.segments)
		std::sort(seg.points.begin(), seg.points.end());
	std::sort(geom.segments.begin(), geom.segments.end());
}

void Omesh2d::create_simple_slice_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype)
{
	bool slices_x = !ctype.slice_at[OMESH2D_TOP].empty() || !ctype.slice_at[OMESH2D_BOTTOM].empty();
	bool slices_y = !ctype.slice_at[OMESH2D_LEFT].empty() || !ctype.slice_at[OMESH2D_RIGHT].empty();
	assert(slices_x != slices_y);

	std::vector<int32_t> slice_a = slices_x ? ctype.slice_at[OMESH2D_TOP] : ctype.slice_at[OMESH2D_LEFT];
	std::vector<int32_t> slice_b = slices_x ? ctype.slice_at[OMESH2D_BOTTOM] : ctype.slice_at[OMESH2D_RIGHT];

	int32_t last_pos = 0;
	slice_a.push_back(OMESH2D_FULL);
	slice_b.push_back(OMESH2D_FULL);

	size_t ia = 0, ib = 0;
	while (ia < slice_a.size() && ib < slice_b.size())
	{
		if (slice_a[ia] == slice_b[ib])
		{
			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, slice_a[ia]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, slice_a[ia]));
			last_pos = slice_a[ia], ia++, ib++;
			continue;
		}

		if (slice_a[ia] < slice_b[ib])
		{
			int32_t tmp = slice_a[ia++];
			assert(slice_a[ia] == slice_b[ib]);

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, slice_a[ia]));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, slice_a[ia]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, slice_a[ia]));

			last_pos = slice_a[ia], ia++, ib++;
			continue;
		}

		if (slice_a[ia] > slice_b[ib])
		{
			int32_t tmp = slice_b[ib++];
			assert(slice_a[ia] == slice_b[ib]);

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, slice_a[ia]));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, slice_a[ia]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, tmp));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, slice_a[ia]));

			last_pos = slice_a[ia], ia++, ib++;
			continue;
		}
	}

	assert(ia == slice_a.size());
	assert(ib == slice_b.size());

	if (slices_x)
		for (auto &seg : geom.segments)
		for (auto &p : seg.points) {
			int32_t tmp = p.first;
			p.first = p.second;
			p.second = tmp;
		}

	for (auto &seg : geom.segments)
		std::sort(seg.points.begin(), seg.points.end());
	std::sort(geom.segments.begin(), geom.segments.end());
}

void Omesh2d::create_split_slice_geometry(Omesh2d::Geometry &geom, const Omesh2d::CellType &ctype)
{
	bool slices_x = !ctype.slice_at[OMESH2D_TOP].empty() || !ctype.slice_at[OMESH2D_BOTTOM].empty();
	bool slices_y = !ctype.slice_at[OMESH2D_LEFT].empty() || !ctype.slice_at[OMESH2D_RIGHT].empty();
	assert(slices_x != slices_y);

	std::vector<int32_t> slice_a = slices_x ? ctype.slice_at[OMESH2D_TOP] : ctype.slice_at[OMESH2D_LEFT];
	std::vector<int32_t> slice_b = slices_x ? ctype.slice_at[OMESH2D_BOTTOM] : ctype.slice_at[OMESH2D_RIGHT];

	if (slices_x ? ctype.split_at[OMESH2D_RIGHT] : ctype.split_at[OMESH2D_BOTTOM]) {
		for (auto &s : slice_a) s = OMESH2D_FULL - s;
		for (auto &s : slice_b) s = OMESH2D_FULL - s;
		std::sort(slice_a.begin(), slice_a.end());
		std::sort(slice_b.begin(), slice_b.end());
	}

	std::vector<int32_t> indent_a, indent_b;
	if (!calc_indent_chain(indent_a, indent_b, slice_a, slice_b, config_split_slice_aggressiveness))
		return;

	int32_t last_pos = 0;
	slice_a.push_back(OMESH2D_FULL);
	slice_b.push_back(OMESH2D_FULL);

	int32_t last_indent_a = OMESH2D_HALF;
	int32_t last_indent_b = OMESH2D_HALF;
	indent_a.push_back(0);
	indent_b.push_back(0);

	size_t ia = 0, ib = 0;
	while (ia < slice_a.size() && ib < slice_b.size())
	{
		if (slice_a[ia] == slice_b[ib])
		{
			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, slice_a[ia]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(indent_a[ia], slice_a[ia]));

			if (last_indent_a != indent_a[ia]) {
				geom.segments.push_back(Omesh2d::Segment());
				geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(0, last_pos));
				geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(indent_a[ia], slice_a[ia]));
			}
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(last_indent_a, last_pos));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(last_indent_a, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(indent_a[ia], slice_a[ia]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL - last_indent_b, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL - indent_b[ib], slice_b[ib]));

			geom.segments.push_back(Omesh2d::Segment());
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, slice_b[ib]));
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL - indent_b[ib], slice_b[ib]));

			if (last_indent_b != indent_b[ib]) {
				geom.segments.push_back(Omesh2d::Segment());
				geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL, last_pos));
				geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL - indent_b[ib], slice_b[ib]));
			}
			geom.segments.back().points.push_back(std::pair<int32_t, int32_t>(OMESH2D_FULL - last_indent_b, last_pos));

			last_indent_a = indent_a[ia];
			last_indent_b = indent_b[ib];
			last_pos = slice_a[ia], ia++, ib++;
			continue;
		}

		if (slice_a[ia] < slice_b[ib])
		{
			// FIXME
			geom.segments.clear();
			return;
		}

		if (slice_a[ia] > slice_b[ib])
		{
			// FIXME
			geom.segments.clear();
			return;
		}
	}

	assert(ia == slice_a.size());
	assert(ib == slice_b.size());

	if (slices_x ? ctype.split_at[OMESH2D_RIGHT] : ctype.split_at[OMESH2D_BOTTOM])
		for (auto &seg : geom.segments)
		for (auto &p : seg.points)
			p.second = OMESH2D_FULL - p.second;

	if (slices_x)
		for (auto &seg : geom.segments)
		for (auto &p : seg.points) {
			int32_t tmp = p.first;
			p.first = p.second;
			p.second = tmp;
		}

	std::vector<Omesh2d::Segment> new_segments;
	for (auto &seg : geom.segments) {
		std::sort(seg.points.begin(), seg.points.end());
		std::vector<std::pair<int32_t, int32_t>> new_points;
		for (auto &p : seg.points)
			if (new_points.empty() || p != new_points.back())
				new_points.push_back(p);
		seg.points.swap(new_points);
		if (seg.points.size() >= 3)
			new_segments.push_back(seg);
	}
	std::sort(new_segments.begin(), new_segments.end());
	geom.segments.swap(new_segments);
}


/*************************************************************************
 *                      various helper functions                         *
 *************************************************************************/

bool Omesh2d::run_slice_refine(const Omesh2d::Coordinate &coord, char axis, std::vector<int32_t> &slice_at, int32_t start, int shift)
{
	bool refine_x = false, refine_y = false;
	int32_t minor_x = coord.minor_x, size_x = coord.size;
	int32_t minor_y = coord.minor_y, size_y = coord.size;

	if (axis == 'x') {
		size_x = size_x >> shift;
		minor_x += start;
	} else {
		size_y = size_y >> shift;
		minor_y += start;
	}

	refine(coord.major_x, coord.major_y, minor_x, minor_y, size_x, size_y, refine_x, refine_y);

	if (axis == 'x') {
		if (refine_y || (refine_x && !run_slice_refine(coord, axis, slice_at, start, shift+1)))
			return false;
		minor_x += size_x;
	} else {
		if (refine_x || (refine_y && !run_slice_refine(coord, axis, slice_at, start, shift+1)))
			return false;
		minor_y += size_y;
	}

	int32_t s = start + (coord.size >> shift);
	for (int32_t i = coord.size; i < OMESH2D_FULL; i = i << 1)
		s = s << 1;

	slice_at.push_back(s);

	refine(coord.major_x, coord.major_y, minor_x, minor_y, size_x, size_y, refine_x, refine_y);

	if (axis == 'x') {
		if (refine_y || (refine_x && !run_slice_refine(coord, axis, slice_at, start + size_x, shift+1)))
			return false;
	} else {
		if (refine_x || (refine_y && !run_slice_refine(coord, axis, slice_at, start + size_y, shift+1)))
			return false;
	}

	return true;
}

void Omesh2d::propagate_into(std::set<Omesh2d::Coordinate> &prop_queue, std::vector<Omesh2d::Coordinate> &split_queue, const Omesh2d::Coordinate &coord)
{
	Omesh2d::CellType ctype = grid.at(coord);

	std::set<Omesh2d::Coordinate> neigh;
	add_neighbourhood_to_queue(neigh, coord);

	for (auto &n : neigh)
	{
		const Omesh2d::CellType &nct = grid.at(n);
		bool is_above = n.is_above(coord);
		bool is_below = n.is_below(coord);
		bool is_left_of = n.is_left_of(coord);
		bool is_right_of = n.is_right_of(coord);

		if (n.size < coord.size)
		{
			if (is_above)    ctype.split_at[OMESH2D_TOP]    = true;
			if (is_below)    ctype.split_at[OMESH2D_BOTTOM] = true;
			if (is_left_of)  ctype.split_at[OMESH2D_LEFT]   = true;
			if (is_right_of) ctype.split_at[OMESH2D_RIGHT]  = true;
		}

		if (is_above || is_below)
		{
			int32_t last_step = 0;
			if (!nct.slice_at[OMESH2D_LEFT].empty()) {
				int32_t sl = is_above ? OMESH2D_FULL - nct.slice_at[OMESH2D_LEFT].back() : nct.slice_at[OMESH2D_LEFT].front();
				last_step = std::max(last_step, coord.scale_from(n, sl));
			}
			if (!nct.slice_at[OMESH2D_RIGHT].empty()) {
				int32_t sl = is_above ? OMESH2D_FULL - nct.slice_at[OMESH2D_RIGHT].back() : nct.slice_at[OMESH2D_RIGHT].front();
				last_step = std::max(last_step, coord.scale_from(n, sl));
			}

			if (0 < last_step && last_step < OMESH2D_HALF) {
				int32_t s = is_above ? (last_step << 1) : OMESH2D_FULL - (last_step << 1);
				ctype.slice_at[OMESH2D_LEFT].push_back(s);
				ctype.slice_at[OMESH2D_RIGHT].push_back(s);
			}

			int from_slice_at = is_above ? OMESH2D_BOTTOM : OMESH2D_TOP;
			int to_slice_at = is_above ? OMESH2D_TOP : OMESH2D_BOTTOM;

			for (auto s : nct.slice_at[from_slice_at])
				ctype.slice_at[to_slice_at].push_back(coord.scale_x_from(n, s));

			if (nct.split_at[from_slice_at])
				ctype.slice_at[to_slice_at].push_back(coord.scale_x_from(n, OMESH2D_HALF));
		}
		else
		if (is_left_of || is_right_of)
		{
			int32_t last_step = 0;
			if (!nct.slice_at[OMESH2D_TOP].empty()) {
				int32_t sl = is_left_of ? OMESH2D_FULL - nct.slice_at[OMESH2D_TOP].back() : nct.slice_at[OMESH2D_TOP].front();
				last_step = std::max(last_step, coord.scale_from(n, sl));
			}
			if (!nct.slice_at[OMESH2D_BOTTOM].empty()) {
				int32_t sl = is_left_of ? OMESH2D_FULL - nct.slice_at[OMESH2D_BOTTOM].back() : nct.slice_at[OMESH2D_BOTTOM].front();
				last_step = std::max(last_step, coord.scale_from(n, sl));
			}

			if (0 < last_step && last_step < OMESH2D_HALF) {
				int32_t s = is_left_of ? (last_step << 1) : OMESH2D_FULL - (last_step << 1);
				ctype.slice_at[OMESH2D_TOP].push_back(s);
				ctype.slice_at[OMESH2D_BOTTOM].push_back(s);
			}

			int from_slice_at = is_left_of ? OMESH2D_RIGHT : OMESH2D_LEFT;
			int to_slice_at = is_left_of ? OMESH2D_LEFT : OMESH2D_RIGHT;

			for (auto s : nct.slice_at[from_slice_at])
				ctype.slice_at[to_slice_at].push_back(coord.scale_y_from(n, s));

			if (nct.split_at[from_slice_at])
				ctype.slice_at[to_slice_at].push_back(coord.scale_y_from(n, OMESH2D_HALF));
		}
		else
			abort();
	}

	if (ctype.split_at[OMESH2D_TOP] && ctype.split_at[OMESH2D_BOTTOM] && ctype.split_at[OMESH2D_LEFT] && ctype.split_at[OMESH2D_RIGHT]) {
		split_queue.push_back(coord);
		return;
	}

	for (int i = 0; i < 4; i++)
		fixup_slice_vector(ctype.slice_at[i]);

	for (int i = 0; i < 4; i++)
		if (ctype.slice_at[i].size() == 1 && ctype.slice_at[i ^ 1].size() <= 1) {
			ctype.slice_at[i].clear();
			ctype.split_at[i] = true;
		}

	for (int i = 0; i < 4; i++) {
		for (auto s : ctype.slice_at[i])
			if (s == OMESH2D_HALF) ctype.split_at[i] = false;
	}

	bool slices_x = !ctype.slice_at[OMESH2D_TOP].empty() || !ctype.slice_at[OMESH2D_BOTTOM].empty();
	bool slices_y = !ctype.slice_at[OMESH2D_LEFT].empty() || !ctype.slice_at[OMESH2D_RIGHT].empty();

	if (ctype.split_at[OMESH2D_TOP] && ctype.split_at[OMESH2D_BOTTOM]) slices_x = true;
	if (ctype.split_at[OMESH2D_LEFT] && ctype.split_at[OMESH2D_RIGHT]) slices_y = true;

	if ((slices_x && slices_y) || !check_slice_compatibility(ctype.slice_at[OMESH2D_TOP], ctype.slice_at[OMESH2D_BOTTOM]) ||
			!check_slice_compatibility(ctype.slice_at[OMESH2D_LEFT], ctype.slice_at[OMESH2D_RIGHT])) {
		if (config_verbose > 3)
			printf("      putting cell %s on split queue (incompatible slices).\n", coord.to_string().c_str());
		split_queue.push_back(coord);
		return;
	}

	if (!create_geometry(ctype)) {
		if (config_verbose > 3)
			printf("      putting cell %s on split queue (no valid geometry).\n", coord.to_string().c_str());
		split_queue.push_back(coord);
		return;
	}

	if (ctype != grid.at(coord)) {
		grid[coord] = ctype;
		prop_queue.insert(neigh.begin(), neigh.end());
		return;
	}
}

// we are splitting at coord. if this cell has any neighbours larger then itself, we also need
// to split those. The list of neighbours that need splitting is returned in &neigh.
void Omesh2d::find_neighbours_to_split(std::vector<Omesh2d::Coordinate> &neigh, const Omesh2d::Coordinate &coord)
{
	// list of candidates as pairs of x/y offsets (unit length is coord.size)
	static const int neigh_map[][2] = {
		// left
		{ -2, -1 }, { -2,  0 },
		// right
		{  1, -1 }, {  1,  0 },
		// top
		{ -1, -2 }, {  0, -2 },
		// bottom
		{ -1,  1 }, {  0,  1 }
	};

	for (size_t i = 0; i < sizeof(neigh_map)/sizeof(neigh_map[0]); i++) {
		Omesh2d::Coordinate n = coord;
		n.minor_x += neigh_map[i][0] * n.size;
		n.minor_y += neigh_map[i][1] * n.size;
		n.normalize();
		n.size *= 2;
		if (grid.count(n) > 0) {
			if (config_verbose > 3)
				printf("      putting cell %s on split queue (neighbour).\n", n.to_string().c_str());
			neigh.push_back(n);
		}
	}
}

// add ALL the neighbours
void Omesh2d::add_neighbourhood_to_queue(std::set<Omesh2d::Coordinate> &neigh, const Omesh2d::Coordinate &coord)
{
	// list of candidates as tuples of x_off, y_off, size (unit length is half coord.size)
	static const int neigh_map[][3] =
	{
		// double size
			/* left   */ { -4, -2,  4 }, { -4,  0,  4 },
			/* right  */ {  2, -2,  4 }, {  2,  0,  4 },
			/* top    */ { -2, -4,  4 }, {  0, -4,  4 },
			/* bottom */ { -2,  2,  4 }, {  0,  2,  4 },

		// same size
			/* left   */ { -2,  0,  2 },
			/* right  */ {  2,  0,  2 },
			/* top    */ {  0, -2,  2 },
			/* bottom */ {  0,  2,  2 },

		// half size
			/* left   */ { -1,  0,  1 }, { -1,  1,  1 },
			/* right  */ {  2,  0,  1 }, {  2,  1,  1 },
			/* top    */ {  0, -1,  1 }, {  1, -1,  1 },
			/* bottom */ {  0,  2,  1 }, {  1,  2,  1 }
	};

	for (size_t i = 0; i < sizeof(neigh_map)/sizeof(neigh_map[0]); i++) {
		Omesh2d::Coordinate n = coord;
		n.minor_x += neigh_map[i][0] * (n.size >> 1);
		n.minor_y += neigh_map[i][1] * (n.size >> 1);
		n.size = neigh_map[i][2] * (n.size >> 1);
		n.normalize();
		if (grid.count(n) > 0)
			neigh.insert(n);
	}
}


/*************************************************************************
 *                    various static helper functions                    *
 *************************************************************************/

int Omesh2d::count_trailing_zeros(int32_t v)
{
	if (v == 0)
		return 32;
	int zeros = 0;
	for (int i = 16; i != 0; i = i >> 1)
		if ((v & ((1 << i) - 1)) == 0) {
			zeros += i;
			v = v >> i;
		}
	return zeros;
}

void Omesh2d::fixup_slice_vector(std::vector<int32_t> &slice_at)
{
	// STEP 1: Make sure the slice hierarchy is correct

	std::set<int32_t> slices;

	for (auto s : slice_at)
		add_slice_hierarchy(slices, s);

	slice_at.clear();
	slice_at.insert(slice_at.end(), slices.begin(), slices.end());

	// STEP 2: Make sure we don't change step sizes by more than a factor two

	std::list<int32_t> steps;

	for (size_t i = 0; i <= slice_at.size(); i++) {
		int32_t from = i > 0 ? slice_at[i-1] : 0;
		int32_t to = i < slice_at.size() ? slice_at[i] : OMESH2D_FULL;
		if (to != from)
			steps.push_back(to - from);
	}

	int32_t position = 0;
	for (std::list<int32_t>::iterator it = steps.begin(); it != steps.end(); it++)
	{
		if (it != steps.begin()) {
			std::list<int32_t>::iterator prev = it;
			int32_t max_step = *(--prev);
			if (count_trailing_zeros(position + 2*max_step) > count_trailing_zeros(position + max_step))
				max_step += max_step;
			if (*it > max_step) {
				position += max_step;
				steps.insert(it, max_step);
				*(it--) -= max_step;
				continue;
			}
		}

		std::list<int32_t>::iterator next = it;
		if (++next != steps.end()) {
			int32_t max_step = *next;
			if (count_trailing_zeros(position + *it - 2*max_step) > count_trailing_zeros(position + *it - max_step))
				max_step += max_step;
			if (*it > max_step) {
				steps.insert(next, max_step);
				*(it--) -= max_step;
				continue;
			}
		}

		position += *it;
	}

	slice_at.clear();
	for (auto it : steps)
		slice_at.push_back(slice_at.empty() ? it : slice_at.back() + it);
	assert(slice_at.back() == OMESH2D_FULL);
	slice_at.pop_back();
}

void Omesh2d::add_slice_hierarchy(std::set<int32_t> &slices, int32_t s)
{
	int32_t p = 0;
	int32_t step = OMESH2D_FULL;

	if (s == 0 || s == OMESH2D_FULL) {
		slices.insert(s);
		return;
	}

	while (p != s) {
		step = step >> 1;
		p += s < p ? -step : +step;
		slices.insert(p);
	}
}

void Omesh2d::create_reduced_slices_copy(std::vector<int32_t> &to, const std::vector<int32_t> &from)
{
	to.clear();
	for (size_t i = 0; i < from.size(); i++)
	{
		int32_t prev = i == 0 ? 0 : from[i-1];
		int32_t next = i+1 < from.size() ? from[i+1] : OMESH2D_FULL;
		int this_trailing_zeros = count_trailing_zeros(from[i]);

		if (count_trailing_zeros(prev) < this_trailing_zeros || count_trailing_zeros(next) < this_trailing_zeros)
			to.push_back(from[i]);
	}
}

bool Omesh2d::check_slice_compatibility(const std::vector<int32_t> &sliceA, const std::vector<int32_t> &sliceB)
{
	std::vector<int32_t> reducedA, reducedB;
	create_reduced_slices_copy(reducedA, sliceA);
	create_reduced_slices_copy(reducedB, sliceB);

	size_t i = 0;
	for (auto s : reducedA) {
		while (i < sliceB.size() && sliceB[i] < s) i++;
		if (sliceB[i] != s)
			return false;
	}

	i = 0;
	for (auto s : reducedB) {
		while (i < sliceA.size() && sliceA[i] < s) i++;
		if (sliceA[i] != s)
			return false;
	}

	return true;
}

int32_t Omesh2d::next_indent(int32_t last_indent, int32_t last_gap, int32_t this_gap, int32_t next_gap)
{
	// l .. width of half cell
	// u .. height of this gap
	// a .. horiz. distance of upper point from center line
	// b .. horiz. distance of lower point from center line
	// h .. max vert. distance of circumcenter of gap center

	// Maxima Fomulas:
	//   solve((a+b)*(b-a)/(2*u) = h, b);
	//   solve((l-b)*(b-a)/(2*u) = h, b);

	double l = OMESH2D_HALF;
	double u = this_gap;
	double a = OMESH2D_HALF - last_indent;

	// case 1: limit from below
	double h1 = (this_gap + next_gap) / 2.0;
	double b1 = sqrt(a*a + 2*h1*u);

	// case 2: limit from above
	double h2 = (this_gap + last_gap) / 2.0;
	double b2 = (l + a - sqrt(a*a - 2*a*l + l*l - 8*h2*u)) / 2.0;

	return round(std::max(std::max(l - b1, l - b2), 0.0));
}

int32_t Omesh2d::prev_indent(int32_t this_indent, int32_t last_gap, int32_t this_gap, int32_t next_gap)
{
	double l = OMESH2D_HALF;
	double u = this_gap;
	double b = OMESH2D_HALF - this_indent;

	double h1 = (this_gap + next_gap) / 2.0;
	double a1 = sqrt(b*b - 2*h1*u);

	double h2 = (this_gap + last_gap) / 2.0;
	double a2 = (2*h2*u - b*l + b*b) / (b-l);

	return round(std::max(std::max(l - a1, l - a2), 0.0));
}

bool Omesh2d::calc_indent_chain(std::vector<int32_t> &indents, const std::vector<int32_t> &slices, int aggressiveness)
{
	// the trivial cases
	if (slices.size() > 0 && slices.at(0) >= OMESH2D_HALF) {
		indents.insert(indents.end(), slices.size(), 0);
		return true;
	}
	if (slices.size() > 1 && slices.at(slices.size()-2) <= OMESH2D_HALF) {
		indents.insert(indents.end(), slices.size()-1, OMESH2D_HALF);
		indents.push_back(0);
		return true;
	}

	if (aggressiveness == 0)
		return false;

	int32_t last_pos = 0, last_gap = 0, last_indent = OMESH2D_HALF;
	for (size_t i = 0; i < slices.size(); i++)
	{
		int32_t this_gap = slices[i] - last_pos;
		int32_t next_gap = i+1 < slices.size() ? slices[i+1] - slices[i] : 0;
		int32_t indent = last_indent;
		if (aggressiveness == 2 || this_gap >= OMESH2D_FULL / 4) {
			indent = next_indent(last_indent, last_gap * 0.9, this_gap, next_gap * (aggressiveness == 2 ? 1.0 : 1.1));
			if (indent == 0)
				indent = next_indent(last_indent, last_gap * 0.9, this_gap, next_gap * 0.9);
		}
		if (indent == 0 && last_indent != 0 && !indents.empty()) {
			int32_t upd_last_indent = prev_indent(indent, last_gap * 0.9, this_gap, next_gap * 0.9);
			if (upd_last_indent > last_indent)
				indents.back() = upd_last_indent;
		}
		indents.push_back(indent);
		last_pos = slices[i];
		last_gap = this_gap;
		last_indent = indent;
	}
	return indents.back() == 0;
}

bool Omesh2d::calc_indent_chain(std::vector<int32_t> &indents_a, std::vector<int32_t> &indents_b, const std::vector<int32_t> &slices_a, const std::vector<int32_t> &slices_b, int aggressiveness)
{
	std::vector<int32_t> slices(slices_a.size() + slices_b.size());
	std::merge(slices_a.begin(), slices_a.end(), slices_b.begin(), slices_b.end(), slices.begin());
	slices.erase(std::unique(slices.begin(), slices.end()), slices.end());
	slices.push_back(OMESH2D_FULL);

	std::vector<int32_t> indents;
	if (!calc_indent_chain(indents, slices, std::min(aggressiveness, 1))) {
		if (aggressiveness < 2 || !calc_indent_chain(indents, slices, 2))
			return false;
	}

	assert(indents.back() == 0);
	assert(indents.size() == slices.size());

	for (size_t i = 0, ia = 0, ib = 0; i < slices.size(); i++) {
		while (ia < slices_a.size() && slices_a[ia] < slices[i]) ia++;
		while (ib < slices_b.size() && slices_b[ib] < slices[i]) ib++;
		if (ia < slices_a.size() && slices_a[ia] == slices[i]) indents_a.push_back(indents[i]);
		if (ib < slices_b.size() && slices_b[ib] == slices[i]) indents_b.push_back(indents[i]);
	}

	assert(indents_a.size() == slices_a.size());
	assert(indents_b.size() == slices_b.size());
	return true;
}

bool Omesh2d::boundary_check(const std::pair<int32_t, int32_t> &p1, const std::pair<int32_t, int32_t> &p2, const Omesh2d::Segment &seg)
{
	bool found_negative = false, found_positive = false;
	int64_t k1 = p1.second - p2.second, k2 = p2.first - p1.first;

	for (const auto &p : seg.points) {
		int64_t v = (p.first - p1.first) * k1 + (p.second - p1.second) * k2;
		if (v < 0) found_negative = true;
		if (v > 0) found_positive = true;
	}

	return found_negative != found_positive;
}

void Omesh2d::triangle_circumcenter(double &x, double &y, double x1, double y1, double x2, double y2, double x3, double y3)
{
	Vector2f A(x1, y1), B(x2, y2), C(x3, y3);
	Vector2f AB = B - A, AC = C - A;
	Vector2f ABx(AB[1], -AB[0]);
	Vector2f ACx(AC[1], -AC[0]);

	Matrix2f M;
	M << ABx[0], -ACx[0],
	     ABx[1], -ACx[1];

	Vector2f r = 0.5*AC - 0.5*AB;
	Vector2f p;

	if (!M.lu().solve(r, &p)) {
		fprintf(stderr, "*** not really a triangle: %f,%f %f,%f %f,%f\n", x1, y1, x2, y2, x3, y3);
		abort();
	}

	Vector2f X = A + 0.5*AB + ABx*p[0];
	x = X[0], y = X[1];
}

double Omesh2d::triangle_area(double x1, double y1, double x2, double y2, double x3, double y3)
{
	Vector2f A(x1, y1), B(x2, y2), C(x3, y3);

	// project C onto AB
	Vector2f X = A + (B-A) * ((C-A).dot(B-A) / (B-A).squaredNorm());

	// calculate area with height on C
	return ((C-X).norm() * (A-B).norm()) / 2.0;
}

double Omesh2d::polygon_area(const std::vector<double> &xy_data)
{
	Vector2f C(0, 0);
	std::vector<Vector2f> Q;

	for (size_t i = 0; i < xy_data.size(); i += 2) {
		C += Vector2f(xy_data[i], xy_data[i+1]) / (xy_data.size()/2);
		Q.push_back(Vector2f(xy_data[i], xy_data[i+1]));
	}

	std::sort(Q.begin(), Q.end(), [&](const Vector2f &A, const Vector2f &B) -> bool { return atan2(A[1] - C[1], A[0] - C[0]) < atan2(B[1] - C[1], B[0] - C[0]); });

	double area = 0;
	for (size_t i = 1; i+1 < Q.size(); i++)
		area += triangle_area(Q[0][0], Q[0][1], Q[i][0], Q[i][1], Q[i+1][0], Q[i+1][1]);
	return area;
}


/*************************************************************************
 *                    SVG and HTML output functions                      *
 *************************************************************************/

void Omesh2d::svg_write_delaunay(FILE *f, const Omesh2d::Geometry &geom)
{
	for (const auto &seg : geom.segments)
	{
		Vector2f C(0, 0);
		std::vector<Vector2f> Q;
		for (auto &p : seg.points) {
			C += Vector2f(p.first * OMESH2D_STEP, p.second * OMESH2D_STEP) / seg.points.size();
			Q.push_back(Vector2f(p.first * OMESH2D_STEP, p.second * OMESH2D_STEP));
		}

		std::sort(Q.begin(), Q.end(), [&](const Vector2f &A, const Vector2f &B) -> bool { return atan2(A[1] - C[1], A[0] - C[0]) < atan2(B[1] - C[1], B[0] - C[0]); });

		fprintf(f, "<polygon points='");
		for (size_t i = 0; i < Q.size(); i++)
			fprintf(f, "%s%f,%f", i ? " " : "", Q[i][0], Q[i][1]);
		fprintf(f, "'/>\n");
	}
}

void Omesh2d::svg_write_voronoi(FILE *f, const Omesh2d::Geometry &geom)
{
	for (const auto &seg : geom.segments)
	{
		assert(seg.points.size() >= 3);

		double cent_x, cent_y;
		triangle_circumcenter(cent_x, cent_y, seg.points[0].first, seg.points[0].second,
				seg.points[1].first, seg.points[1].second, seg.points[2].first, seg.points[2].second);

		for (const auto p1 : seg.points)
		for (const auto p2 : seg.points)
		{
			if (p1 >= p2 || !boundary_check(p1, p2, seg))
				continue;

			int32_t x1 = cent_x, y1 = cent_y;
			int32_t x2 = (p1.first >> 1) + (p2.first >> 1);
			int32_t y2 = (p1.second >> 1) + (p2.second >> 1);

			fprintf(f, "<line x1='%f' y1='%f' x2='%f' y2='%f' />\n", x1 * OMESH2D_STEP, y1 * OMESH2D_STEP, x2 * OMESH2D_STEP, y2 * OMESH2D_STEP);
		}
	}
}

void Omesh2d::svg_write_delaunay(FILE *f, const Omesh2d::Coordinate &coord)
{
	Omesh2d::Geometry &geom = geometries.at(grid.at(coord));
	int shift = OMESH2D_BITS - count_trailing_zeros(coord.size);

	for (const auto &seg : geom.segments)
	{
		Vector2f C(0, 0);
		std::vector<Vector2f> Q;
		for (auto &p : seg.points) {
			C += Vector2f(p.first * OMESH2D_STEP, p.second * OMESH2D_STEP) / seg.points.size();
			Q.push_back(Vector2f(p.first * OMESH2D_STEP, p.second * OMESH2D_STEP));
		}

		std::sort(Q.begin(), Q.end(), [&](const Vector2f &A, const Vector2f &B) -> bool { return atan2(A[1] - C[1], A[0] - C[0]) < atan2(B[1] - C[1], B[0] - C[0]); });

		std::vector<int32_t> minor_xy;
		for (auto &p : seg.points) {
			minor_xy.push_back((p.first >> shift) + coord.minor_x);
			minor_xy.push_back((p.second >> shift) + coord.minor_y);
		}

		fprintf(f, "<polygon fill='%s' transform='translate(%f, %f) scale(%f, %f)' points='",
				getcolor(coord.major_x, coord.major_y, minor_xy).c_str(),
				coord.major_x + coord.minor_x * OMESH2D_STEP,
				coord.major_y + coord.minor_y * OMESH2D_STEP,
				coord.size * OMESH2D_STEP, coord.size * OMESH2D_STEP);
		for (size_t i = 0; i < Q.size(); i++)
			fprintf(f, "%s%f,%f", i ? " " : "", Q[i][0], Q[i][1]);
		fprintf(f, "'/>\n");
	}
}

void Omesh2d::svg_write_grid(FILE *f, double scale, bool container)
{
	if (container) {
		fprintf(f, "<?xml version='1.0' encoding='UTF-8'?>\n");
		fprintf(f, "<svg xmlns='http://www.w3.org/2000/svg' version='1.2' baseProfile='tiny' width='%f' height='%f'>\n", scale*grid_w + 20, scale*grid_h + 20);
		fprintf(f, "<g transform='translate(10, 10)'>\n");
	}

	for (auto &it : grid)
	{
		auto &coord = it.first;
		fprintf(f, "<g stroke='black' stroke-width='%f' transform='scale(%f, %f)'>\n",
				(0.5 / scale) / (it.first.size * OMESH2D_STEP), scale, scale);
		svg_write_delaunay(f, coord);
		fprintf(f, "</g>\n");
	}

	for (int32_t x = 0; x <= grid_w; x++)
		fprintf(f, "<line stroke='darkred' stroke-width='0.5' x1='%f' y1='0' x2='%f' y2='%f' />\n", scale*x, scale*x, scale*grid_h);
	for (int32_t y = 0; y <= grid_h; y++)
		fprintf(f, "<line stroke='darkred' stroke-width='0.5' x1='0' y1='%f' x2='%f' y2='%f' />\n", scale*y, scale*grid_w, scale*y);

#if 0
	for (auto &it : grid)
	{
		auto &coord = it.first;
		auto &geom = geometries.at(it.second);

		fprintf(f, "<g fill='none' stroke='blue' stroke-width='%f' transform='scale(%f, %f) translate(%f, %f) scale(%f, %f)'>\n",
				(0.2 / scale) / (coord.size * OMESH2D_STEP), scale, scale,
				coord.major_x + coord.minor_x * OMESH2D_STEP, coord.major_y + coord.minor_y * OMESH2D_STEP,
				coord.size * OMESH2D_STEP, coord.size * OMESH2D_STEP);
		svg_write_voronoi(f, geom);
		fprintf(f, "</g>\n");
	}
#endif

	if (container)
		fprintf(f, "</g></svg>\n");
}

void Omesh2d::html_write_celltype(FILE *f, const Omesh2d::CellType &ctype)
{
	fprintf(f, "<tr><td valign='top'>\n");

	int count_split = 0, count_slice = 0;
	for (int i = 0; i < 4; i ++) {
		count_split += ctype.split_at[i] ? 1 : 0;
		count_slice += ctype.slice_at[i].size();
	}

	fprintf(f, "<b>CellType with %d split- and %d slice-constraints: %zd segments</b><br/>",
			count_split, count_slice, geometries.at(ctype).segments.size());

	static const char *edge_names[] = {
		"left", "right", "top", "bottom"
	};

	for (int i = 0; i < 4; i ++)
	{
		if (ctype.split_at[i])
			fprintf(f, "Split at the %s edge.<br/>\n", edge_names[i]);

		if (!ctype.slice_at[i].empty())
		{
			fprintf(f, "Slice at the %s edge:", edge_names[i]);
			for (auto s : ctype.slice_at[i]) {
				int shift = OMESH2D_BITS - count_trailing_zeros(s);
				fprintf(f, " %d/%d", s / (OMESH2D_FULL >> shift), 1 << shift);
			}
			fprintf(f, "<br/>\n");
		}
	}

	fprintf(f, "</td><td>\n");

	fprintf(f, "<svg width='300' height='300'>\n");
	fprintf(f, "<g stroke-width='0.005' transform='translate(10, 10) scale(280, 280)'>\n");
	fprintf(f, "<g fill='gray' stroke='black'>\n");
	svg_write_delaunay(f, geometries.at(ctype));
	fprintf(f, "</g></g></svg>\n");

	fprintf(f, "</td><td>\n");

	fprintf(f, "<svg width='300' height='300'>\n");
	fprintf(f, "<g stroke-width='0.005' transform='translate(10, 10) scale(280, 280)'>\n");
	fprintf(f, "<g fill='gray' stroke='black'>\n");
	svg_write_delaunay(f, geometries.at(ctype));
	fprintf(f, "</g><g stroke='blue'>\n");
	svg_write_voronoi(f, geometries.at(ctype));
	fprintf(f, "</g></g></svg>\n");

	fprintf(f, "</td></tr>\n");
}


/*************************************************************************
 *                            data structures                            *
 *************************************************************************/

Omesh2d::CellType::CellType()
{
	for (int i = 0; i < 4; i++)
		split_at[i] = false;
}

void Omesh2d::CellType::check() const
{
#ifndef NDEBUG
	for (int i = 0; i < 4; i++) {
		int32_t last_s = 0;
		for (auto s : slice_at[i]) {
			assert(last_s < s);
			last_s = s;
		}
		assert(last_s < OMESH2D_FULL);
	}
#endif
}

bool Omesh2d::CellType::operator<(const Omesh2d::CellType &other) const
{
#ifndef NDEBUG
	check();
	other.check();
#endif
	
	for (int i = 0; i < 4; i++) {
		if (split_at[i] != other.split_at[i])
			return split_at[i] < other.split_at[i];
		if (slice_at[i] != other.slice_at[i])
			return slice_at[i] < other.slice_at[i];
	}
	return false;
}

bool Omesh2d::CellType::operator!=(const Omesh2d::CellType &other) const
{
	for (int i = 0; i < 4; i++) {
		if (split_at[i] != other.split_at[i])
			return true;
		if (slice_at[i] != other.slice_at[i])
			return true;
	}
	return false;
}

void Omesh2d::Segment::check() const
{
#ifndef NDEBUG
	std::pair<int32_t, int32_t> last_point(-1, -1);
	assert(points.size() >= 3);
	for (auto &p : points) {
		assert(last_point < p);
		assert(0 <= p.first);
		assert(0 <= p.second);
		assert(OMESH2D_FULL >= p.first);
		assert(OMESH2D_FULL >= p.second);
		last_point = p;
	}

	// all points must sit on the same circumcircle
	double ref_x = -1, ref_y = -1;
	for (size_t i = 0; i < points.size(); i++)
	for (size_t j = i+1; j < points.size(); j++)
	for (size_t k = j+1; k < points.size(); k++)
	{
		auto &p1 = points[i];
		auto &p2 = points[j];
		auto &p3 = points[k];

		double this_x, this_y;
		Omesh2d::triangle_circumcenter(this_x, this_y, p1.first * OMESH2D_STEP, p1.second * OMESH2D_STEP,
				p2.first * OMESH2D_STEP, p2.second * OMESH2D_STEP, p3.first * OMESH2D_STEP, p3.second * OMESH2D_STEP);

		if (i == 0 && j == 1 && k == 2)
			ref_x = this_x, ref_y = this_y;
		else
			assert(fabs(ref_x - this_x) < 1e-6 && fabs(ref_y - this_y) < 1e-6);
	}

	// the circumcenter must be within the cell borders
	assert(0.0 <= ref_x && ref_x <= 1.0);
	assert(0.0 <= ref_y && ref_y <= 1.0);
#endif
}

void Omesh2d::Geometry::check() const
{
#ifndef NDEBUG
	double total_area = 0;
	for (auto &seg : segments) {
		seg.check();
		std::vector<double> xy_data;
		for (auto &p : seg.points) {
			xy_data.push_back(p.first * OMESH2D_STEP);
			xy_data.push_back(p.second * OMESH2D_STEP);
		}
		total_area += Omesh2d::polygon_area(xy_data);
	}
	assert(fabs(total_area - 1.0) < 1e-6);

	// check for delaunay within cell
	for (size_t i = 0; i < segments.size(); i++)
	for (auto &seg : segments)
	{
		double x, y, r;
		Omesh2d::triangle_circumcenter(x, y,
				seg.points[0].first * OMESH2D_STEP, seg.points[0].second * OMESH2D_STEP,
				seg.points[1].first * OMESH2D_STEP, seg.points[1].second * OMESH2D_STEP,
				seg.points[2].first * OMESH2D_STEP, seg.points[2].second * OMESH2D_STEP);
		r = sqrt(pow(seg.points[0].first * OMESH2D_STEP - x, 2) + pow(seg.points[0].second * OMESH2D_STEP - y, 2));

		std::set<std::pair<int32_t, int32_t>> sp;
		sp.insert(seg.points.begin(), seg.points.end());

		for (auto &s : segments)
		for (auto &p : s.points)
			if (sp.count(p) == 0)
				assert(sqrt(pow(p.first * OMESH2D_STEP - x, 2) + pow(p.second * OMESH2D_STEP - y, 2)) > r);
	}
#endif
}

void Omesh2d::Coordinate::check() const
{
	assert(0 < size);
	assert(0 <= major_x);
	assert(0 <= major_y);
	assert(0 <= minor_x);
	assert(0 <= minor_y);

	assert(0 <= OMESH2D_FULL);
	assert(minor_x + size <= OMESH2D_FULL);
	assert(minor_y + size <= OMESH2D_FULL);

	assert(minor_x % size == 0);
	assert(minor_y % size == 0);
}

void Omesh2d::Coordinate::normalize()
{
	if (minor_x < 0) {
		minor_x += OMESH2D_FULL;
		major_x--;
	}
	if (minor_y < 0) {
		minor_y += OMESH2D_FULL;
		major_y--;
	}
	if (minor_x >= OMESH2D_FULL) {
		minor_x -= OMESH2D_FULL;
		major_x++;
	}
	if (minor_y >= OMESH2D_FULL) {
		minor_y -= OMESH2D_FULL;
		major_y++;
	}
}

bool Omesh2d::Coordinate::operator<(const Omesh2d::Coordinate &other) const
{
	if (major_x != other.major_x)
		return major_x < other.major_x;
	if (major_y != other.major_y)
		return major_y < other.major_y;
	if (minor_x != other.minor_x)
		return minor_x < other.minor_x;
	if (minor_y != other.minor_y)
		return minor_y < other.minor_y;
	if (size != other.size)
		return size < other.size;
	return false;
}

bool Omesh2d::Coordinate::is_above(const Omesh2d::Coordinate &other) const
{
	Omesh2d::Coordinate tmp = *this;
	tmp.minor_y += tmp.size;
	tmp.normalize();
	return tmp.major_y == other.major_y && tmp.minor_y == other.minor_y;
}

bool Omesh2d::Coordinate::is_below(const Omesh2d::Coordinate &other) const
{
	return other.is_above(*this);
}

bool Omesh2d::Coordinate::is_left_of(const Omesh2d::Coordinate &other) const
{
	Omesh2d::Coordinate tmp = *this;
	tmp.minor_x += tmp.size;
	tmp.normalize();
	return tmp.major_x == other.major_x && tmp.minor_x == other.minor_x;
}

bool Omesh2d::Coordinate::is_right_of(const Omesh2d::Coordinate &other) const
{
	return other.is_left_of(*this);
}

int32_t Omesh2d::Coordinate::scale_from(const Omesh2d::Coordinate &other, int32_t v) const
{
	if (size == other.size) {
		return v;
	} else
	if (size == other.size >> 1) {
		return v < OMESH2D_HALF ? v << 1 : OMESH2D_FULL;
	} else
	if (size >> 1 == other.size) {
		return v >> 1;
	} else
		abort();
}

int32_t Omesh2d::Coordinate::scale_x_from(const Omesh2d::Coordinate &other, int32_t v) const
{
	assert(major_x == other.major_x);

	if (size == other.size) {
		return v;
	} else
	if (size == other.size >> 1) {
		assert(minor_x == other.minor_x || minor_x == other.minor_x + size);
		if (minor_x != other.minor_x)
			v -= OMESH2D_HALF;
		return v < 0 ? 0 : v < OMESH2D_HALF ? v << 1 : OMESH2D_FULL;
	} else
	if (size >> 1 == other.size) {
		if (minor_x != other.minor_x)
			return (v >> 1) + OMESH2D_HALF;
		return v >> 1;
	} else
		abort();
}

int32_t Omesh2d::Coordinate::scale_y_from(const Omesh2d::Coordinate &other, int32_t v) const
{
	assert(major_y == other.major_y);

	if (size == other.size) {
		return v;
	} else
	if (size == other.size >> 1) {
		assert(minor_y == other.minor_y || minor_y == other.minor_y + size);
		if (minor_y != other.minor_y)
			v -= OMESH2D_HALF;
		return v < 0 ? 0 : v < OMESH2D_HALF ? v << 1 : OMESH2D_FULL;
	} else
	if (size >> 1 == other.size) {
		if (minor_y != other.minor_y)
			return (v >> 1) + OMESH2D_HALF;
		return v >> 1;
	} else
		abort();
}

std::string Omesh2d::Coordinate::to_string() const
{
	char buffer[1024], *p = buffer;
	int mxs = OMESH2D_BITS - count_trailing_zeros(minor_x);
	int mys = OMESH2D_BITS - count_trailing_zeros(minor_y);
	int szs = OMESH2D_BITS - count_trailing_zeros(size);

	if (minor_x == 0) mxs = 0;
	if (minor_y == 0) mys = 0;

	snprintf(p, sizeof(buffer), "(%d %d | %d/%d %d/%d | 1/%d)", int(major_x), int(major_y),
			minor_x / (OMESH2D_FULL >> mxs), 1 << mxs,
			minor_y / (OMESH2D_FULL >> mys), 1 << mys,
			1 << szs);

	return std::string(buffer);
}

