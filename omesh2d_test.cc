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
#include <math.h>

struct DemoMesh : Omesh2d
{
	DemoMesh() : Omesh2d(10, 5)
	{
		config_verbose = 10;
	}

	bool bbox_intersect_rect(double bbox[4], double x1, double y1, double x2, double y2, double eps = 0)
	{
		if (bbox[OMESH2D_RIGHT] <= x1 + eps || x2 - eps <= bbox[OMESH2D_LEFT])
			return false;
		if (bbox[OMESH2D_BOTTOM] <= y1 + eps || y2 - eps <= bbox[OMESH2D_TOP])
			return false;
		return true;
	}

	bool bbox_intersect_poly2(double bbox[4], double a, double b, double c)
	{
		for (double x = bbox[OMESH2D_LEFT]; x <= bbox[OMESH2D_RIGHT]; x += (bbox[OMESH2D_RIGHT] - bbox[OMESH2D_LEFT]) / 16.0) {
			double y = a*x*x + b*x + c;
			if (bbox[OMESH2D_TOP] <= y && y <= bbox[OMESH2D_BOTTOM])
				return true;
		}
		return false;
	}

	virtual void refine(int32_t major_x, int32_t major_y, int32_t minor_x, int32_t minor_y, int32_t size_x, int32_t size_y, bool &refine_x, bool &refine_y)
	{
		double max_x = 1.0, max_y = 1.0;
		double bbox[4];

		// the real bounding box of the grid is x in [-5 .. +5] and y in [-2 .. +3].
		double grid2bbox_scale =  1.0;
		double grid2bbox_off_x = -5.0;
		double grid2bbox_off_y = -2.0;

		bbox[OMESH2D_LEFT]   = grid2bbox_scale * (major_x + minor_x * exp2(-30)) + grid2bbox_off_x;
		bbox[OMESH2D_RIGHT]  = bbox[OMESH2D_LEFT] + grid2bbox_scale * size_x * exp2(-30);

		bbox[OMESH2D_TOP]    = grid2bbox_scale * (major_y + minor_y * exp2(-30)) + grid2bbox_off_y;
		bbox[OMESH2D_BOTTOM] = bbox[OMESH2D_TOP] + grid2bbox_scale * size_y * exp2(-30);

		// left diffusion well
		if (0 <= bbox[OMESH2D_TOP] && bbox_intersect_poly2(bbox, -2, -12, -16))
			max_x = std::min(max_x, 0.1), max_y = std::min(max_y, 0.1);

		// right diffusion well
		if (0 <= bbox[OMESH2D_TOP] && bbox_intersect_poly2(bbox, -2, 12, -16))
			max_x = std::min(max_x, 0.1), max_y = std::min(max_y, 0.1);

		// disturbance
		if (bbox_intersect_poly2(bbox, 0, 2, 0.5))
			max_x = std::min(max_x, 0.2), max_y = std::min(max_y, 0.2);

		// channel
		if (bbox_intersect_rect(bbox, -2, 0, 2, 0.05))
			max_x = std::min(max_x, 0.3), max_y = std::min(max_y, 0.01);

		// extra refinement test (right side of screen)
		if (bbox_intersect_rect(bbox, 2.7, -1.3, 3.0, -1.0, 0.01)) {
			max_x = std::min(max_x, std::max(0.01,  2.95 - bbox[OMESH2D_RIGHT]));
			max_y = std::min(max_y, std::max(0.01, -1.05 - bbox[OMESH2D_BOTTOM]));
		}
		if (bbox_intersect_rect(bbox, 3.0, -1.3, 3.3, -1.0, 0.01)) {
			max_x = std::min(max_x, std::max(0.01, -3.05 + bbox[OMESH2D_LEFT]));
			max_y = std::min(max_y, std::max(0.01, -1.05 - bbox[OMESH2D_BOTTOM]));
		}
		if (bbox_intersect_rect(bbox, 2.5, -1.0, 3.0, -0.9, 0.01)) {
			max_x = std::min(max_x, std::max(0.02,  2.95 - bbox[OMESH2D_RIGHT]));
			max_y = std::min(max_y, std::max(0.02,  0.95 + bbox[OMESH2D_BOTTOM]));
		}
		if (bbox_intersect_rect(bbox, 3.0, -1.0, 3.5, -0.9, 0.01)) {
			max_x = std::min(max_x, std::max(0.05, -3.05 + bbox[OMESH2D_LEFT]));
			max_y = std::min(max_y, std::max(0.05,  0.95 + bbox[OMESH2D_BOTTOM]));
		}

		// extra refinement test (left side of screen)
		if (bbox_intersect_rect(bbox, -4.5, -1.5, -4.0, -1.0, 0.01)) {
			max_x = std::min(max_x, 0.05);
			max_y = std::min(max_y, 0.5);
		}
		if (bbox_intersect_rect(bbox, -4.0, -1.5, -3.5, -1.0, 0.01)) {
			max_x = std::min(max_x, 0.5);
			max_y = std::min(max_y, 0.05);
		}

		// extra refinement test (bottom of screen)
		if (bbox_intersect_rect(bbox, -2.0, 2.5, 0.0, 3.0, 0.01)) {
			max_x = std::min(max_x, 0.5);
			max_y = std::min(max_y, std::max(0.01, (2.0 + bbox[OMESH2D_LEFT]) * 0.2));
		}

		refine_x = bbox[OMESH2D_RIGHT]  - bbox[OMESH2D_LEFT] > max_x;
		refine_y = bbox[OMESH2D_BOTTOM] - bbox[OMESH2D_TOP]  > max_y;
	}
};

int main()
{
	FILE *f;

	DemoMesh mesher;
	mesher.run();

	Omesh2d celltest(1, 1);
	for (int i = 0; i < 15; i++) {
		Omesh2d::CellType ctype;
		ctype.split_at[0] = (i & 1) != 0;
		ctype.split_at[1] = (i & 2) != 0;
		ctype.split_at[2] = (i & 4) != 0;
		ctype.split_at[3] = (i & 8) != 0;
		celltest.create_geometry(ctype);
	}
	for (int i = 1; i < 256; i++)
	for (int j = 1; j < 256; j++) {
		Omesh2d::CellType ctype;
		for (int k = 1; k < 8; k++) {
			if ((i & (1 << k)) != 0) ctype.slice_at[OMESH2D_LEFT].push_back((OMESH2D_FULL >> 3) * k);
			if ((j & (1 << k)) != 0) ctype.slice_at[OMESH2D_RIGHT].push_back((OMESH2D_FULL >> 3) * k);
		}
		celltest.fixup_slice_vector(ctype.slice_at[OMESH2D_LEFT]);
		celltest.fixup_slice_vector(ctype.slice_at[OMESH2D_RIGHT]);
		if (celltest.check_slice_compatibility(ctype.slice_at[OMESH2D_LEFT], ctype.slice_at[OMESH2D_RIGHT])) {
			celltest.create_geometry(ctype);
			ctype.split_at[OMESH2D_TOP] = true;
			celltest.create_geometry(ctype);
		}
	}

	f = fopen("omesh2d_test.svg", "wt");
	mesher.svg_write_grid(f);
	fclose(f);

	f = fopen("omesh2d_test.html", "wt");

	fprintf(f, "<h3>Generic Cell Type Tests</h3>\n");
	fprintf(f, "<table width='1200' border cellpadding='5'>\n");
	fprintf(f, "<tr><th width='80%%'>Cell Description</th><th>Delaunay Map</th><th>Delaunay + Voronoi Map</th></tr>\n");
	for (auto &it : celltest.geometries)
		if (!it.second.segments.empty())
			celltest.html_write_celltype(f, it.first);
	fprintf(f, "</table>\n");

	fprintf(f, "<h3>Cell Types From The Example</h3>\n");
	fprintf(f, "<table width='1200' border cellpadding='5'>\n");
	fprintf(f, "<tr><th width='80%%'>Cell Description</th><th>Delaunay Map</th><th>Delaunay + Voronoi Map</th></tr>\n");
	for (auto &it : mesher.geometries)
		if (!it.second.segments.empty())
			mesher.html_write_celltype(f, it.first);
	fprintf(f, "</table>\n");

	fclose(f);

	return 0;
}

