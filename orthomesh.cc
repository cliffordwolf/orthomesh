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
#include "gjk.h"

#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <string>

struct orthomes_options_t
{
	std::string format;
	double scale;
	bool plain;
	int verbose;

	orthomes_options_t() {
		scale = 100;
		plain = false;
		verbose = 0;
	}
};

/*************************************************************************
 *                              lexer API                                *
 *************************************************************************/

static int fgetc_skip_comment(FILE *f)
{
	int ch = fgetc(f);

	if (ch == '#') {
		do {
			ch = fgetc(f);
		} while (ch != '\r' && ch != '\n' && ch > 0);
	}

	return ch;
}

static std::string next_token(FILE *f, int linec, bool except_eol = false)
{
	std::string tok;

	int ch = fgetc_skip_comment(f);
	while (ch == ' ' || ch == '\t')
		ch = fgetc_skip_comment(f);

	if (ch == '\r' || ch == '\n') {
		if (!except_eol) {
			fprintf(stderr, "Unexpected end of line in line %d!\n", linec);
			exit(1);
		}
		ungetc(ch, f);
		return tok;
	}

	do {
		if (ch < 0) {
			fprintf(stderr, "Unexpected end of file in line %d!\n", linec);
			exit(1);
		}

		tok += ch;
		ch = fgetc_skip_comment(f);
	} while (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n');

	ungetc(ch, f);
	return tok;
}

static void next_line(FILE *f, int &linec)
{
	int count_r = 0, count_n = 0;
	int ch = fgetc_skip_comment(f);

	while (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n') {
		if (ch == '\r')
			count_r++;
		if (ch == '\n')
			count_n++;
		ch = fgetc_skip_comment(f);
	}

	if (count_r == 0 && count_n == 0) {
		fprintf(stderr, "Expected end of line in line %d, but got more input!\n", linec);
		exit(1);
	}

	ungetc(ch, f);
	linec += std::max(count_r, count_n);
}

static int32_t parse_int(int linec, std::string token)
{
	const char *p = token.c_str();
	char *endptr;
	long long v = strtoll(p, &endptr, 0);

	if (*p == 0 || *endptr != 0) {
		fprintf(stderr, "Expected integer in line %d, but got '%s' instead!\n", linec, p);
		exit(1);
	}

	return v;
}

static double parse_float(int linec, std::string token)
{
	const char *p = token.c_str();
	char *endptr;
	double v = strtod(p, &endptr);

	if (*p == 0 || *endptr != 0) {
		fprintf(stderr, "Expected real number in line %d, but got '%s' instead!\n", linec, p);
		exit(1);
	}

	return v;
}

static void other_error(int linec, std::string errmsg)
{
	fprintf(stderr, "Error '%s' in line %d!\n", errmsg.c_str(), linec);
	exit(1);
}

static void syntax_error(int linec, std::string token = std::string(), std::string expected = std::string())
{
	if (token.empty())
		fprintf(stderr, "Syntax error in line %d!\n", linec);
	else if (expected.empty())
		fprintf(stderr, "Syntax error near token '%s' in line %d!\n", token.c_str(), linec);
	else
		fprintf(stderr, "Syntax error near token '%s' in line %d, expected %s!\n", token.c_str(), linec, expected.c_str());
	exit(1);
}


/*************************************************************************
 *                            orthomesh 2d                               *
 *************************************************************************/

struct geometry_2d_data_t
{
	std::string tag;
	double resolution_x, resolution_y;
	bool default_geometry;

	std::vector<std::vector<double>*> gjk_data;
	std::vector<GJK<double>*> gjk_tmp_supports;
	std::vector<GJK<double>*> gjk_supports;

	geometry_2d_data_t()
	{
		resolution_x = 1.0;
		resolution_y = 1.0;
		default_geometry = false;
	}

	void free_gjk_stuff()
	{
		for (auto ptr : gjk_data)
			delete ptr;
		for (auto ptr : gjk_tmp_supports)
			delete ptr;
		for (auto ptr : gjk_supports)
			delete ptr;
	}
};

struct Omesh2d_GJK : Omesh2d
{
	std::vector<geometry_2d_data_t> &geometries;

	Omesh2d_GJK(int32_t major_width, int32_t major_height, std::vector<geometry_2d_data_t> &geometries) :
			Omesh2d(major_width, major_height), geometries(geometries) { }

	virtual void refine(int32_t major_x, int32_t major_y, int32_t minor_x, int32_t minor_y, int32_t size_x, int32_t size_y, bool &refine_x, bool &refine_y)
	{
		double max_x = 1.0, max_y = 1.0;
		
		double x1 = major_x + minor_x * OMESH2D_STEP;
		double y1 = major_y + minor_y * OMESH2D_STEP;
		double x2 = x1 + size_x * OMESH2D_STEP;
		double y2 = y1 + size_y * OMESH2D_STEP;

		double gjk_points[8] = { x1, y1, x2, y1, x2, y2, x1, y2 };
		GJK_Hull2D<double> gjk_box(gjk_points, 4);

		for (auto &g : geometries)
		{
			for (auto gjk_ptr : g.gjk_supports) {
				GJK_MinkowskiDifference<double> gjk(&gjk_box, gjk_ptr);
				if (gjk.gjk_analyze())
					goto gjk_matched;
			}

			if (!g.default_geometry)
				continue;

		gjk_matched:
			max_x = std::min(max_x, g.resolution_x);
			max_y = std::min(max_y, g.resolution_y);
		}

		refine_x = max_x < size_x * OMESH2D_STEP;
		refine_y = max_y < size_y * OMESH2D_STEP;
	}
};

static void orthomesh_2d(FILE *f, int &linec, FILE *f_out, orthomes_options_t &options)
{
	std::string tok;

	bool set_verbose = false;
	bool set_split_slice = false;

	int verbose = 0;
	int split_slice = 0;

	bool set_grid_w = false;
	bool set_grid_h = false;

	int grid_w = 1, grid_h = 1;
	double offset_x = 0, offset_y = 0;
	double scale = 1, grid_scale = 1;
	std::vector<geometry_2d_data_t> geometries;

	while (1)
	{
		next_line(f, linec);
		tok = next_token(f, linec);

		if (tok == "verbose") {
			if (set_verbose)
				other_error(linec, "duplicate verbose");
			verbose = parse_int(linec, next_token(f, linec));
			set_verbose = true;
			continue;
		}

		if (tok == "split_slice") {
			if (set_split_slice)
				other_error(linec, "duplicate split_slice");
			split_slice = parse_int(linec, next_token(f, linec));
			set_split_slice = true;
			continue;
		}

		if (tok == "grid_w") {
			if (set_grid_w)
				other_error(linec, "duplicate grid_w");
			grid_w = parse_int(linec, next_token(f, linec));
			set_grid_w = true;
			continue;
		}

		if (tok == "grid_h") {
			if (set_grid_h)
				other_error(linec, "duplicate grid_h");
			grid_h = parse_int(linec, next_token(f, linec));
			set_grid_h = true;
			continue;
		}

		if (tok == "offset_x") {
			offset_x = parse_float(linec, next_token(f, linec));
			continue;
		}

		if (tok == "offset_y") {
			offset_y = parse_float(linec, next_token(f, linec));
			continue;
		}

		if (tok == "scale") {
			scale = parse_float(linec, next_token(f, linec));
			continue;
		}

		if (tok == "grid_scale") {
			grid_scale = parse_float(linec, next_token(f, linec));
			continue;
		}

		if (tok == "geometry")
		{
			std::string mode;
			bool set_radius = false;
			double radius = 0;

			geometry_2d_data_t g;
			std::vector<double> xy_data;

			while (1)
			{
				next_line(f, linec);
				tok = next_token(f, linec);

				if (tok == "default" || tok == "points" || tok == "lines" || tok == "triangles") {
					if (!mode.empty())
						other_error(linec, "duplicate mode");
					mode = tok;
					continue;
				}

				if (tok == "point") {
					double x = parse_float(linec, next_token(f, linec));
					double y = parse_float(linec, next_token(f, linec));
					xy_data.push_back((x * scale + offset_x) / grid_scale);
					xy_data.push_back((y * scale + offset_y) / grid_scale);
					continue;
				}

				if (tok == "radius") {
					if (set_radius)
						other_error(linec, "duplicate radius");
					double r = parse_float(linec, next_token(f, linec));
					radius = r * scale / grid_scale;
					set_radius = true;
					continue;
				}

				if (tok == "resolution_x") {
					double r = parse_float(linec, next_token(f, linec));
					g.resolution_x = std::min(g.resolution_x, r * scale / grid_scale);
					continue;
				}

				if (tok == "resolution_y") {
					double r = parse_float(linec, next_token(f, linec));
					g.resolution_y = std::min(g.resolution_y, r * scale / grid_scale);
					continue;
				}

				if (tok == "tag") {
					if (!g.tag.empty())
						other_error(linec, "duplicate tag");
					g.tag = next_token(f, linec);
					continue;
				}

				if (tok == "endgeometry")
					break;

				syntax_error(linec, tok, "'default', 'points', 'lines', 'triangles', 'point', "
						"'radius', 'resolution_x', 'resolution_y', 'tag', or 'endgeometry'");
			}

			if (mode == "default")
			{
				if (!xy_data.empty())
					other_error(linec, "default gemotries may not contain any points");
				if (set_radius)
					other_error(linec, "default gemotries may not contain a radius");
				g.default_geometry = true;
			}
			else
			{
				if (xy_data.empty())
					other_error(linec, "need at least one point");

				size_t cursor = 0;

				if (mode == "lines")
					cursor += 2;

				if (mode == "triangles")
					cursor += 4;

				while (cursor < xy_data.size())
				{
					g.gjk_data.push_back(new std::vector<double>());

					if (mode == "points") {
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
					} else
					if (mode == "lines") {
						g.gjk_data.back()->push_back(xy_data.at(cursor-2));
						g.gjk_data.back()->push_back(xy_data.at(cursor-1));
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
					} else
					if (mode == "triangles") {
						g.gjk_data.back()->push_back(xy_data.at(cursor-4));
						g.gjk_data.back()->push_back(xy_data.at(cursor-3));
						g.gjk_data.back()->push_back(xy_data.at(cursor-2));
						g.gjk_data.back()->push_back(xy_data.at(cursor-1));
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
						g.gjk_data.back()->push_back(xy_data.at(cursor++));
					} else
						g.gjk_data.back()->swap(xy_data);

					g.gjk_supports.push_back(new GJK_Hull2D<double>(&g.gjk_data.back()->at(0), g.gjk_data.back()->size()/2));

					if (set_radius)
					{
						GJK<double> *tmp1 = g.gjk_supports.back();
						g.gjk_tmp_supports.push_back(tmp1);

						GJK<double> *tmp2 = new GJK_Sphere<double>(0, 0, 0, radius);
						g.gjk_tmp_supports.push_back(tmp2);

						g.gjk_supports.back() = new GJK_MinkowskiDifference<double>(tmp1, tmp2);
					}
				}
			}

			geometries.push_back(g);
			continue;
		}

		if (tok == "end")
			break;

		syntax_error(linec, tok, "'grid_w', 'grid_h', 'offset_x', 'offset_y', 'scale', 'grid_scale', 'geometry', or 'end'");
	}

	Omesh2d_GJK mesher(grid_w, grid_h, geometries);
	mesher.config_verbose = std::max(options.verbose, verbose);
	mesher.config_split_slice_aggressiveness = split_slice;
	mesher.run();

	if (options.format.empty() || options.format == "svg") {
		mesher.svg_write_grid(f_out, options.scale, !options.plain);
	} else {
		fprintf(stderr, "Unsupported format '%s' for 2d mesher!\n", options.format.c_str());
	}

	for (auto &g : geometries)
		g.free_gjk_stuff();
}

/*************************************************************************
 *                                 main                                  *
 *************************************************************************/

int main(int argc, char **argv)
{
	orthomes_options_t options;
	FILE *f_out = stdout;
	FILE *f_in = stdin;
	int linec = 1;

	int opt;
	while ((opt = getopt(argc, argv, "vo:f:ps:i:")) != -1)
	{
		switch (opt)
		{
		case 'v':
			options.verbose++;
			break;

		case 'o':
			f_out = fopen(optarg, "wt");
			if (f_out != NULL)
				break;
			fprintf(stderr, "Can't open output file: %s\n", strerror(errno));
			return 1;

		case 'f':
			options.format = optarg;
			break;

		case 'p':
			options.plain = true;
			break;

		case 's':
			options.scale = atof(optarg);
			if (options.scale > 0)
				break;
			goto print_help;

		case 'i':
			f_in = fopen(optarg, "rt");
			if (f_in != NULL)
				break;
			fprintf(stderr, "Can't open output file: %s\n", strerror(errno));
			return 1;

		default:
			goto print_help;
		}
	}

	if (optind < argc) {
print_help:
		fprintf(stderr, "Usage: %s [-v..] [-o <outfile>] [-f {svg|omm}] [-p] [-s <scale>] [-i <infile>]\n", argv[0]);
		return 1;
	}

	std::string initial_token = next_token(f_in, linec, true);
	if (initial_token.empty()) {
		next_line(f_in, linec);
		initial_token = next_token(f_in, linec);
	}

	if (initial_token == "orthomesh_2d")
		orthomesh_2d(f_in, linec, f_out, options);
	else
		syntax_error(linec, initial_token, "'orthomesh_2d'");

	return 0;
}

