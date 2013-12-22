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
#include <list>

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

struct Lexer
{
	FILE *f;
	int linec;

	Lexer(FILE *f) : f(f), linec(1) { }

	std::map<std::string, std::vector<std::string>> macros;
	std::list<std::string> token_queue;
	std::string command_token;

	int fgetc_skip_comment()
	{
		int ch = fgetc(f);

		if (ch == '#') {
			do {
				ch = fgetc(f);
			} while (ch != '\r' && ch != '\n' && ch > 0);
		}

		return ch;
	}

	std::string next_token(bool except_eol = false)
	{
		int ch;
		std::string tok;

		if (!token_queue.empty()) {
			tok = token_queue.front();
			token_queue.pop_front();
			goto return_tok;
		}

		ch = fgetc_skip_comment();
		while (ch == ' ' || ch == '\t')
			ch = fgetc_skip_comment();

		if (ch == '\r' || ch == '\n' || ch == ',') {
			if (!except_eol) {
				fprintf(stderr, "Unexpected end of line or ',' in line %d!\n", linec);
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
			ch = fgetc_skip_comment();
		} while (ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n' && ch != ',');
		ungetc(ch, f);

		if (macros.count(tok) != 0) {
			token_queue.insert(token_queue.end(), macros.at(tok).begin(), macros.at(tok).end());
			return next_token(except_eol);
		}

	return_tok:
		if (command_token.empty())
			command_token = tok;
		return tok;
	}

	void next_line()
	{
		int count_r = 0, count_n = 0;
		int ch = fgetc_skip_comment();

		while (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n') {
			if (ch == '\r')
				count_r++;
			if (ch == '\n')
				count_n++;
			ch = fgetc_skip_comment();
		}

		if (ch == ',') {
			token_queue.insert(token_queue.end(), command_token);
			return;
		}

		if (count_r == 0 && count_n == 0) {
			fprintf(stderr, "Expected end of line in line %d, but got more input!\n", linec);
			exit(1);
		}

		ungetc(ch, f);
		linec += std::max(count_r, count_n);
		command_token = std::string();
	}

	int32_t parse_int(std::string token)
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

	double parse_float(std::string token)
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

	void other_error(std::string errmsg)
	{
		fprintf(stderr, "Error '%s' in line %d!\n", errmsg.c_str(), linec);
		exit(1);
	}

	void syntax_error(std::string token = std::string(), std::string expected = std::string())
	{
		if (token.empty())
			fprintf(stderr, "Syntax error in line %d!\n", linec);
		else if (expected.empty())
			fprintf(stderr, "Syntax error near token '%s' in line %d!\n", token.c_str(), linec);
		else
			fprintf(stderr, "Syntax error near token '%s' in line %d, expected %s!\n", token.c_str(), linec, expected.c_str());
		exit(1);
	}
};


/*************************************************************************
 *                            orthomesh 2d                               *
 *************************************************************************/

struct geometry_2d_data_t
{
	std::string tag, color;
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
		GJK_Hull2D<double> gjk_box(gjk_points, 8);

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

	virtual std::string getcolor(int32_t major_x, int32_t major_y, std::vector<int32_t> &minor_xy)
	{
		std::string default_color = "gray";

		std::vector<double> gjk_points;
		for (size_t i = 0; i < minor_xy.size(); i += 2) {
			gjk_points.push_back(major_x + minor_xy[i+0] * OMESH2D_STEP);
			gjk_points.push_back(major_y + minor_xy[i+1] * OMESH2D_STEP);
		}

		GJK_Hull2D<double> gjk_hull(&gjk_points[0], gjk_points.size());

		for (auto &g : geometries)
		{
			if (g.color.empty())
				continue;

			if (g.default_geometry)
				default_color = g.color;

			for (auto gjk_ptr : g.gjk_supports) {
				GJK_MinkowskiDifference<double> gjk(&gjk_hull, gjk_ptr);
				if (gjk.gjk_analyze())
					return g.color;
			}
		}

		return default_color;
	}
};

static void orthomesh_2d(Lexer &lex, FILE *f_out, orthomes_options_t &options)
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
		lex.next_line();
		tok = lex.next_token();

		if (tok == "define") {
			std::string macro_name = lex.next_token();
			std::vector<std::string> macro_body;
			while (1) {
				tok = lex.next_token(true);
				if (tok.empty())
					break;
				macro_body.push_back(tok);
			}
			lex.macros[macro_name] = macro_body;
			continue;
		}

		if (tok == "verbose") {
			if (set_verbose)
				lex.other_error("duplicate verbose");
			verbose = lex.parse_int(lex.next_token());
			set_verbose = true;
			continue;
		}

		if (tok == "split_slice") {
			if (set_split_slice)
				lex.other_error("duplicate split_slice");
			split_slice = lex.parse_int(lex.next_token());
			set_split_slice = true;
			continue;
		}

		if (tok == "grid_w") {
			if (set_grid_w)
				lex.other_error("duplicate grid_w");
			grid_w = lex.parse_int(lex.next_token());
			set_grid_w = true;
			continue;
		}

		if (tok == "grid_h") {
			if (set_grid_h)
				lex.other_error("duplicate grid_h");
			grid_h = lex.parse_int(lex.next_token());
			set_grid_h = true;
			continue;
		}

		if (tok == "offset_x") {
			offset_x = lex.parse_float(lex.next_token());
			continue;
		}

		if (tok == "offset_y") {
			offset_y = lex.parse_float(lex.next_token());
			continue;
		}

		if (tok == "scale") {
			scale = lex.parse_float(lex.next_token());
			continue;
		}

		if (tok == "grid_scale") {
			grid_scale = lex.parse_float(lex.next_token());
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
				lex.next_line();
				tok = lex.next_token();

				if (tok == "default" || tok == "points" || tok == "lines" || tok == "triangles") {
					if (!mode.empty())
						lex.other_error("duplicate mode");
					mode = tok;
					continue;
				}

				if (tok == "point") {
					double x = lex.parse_float(lex.next_token());
					double y = lex.parse_float(lex.next_token());
					xy_data.push_back((x * scale + offset_x) / grid_scale);
					xy_data.push_back((y * scale + offset_y) / grid_scale);
					continue;
				}

				if (tok == "radius") {
					if (set_radius)
						lex.other_error("duplicate radius");
					double r = lex.parse_float(lex.next_token());
					radius = r * scale / grid_scale;
					set_radius = true;
					continue;
				}

				if (tok == "resolution_x") {
					double r = lex.parse_float(lex.next_token());
					g.resolution_x = std::min(g.resolution_x, r * scale / grid_scale);
					continue;
				}

				if (tok == "resolution_y") {
					double r = lex.parse_float(lex.next_token());
					g.resolution_y = std::min(g.resolution_y, r * scale / grid_scale);
					continue;
				}

				if (tok == "color") {
					if (!g.color.empty())
						lex.other_error("duplicate color");
					g.color = lex.next_token();
					continue;
				}

				if (tok == "tag") {
					if (!g.tag.empty())
						lex.other_error("duplicate tag");
					g.tag = lex.next_token();
					continue;
				}

				if (tok == "endgeometry")
					break;

				lex.syntax_error("'default', 'points', 'lines', 'triangles', 'point', "
						"'radius', 'resolution_x', 'resolution_y', 'color', 'tag', or 'endgeometry'");
			}

			if (mode == "default")
			{
				if (!xy_data.empty())
					lex.other_error("default gemotries may not contain any points");
				if (set_radius)
					lex.other_error("default gemotries may not contain a radius");
				g.default_geometry = true;
			}
			else
			{
				if (xy_data.empty())
					lex.other_error("need at least one point");

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

					g.gjk_supports.push_back(new GJK_Hull2D<double>(&g.gjk_data.back()->at(0), g.gjk_data.back()->size()));

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

		lex.syntax_error("'grid_w', 'grid_h', 'offset_x', 'offset_y', 'scale', 'grid_scale', 'geometry', or 'end'");
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

	Lexer lex(f_in);

	std::string initial_token = lex.next_token(true);
	if (initial_token.empty()) {
		lex.next_line();
		initial_token = lex.next_token();
	}

	if (initial_token == "orthomesh_2d")
		orthomesh_2d(lex, f_out, options);
	else
		lex.syntax_error(initial_token, "'orthomesh_2d'");

	return 0;
}

