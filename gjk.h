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

// implementation of the Gilbert-Johnson-Keerthi (GJK) algorithm,
// based on the description given in the video lecture at http://mollyrocket.com/849
// (used the C# code from https://mollyrocket.com/forums/viewtopic.php?p=7551 as reference)

#ifndef GJK_H
#define GJK_H

template<typename T>
static inline void GJK_neg(T a[3])
{
	a[0] *= -1;
	a[1] *= -1;
	a[2] *= -1;
}

template<typename T>
static inline void GJK_sub(T a[3], T b[3])
{
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

template<typename T>
static inline void GJK_copy(T out[3], T a[3])
{
	out[0] = a[0];
	out[1] = a[1];
	out[2] = a[2];
}

template<typename T>
static inline T GJK_dot(T a[3], T b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

template<typename T>
static inline void GJK_cross(T a[3], T b[3])
{
	T buf[3] = {
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};

	a[0] = buf[0];
	a[1] = buf[1];
	a[2] = buf[2];
}

template<typename T>
static inline bool GJK_simplex(T points[4][3], int &points_n, T direction[3])
{
	if (points_n == 2)
	{
		T A[3], B[3];
		GJK_copy(A, points[1]);
		GJK_copy(B, points[0]);

		T AB[3];
		GJK_copy(AB, B);
		GJK_sub(AB, A);

		T AO[3];
		GJK_copy(AO, A);
		GJK_neg(AO);

		if (GJK_dot(AB, AO) > 0) {
			GJK_copy(direction, AB);
			GJK_cross(direction, AO);
			GJK_cross(direction, AB);
		} else {
			GJK_copy(direction, AO);
		}
	}
	else
	if (points_n == 3)
	{
		T A[3], B[3], C[3];
		GJK_copy(A, points[2]);
		GJK_copy(B, points[1]);
		GJK_copy(C, points[0]);

		T AO[3];
		GJK_copy(AO, A);
		GJK_neg(AO);

		T AB[3];
		GJK_copy(AB, B);
		GJK_sub(AB, A);

		T AC[3];
		GJK_copy(AC, C);
		GJK_sub(AC, A);

		T ABC[3];
		GJK_copy(ABC, AB);
		GJK_cross(ABC, AC);

		T ABC_AC[3];
		GJK_copy(ABC_AC, ABC);
		GJK_cross(ABC_AC, AC);

		T AB_ABC[3];
		GJK_copy(AB_ABC, AB);
		GJK_cross(AB_ABC, ABC);

		if (GJK_dot(ABC_AC, AO) > 0)
		{
			if (GJK_dot(AC, AO) > 0)
			{
				points_n = 2;
				GJK_copy(points[0], C);
				GJK_copy(points[1], A);
				GJK_copy(direction, AC);
				GJK_cross(direction, AO);
				GJK_cross(direction, AC);
			}
			else
			if (GJK_dot(AB, AO) > 0)
			{
				points_n = 2;
				GJK_copy(points[0], B);
				GJK_copy(points[1], A);
				GJK_copy(direction, AB);
				GJK_cross(direction, AO);
				GJK_cross(direction, AB);
			}
			else
			{
				points_n = 1;
				GJK_copy(points[0], A);
				GJK_copy(direction, AO);
			}
		}
		else
		if (GJK_dot(AB_ABC, AO) > 0)
		{
			if (GJK_dot(AB, AO) > 0)
			{
				points_n = 2;
				GJK_copy(points[0], B);
				GJK_copy(points[1], A);
				GJK_copy(direction, AB);
				GJK_cross(direction, AO);
				GJK_cross(direction, AB);
			}
			else
			{
				points_n = 1;
				GJK_copy(points[0], A);
				GJK_copy(direction, AO);
			}
		}
		else
		if (GJK_dot(ABC, AO) > 0)
		{
			GJK_copy(direction, ABC);
		}
		else
		{
			GJK_copy(points[0], B);
			GJK_copy(points[1], C);
			GJK_copy(points[2], A);
			GJK_copy(direction, ABC);
			GJK_neg(direction);
		}
	}
	else
	if (points_n == 4)
	{
		T A[3], B[3], C[3], D[3];
		GJK_copy(A, points[3]);
		GJK_copy(B, points[2]);
		GJK_copy(C, points[1]);
		GJK_copy(D, points[0]);

		T AO[3];
		GJK_copy(AO, A);
		GJK_neg(AO);

		T AB[3];
		GJK_copy(AB, B);
		GJK_sub(AB, A);

		T AC[3];
		GJK_copy(AC, C);
		GJK_sub(AC, A);

		T AD[3];
		GJK_copy(AD, D);
		GJK_sub(AD, A);

		T ABC[3];
		GJK_copy(ABC, AB);
		GJK_cross(ABC, AC);

		T ACD[3];
		GJK_copy(ACD, AC);
		GJK_cross(ACD, AD);

		T ADB[3];
		GJK_copy(ADB, AD);
		GJK_cross(ADB, AB);

		bool B_above_ACD = GJK_dot(ACD, AB) > 0;
		bool C_above_ADB = GJK_dot(ADB, AC) > 0;
		bool D_above_ABC = GJK_dot(ABC, AD) > 0;

		bool B_to_ACD_like_O = (GJK_dot(ACD, AO) > 0) == B_above_ACD;
		bool C_to_ADB_like_O = (GJK_dot(ADB, AO) > 0) == C_above_ADB;
		bool D_to_ABC_like_O = (GJK_dot(ABC, AO) > 0) == D_above_ABC;

		if (B_to_ACD_like_O && C_to_ADB_like_O && D_to_ABC_like_O)
			return true;

		if (!B_to_ACD_like_O)
		{
			points_n = 3;
			GJK_copy(points[2], A);
			GJK_copy(direction, ACD);
			if (B_above_ACD) GJK_neg(direction);
		}
		else
		if (!C_to_ADB_like_O)
		{
			points_n = 3;
			GJK_copy(points[2], A);
			GJK_copy(points[1], B);
			GJK_copy(direction, ADB);
			if (C_above_ADB) GJK_neg(direction);
		}
		else
		{
			points_n = 3;
			GJK_copy(points[2], A);
			GJK_copy(points[1], B);
			GJK_copy(points[0], C);
			GJK_copy(direction, ABC);
			if (D_above_ABC) GJK_neg(direction);
		}

		return GJK_simplex(points, points_n, direction);
	}

	return false;
}

template<typename T>
struct GJK
{
	virtual void gjk_support(T dir[3]) = 0;

	// this returns true if the convex body described by gjk_support() contains
	// the origin point. for intersection tests this is the Minkowski difference
	// of the two bodies.
	inline bool gjk_analyze(int max_iterations = 10, bool *loop_maxed = 0)
	{
		T direction[3] = { 1, 1, 1 }, tmp[3];
		T points[4][3];
		int points_n;

		gjk_support(direction);
		GJK_copy(points[0], direction);
		points_n = 1;

		GJK_neg(direction);

		if (loop_maxed)
			*loop_maxed = false;

		for (int i = 0; i < max_iterations; i++)
		{
			GJK_copy(tmp, direction);
			gjk_support(tmp);

			if (GJK_dot(tmp, direction) < 0)
				return false;

			GJK_copy(points[points_n], tmp);
			points_n++;

			if (GJK_simplex(points, points_n, direction))
				return true;
		}

		if (loop_maxed)
			*loop_maxed = true;

		return false;
	}
};

template<typename T>
struct GJK_MinkowskiDifference : GJK<T>
{
	GJK<T> *gjk_left, *gjk_right;
	GJK_MinkowskiDifference(GJK<T> *gjk_left, GJK<T> *gjk_right) : gjk_left(gjk_left), gjk_right(gjk_right) { }

	inline virtual void gjk_support(T dir[3])
	{
		T left[3];
		GJK_copy(left, dir);
		gjk_left->gjk_support(left);

		T right[3];
		GJK_copy(right, dir);
		GJK_neg(right);
		gjk_right->gjk_support(right);

		GJK_copy(dir, left);
		GJK_sub(dir, right);
	}
};

template<typename T>
struct GJK_Hull3D : GJK<T>
{
	T *gjk_points;
	int gjk_points_n;

	GJK_Hull3D(T *gjk_points, int gjk_points_n) :
			gjk_points(gjk_points), gjk_points_n(gjk_points_n) { }

	inline virtual void gjk_support(T dir[3])
	{
		T best_norm = 0;
		int best_idx = 0;

		for (int i = 0; i < 3*gjk_points_n; i += 3) {
			T this_norm = GJK_dot(dir, gjk_points+i);
			if (i == 0 || this_norm > best_norm)
				best_norm = this_norm, best_idx = i;
		}

		GJK_copy(dir, gjk_points+best_idx);
	}
};

template<typename T>
struct GJK_Hull2D : GJK<T>
{
	T *gjk_points;
	int gjk_points_n;

	GJK_Hull2D(T *gjk_points, int gjk_points_n) :
			gjk_points(gjk_points), gjk_points_n(gjk_points_n) { }

	inline virtual void gjk_support(T dir[3])
	{
		T best_norm = 0;
		int best_idx = 0;

		for (int i = 0; i < 2*gjk_points_n; i += 2) {
			T this_norm = dir[0]*gjk_points[i] + dir[1]*gjk_points[i+1];
			if (i == 0 || this_norm > best_norm)
				best_norm = this_norm, best_idx = i;
		}

		dir[0] = gjk_points[best_idx];
		dir[1] = gjk_points[best_idx+1];
		dir[2] = dir[2] > 0 ? +1 : -1;
	}
};

template<typename T>
struct GJK_Sphere : GJK<T>
{
	T x, y, z, r;

	GJK_Sphere(T x, T y, T z, T r) : x(x), y(y), z(z), r(r) { }

	inline virtual void gjk_support(T dir[3])
	{
		double scale = r / sqrt(GJK_dot(dir, dir));
		dir[0] = x + dir[0] * scale;
		dir[1] = y + dir[1] * scale;
		dir[2] = z + dir[2] * scale;
	}
};

#endif
