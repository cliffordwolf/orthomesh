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

#include "gjk.h"

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <vector>
#include <set>

#include <Eigen/Core>
#include <Eigen/Geometry>
USING_PART_OF_NAMESPACE_EIGEN

#if 0
#  define printf_openscad printf
#else
#  define printf_openscad(...) do { } while (0)
#endif

uint32_t xor128(void)
{
	static uint32_t x = 123456789;
	static uint32_t y = 362436069;
	static uint32_t z = 521288629;
	static uint32_t w = 88675123;
	uint32_t t;
		   
	t = x ^ (x << 11);
	x = y; y = z; z = w;
	return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

template<typename T>
struct MyConvexHull : GJK<T>
{
	std::vector<T> points;

	void add(T x, T y, T z)
	{
		points.push_back(x);
		points.push_back(y);
		points.push_back(z);
	}

	bool empty_volume()
	{
		for (size_t i = 0; i < points.size(); i += 3)
		for (size_t j = i; j < points.size(); j += 3)
		for (size_t k = j; k < points.size(); k += 3)
		for (size_t l = 0; l < points.size(); l += 3) {
			Vector3f A, B, C, D;
			A << points[i+0], points[i+1], points[i+2];
			B << points[j+0], points[j+1], points[j+2];
			C << points[k+0], points[k+1], points[k+2];
			D << points[l+0], points[l+1], points[l+2];
			if (fabs((B-A).cross(C-A).dot(D-A)) > 1)
				return false;
		}
		return true;
	}

	bool share_points(MyConvexHull<T> &other)
	{
		for (size_t i = 0; i < points.size(); i += 3)
		for (size_t j = 0; j < other.points.size(); j += 3)
		{
			if (points[i+0] != other.points[j+0]) continue;
			if (points[i+1] != other.points[j+1]) continue;
			if (points[i+2] != other.points[j+2]) continue;
			return true;
		}
		return false;
	}

	virtual void gjk_support(T p[3])
	{
		GJK_Hull3D<T> worker(&points[0], points.size() / 3);
		worker.gjk_support(p);
	}
};

int my_sign(float v)
{
	if (v > +1e-6) return +1;
	if (v < -1e-6) return -1;
	return 0;
}

bool brute_force_in_tetrahedron(int T[12], float x, float y, float z)
{
	Vector3f A, B, C, D, P;

	A << T[0], T[ 1], T[ 2];
	B << T[3], T[ 4], T[ 5];
	C << T[6], T[ 7], T[ 8];
	D << T[9], T[10], T[11];
	P << x, y, z;

	Vector3f AB = B - A;
	Vector3f AC = C - A;
	Vector3f AD = D - A;
	Vector3f AP = P - A;

	Vector3f ABC = AB.cross(AC);
	Vector3f ACD = AC.cross(AD);
	Vector3f ADB = AD.cross(AB);

	if (my_sign(ABC.dot(AD)) * my_sign(ABC.dot(AP)) <= 0)
		return false;
	if (my_sign(ACD.dot(AB)) * my_sign(ACD.dot(AP)) <= 0)
		return false;
	if (my_sign(ADB.dot(AC)) * my_sign(ADB.dot(AP)) <= 0)
		return false;

	Vector3f BA = A - B;
	Vector3f BC = C - B;
	Vector3f BD = D - B;
	Vector3f BP = P - B;

	Vector3f BCD = BC.cross(BD);

	if ((BCD.dot(BA) > 0) != (BCD.dot(BP) > 0))
		return false;

	return true;
}

bool brute_force_in_convex_hull(MyConvexHull<float> &H, float x, float y, float z)
{
	for (size_t i = 0; i < H.points.size(); i += 3)
	for (size_t j = i+3; j < H.points.size(); j += 3)
	for (size_t k = j+3; k < H.points.size(); k += 3)
	for (size_t l = k+3; l < H.points.size(); l += 3)
	{
		int T[12];
		for (int n = 0; n < 3; n++) {
			T[n] = H.points[i+n];
			T[3+n] = H.points[j+n];
			T[6+n] = H.points[k+n];
			T[9+n] = H.points[l+n];
		}
		if (brute_force_in_tetrahedron(T, x, y, z))
			return true;
	}
	return false;
}

bool brute_force_collide_convex_hull(MyConvexHull<float> &A, MyConvexHull<float> &B)
{
	for (int x = 0; x <= 1000; x += 10)
	for (int y = 0; y <= 1000; y += 10)
	for (int z = 0; z <= 1000; z += 10)
		if (brute_force_in_convex_hull(A, x, y, z) && brute_force_in_convex_hull(B, x, y, z))
			return true;

	return false;
}

int main()
{
	int fail_counter = 0;

	for (int i = 0; i < 500; i++)
	{
		MyConvexHull<float> A, B;

		printf_openscad("--\n");
		printf_openscad("color([0.8, 0.3, 0.3]) hull() {\n");
		for (int j = 0; j < 4; j++) {
			A.points.push_back((xor128() % 11) * 100);
			A.points.push_back((xor128() % 11) * 100);
			A.points.push_back((xor128() % 11) * 100);
			printf_openscad("\ttranslate([ %6.0f, %6.0f, %6.0f ]) sphere(10);\n", A.points.at(A.points.size()-3),
					A.points.at(A.points.size()-2), A.points.at(A.points.size()-1));
		}
		printf_openscad("}\n");

		printf_openscad("color([0.3, 0.8, 0.8]) hull() {\n");
		for (int j = 0; j < 5; j++) {
			B.points.push_back((xor128() % 11) * 100);
			B.points.push_back((xor128() % 11) * 100);
			B.points.push_back((xor128() % 11) * 100);
			printf_openscad("\ttranslate([ %6.0f, %6.0f, %6.0f ]) sphere(10);\n", B.points.at(B.points.size()-3),
					B.points.at(B.points.size()-2), B.points.at(B.points.size()-1));
		}
		printf_openscad("}\n");

		if (i == 17 || i == 20 || i == 78 || i == 87 || i == 135 || i == 141 || i == 178 || i == 265 || i == 268 || i == 286 ||
				i == 315 || i == 318 || i == 330 || i == 388 || i == 406 || i == 430 || i == 482) {
			printf("%5d: disabled (after manual check)\n", i);
			continue;
		}

		if (A.share_points(B)) {
			printf("%5d: skipping (shared point)\n", i);
			continue;
		}

		if (A.empty_volume() || B.empty_volume()) {
			printf("%5d: skipping (no volume)\n", i);
			continue;
		}

		GJK_MinkowskiDifference<float> M(&A, &B);
		bool gjk_result = M.gjk_analyze(), brute_force_result = brute_force_collide_convex_hull(A, B);

		printf("%5d: %5s %5s\n", i, gjk_result ? "true" : "false", brute_force_result ? "true" : "false");

		if (gjk_result != brute_force_result) {
			printf("TEST FAILED!\n");
			fail_counter++;
		}
	}

	if (fail_counter > 0)
		printf("Got %d errors.\n", fail_counter);
	else
		printf("OK.\n");

	return fail_counter > 0;
}

