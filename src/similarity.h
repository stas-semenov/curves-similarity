/*
*  Copyright © 2015-2016 Stanislav Semenov. All rights reserved.
*  Contacts: <stas.semenov@gmail.com>
*  https://github.com/stas-semenov
*
*  This C++ library implements algorithm for computing similarity between curves.
*
*  Source code is provided under the GNU General Public License v3.
*  To access the source code please visit the public code repository.
*/

#pragma once

#define NOMINMAX

#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>

#include "fcm.h"
#include "naive.h"

#define EPSILON (1E-15)
#define INF (std::numeric_limits<double>::infinity())
#define CONTEXT_ANGLE_BIN_SIZE 6
#define CONTEXT_DIST_BIN_SIZE 8
#define CONTEXT_WINDOW_SIZE 32
#define BITMAP_SIZE 144

namespace Similarity
{
	using std::min;
	using std::max;

	typedef std::pair<double, int> ValueKey;
	typedef std::vector<ValueKey> ValueKeyVec;
	typedef std::vector<std::string> Strings;

	static double max_min(double left, double right, double value)
	{
		return max(left, min(right, value));
	}

	static std::string int_str(int value, const char * pref = "", const char * suff = "")
	{
		std::stringstream parser;
		parser << pref << value << suff;
		std::string res;
		parser >> res;

		return res;
	}

	static double calc_dist(double left, double right, int dist_type = 2, int norm_type = 1, double p = 1.)
	{
		double norm = 0;
		switch (norm_type)
		{
		case -1:
			norm = 1;
			break;
		case 1:
			norm = max(max(abs(left), abs(right)), EPSILON);
			break;
		default:
			norm = max(sqrt(pow(left, 2) + pow(right, 2)), EPSILON);
			break;
		}

		switch (dist_type)
		{
		case 1:
			return log(1. + abs(left - right) / norm) / log(2.);
		case 2:
			return (exp(abs(left - right) / norm) - 1.) / (exp(1) - 1.);
		default:
			return pow(abs(left - right), p) / norm;
		}
	}

	class SmartCurve
	{
	public:
		struct Coord;

		struct PolarCoord
		{
			double R, Theta;

			PolarCoord()
			{
				R = 0;
				Theta = 0;
			}

			PolarCoord(double r, double theta)
			{
				R = r;
				Theta = theta;
			}

			Coord to_cartesian() const
			{
				return Coord(exp(R) * cos(Theta), exp(R) * sin(Theta));
			}
		};

		struct Coord
		{
			double X, Y, Value;

			Coord()
			{
				X = 0;
				Y = 0;
				Value = 1;
			}

			Coord(double x, double y, double val = 1)
			{
				X = x;
				Y = y;
				Value = val;
			}

			Coord(const Coord& left, const Coord& right)
			{
				X = right.X - left.X;
				Y = right.Y - left.Y;
				Value = 1;
			}

			PolarCoord to_polar() const
			{
				return PolarCoord(log(sqrt(pow(X, 2) + pow(Y, 2))), atan2(Y, X));
			}

			int get_quadrant()
			{
				if (X >= 0 && Y >= 0)
					return 1;
				else if (X < 0 && Y >= 0)
					return 2;
				else if (X < 0 && Y < 0)
					return 3;
				else 
					return 4;
			}

			double distance_to_point(const Coord& P) const
			{
				return sqrt(pow((X - P.X), 2) + pow((Y - P.Y), 2));
			}

			double distance_to_line(const Coord& a, const Coord& b) const
			{
				double abX = b.X - a.X;
				double abY = b.Y - a.Y;
				double c = a.X * b.Y - b.X * a.Y;
				double ab = sqrt(pow(abX, 2) + pow(abY, 2));

				if (ab > 0)
					return abs((abX * Y - abY * X + c) / ab);
				else
					return 0;
			}
		};

		struct Point;
		struct Segment;
		struct Curve;

		typedef std::vector<double> Vector;
		typedef std::vector<Vector> Matrix;

		typedef std::vector<int> BitmapRow;
		typedef std::vector<BitmapRow> Bitmap;

		typedef std::map<int, int> SparseBitmapRow;
		typedef std::map<int, SparseBitmapRow> SparseBitmap;

		typedef Coord * CoordPtr;
		typedef std::vector<Coord> Coords;
		typedef std::vector<CoordPtr> CoordsPtr;

		typedef Point * PointPtr;
		typedef std::vector<Point> Points;
		typedef std::vector<PointPtr> PointsPtr;
		typedef std::vector<Segment> Segments;

		struct Raster
		{
			Bitmap bm;
			SparseBitmap sbm;
			bool sparse_flag;
			int size;
			Coords coords;

			Raster()
			{
				sparse_flag = false;
				size = BITMAP_SIZE;
			}

			void init(bool sparse = false)
			{
				bm.clear();
				sbm.clear();
				coords.clear();
				sparse_flag = sparse;
				if (!sparse_flag)
					bm.resize(size + 1, BitmapRow(size + 1, 0));
			}

			void fill_coords()
			{
				coords.clear();
				if (sparse_flag) {
					for (SparseBitmap::iterator itc = sbm.begin(); itc != sbm.end(); itc++)
						for (SparseBitmapRow::iterator itr = itc->second.begin(); itr != itc->second.end(); itr++)
							coords.push_back(Coord(itc->first, itr->first));
				} else {
					for (int i = 0; i < size; i++)
						for (int j = 0; j < size; j++)
							if (bm[i][j])
								coords.push_back(Coord(i, j));
				}
			}

			bool is_xy(int x, int y)
			{
				if ((sparse_flag && sbm.find(x) != sbm.end() && sbm[x].find(y) != sbm[x].end()) ||
					(!sparse_flag && x < size && y < size))
					return true;

				return false;
			}

			int get(int x, int y)
			{
				if (sparse_flag)
					return sbm[x][y];
				else
					return bm[x][y];
			}

			bool set(int x, int y, int val = 1)
			{
				if (sparse_flag)
					sbm[x][y] = val;
				else
					bm[x][y] = val;

				return true;
			}
		};

		static double calc_angle_abc(const Coord& a, const Coord& b, const Coord& c)
		{
			double abX = b.X - a.X;
			double abY = b.Y - a.Y;
			double cbX = b.X - c.X;
			double cbY = b.Y - c.Y;

			double dot = abX * cbX + abY * cbY;
			double cross = abX * cbY - abY * cbX;

			return atan2(cross, dot);
		}

		static double calc_angle_xy(const Coord& a, const Coord& b)
		{
			double dX = b.X - a.X;
			double dY = b.Y - a.Y;

			return atan2(dY, dX);
		}

		static double distance_to_features(const Point& left, const Point& right, int nfunc)
		{
			return left.distance_to_features(right, nfunc);
		}

		// Dynamic time warping algorithm
		static double dtw_simple_algorithm(const PointsPtr& left, const PointsPtr& right, int wnd = 0, bool ret_median = true, int nfunc = 0, double (*d_f)(const Point&, const Point&, int) = distance_to_features)
		{
			int N = left.size();
			int M = right.size();
			double ratio = (double)max(N, M) / max((double)min(N, M), EPSILON);

			Vector diff;
			if (wnd == 0)
				wnd = N / 16;
			int window = (int)max_min(abs(N - M), max(N, M), wnd);
			Matrix matrix(N, Vector(M, INF));
			matrix[0][0] = 0;
			for (int i = 0; i < N; i++) {
				for (int j = max(0, i - window); j < min(M, i + window + 1); j++) {
					double cost = INF;
					double dist = (*d_f)(*left[i], *right[j], nfunc);
					if (i == 0 && j == 0)
						matrix[i][j] = dist;
					else {
						if (i > 0)
							cost = matrix[i - 1][j];
						if (j > 0)
							cost = min(cost, matrix[i][j - 1]);
						if (i > 0 && j > 0)
							cost = min(cost, matrix[i - 1][j - 1]);
						matrix[i][j] = cost + dist;
					}
				}
			}

			{
				int i = N - 1;
				int j = M - 1;
				diff.push_back(matrix[N - 1][M - 1]);
				while ((i > 0) || (j > 0)) {
					if ((i > 0) && (j > 0)) {
						double m = min(matrix[i - 1][j], min(matrix[i][j - 1], matrix[i - 1][j - 1]));
						if (abs(m - matrix[i - 1][j - 1]) < EPSILON) {
							i--;
							j--;
						} else if (abs(m - matrix[i - 1][j]) < EPSILON) {
							i--;
						} else if (abs(m - matrix[i][j - 1]) < EPSILON) {
							j--;
						}
					} else if ((i > 0) && (j == 0)) {
						i--;
					} else if ((i == 0) && (j > 0)) {
						j--;
					}
					diff.push_back(matrix[i][j]);
				}
			}

			int nSize = diff.size();
			for (int i = 0; i < nSize - 1; i++)
				diff[i] -= diff[i + 1];
			std::sort(diff.begin(), diff.end());
			double med = 0;
			if (nSize % 2 == 0)
				med = (diff[nSize / 2] + diff[nSize / 2 + 1]) / 2.;
			else
				med = diff[nSize / 2];

			double avg = matrix[N - 1][M - 1] / (double)(nSize + 1);

			if (ret_median)
				return med * ratio;
			else
				return avg * ratio;
		}

		// Histogram intersection
		static double hist_cmp_intersect(const Vector& left, const Vector& right, bool norm = false)
		{
			double sum_l = 0, sum_r = 0;
			if (norm) {
				for (size_t i = 0; i < left.size(); i++) {
					sum_l += left[i];
					sum_r += right[i];
				}
				sum_l = max(sum_l, EPSILON);
				sum_r = max(sum_r, EPSILON);
			}
			double diff = 0;
			for (size_t i = 0; i < left.size(); i++)
				if (norm)
					diff += min(left[i] / sum_l, right[i] / sum_r);
				else
					diff += min(left[i], right[i]);

			return max_min(0., 1., 1. - diff);
		}

		// Ramer–Douglas–Peucker algorithm
		static int rdp_algorithm(const Points& input, Points& output, double eps)
		{
			if (input.empty())
				return 0;
			output.resize(0);
			Points points;
			points.resize(0);
			Point anchor = input.front();
			anchor.index = 0;
			Point floater = input.back();
			floater.index = input.size() - 1;
			output.push_back(anchor);
			points.push_back(floater);
			while (!points.empty()) {
				double max_dist = 0;
				Point farthest = anchor;
				for (int i = anchor.index + 1; i < floater.index; i++) {
					double dist = input[i].coord.distance_to_line(anchor.coord, floater.coord);
					if (dist > max_dist) {
						max_dist = dist;
						farthest = input[i];
						farthest.index = i;
					}
				}
				if (max_dist <= eps) {
					output.push_back(points.back());
					points.pop_back();
					anchor = floater;
					if (!points.empty()) {
						floater = points.back();
						floater.index = points.back().index;
					}
				} else {
					floater = farthest;
					points.push_back(floater);
				}
			}

			return 0;
		}

		// draw virtual plot
		static void plot(const int x, const int y, Raster * pr)
		{
			pr->set(x, y);
		}

		// Bresenham's Line Algorithm
		static void bresenham_algorithm(const Coord& a, const Coord& b, Raster * pr)
		{
			int ax = (int)a.X, ay = (int)a.Y;
			int bx = (int)b.X, by = (int)b.Y;
			int delta_x(bx - ax);
			signed char const ix((delta_x > 0) - (delta_x < 0));
			delta_x = abs(delta_x) << 1;
			int delta_y(by - ay);
			signed char const iy((delta_y > 0) - (delta_y < 0));
			delta_y = abs(delta_y) << 1;
			plot(ax, ay, pr);
			if (delta_x >= delta_y)
			{
				int error(delta_y - (delta_x >> 1));
				while (ax != bx) {
					if ((error >= 0) && (error || (ix > 0))) {
						error -= delta_x;
						ay += iy;
					}
					error += delta_y;
					ax += ix;
					plot(ax, ay, pr);
				}
			} else {
				int error(delta_x - (delta_y >> 1));
				while (ay != by) {
					if ((error >= 0) && (error || (iy > 0))) {
						error -= delta_y;
						ax += ix;
					}
					error += delta_x;
					ay += iy;
					plot(ax, ay, pr);
				}
			}
		}

		struct Point
		{
			struct PointFeatures
			{
				struct Context
				{
					Matrix data;

					void init(int row, int col)
					{
						data.resize(row, Vector(col, 0));
					}

					void normalize()
					{
						double sum = 0;
						for (size_t i = 0; i < data.size(); i++)
							for (size_t j = 0; j < data[i].size(); j++)
								sum += data[i][j];
						sum = max(sum, EPSILON);
						for (size_t i = 0; i < data.size(); i++)
							for (size_t j = 0; j < data[i].size(); j++)
								data[i][j] /= sum;
					}

					double compare(const Matrix& m) const
					{
						double diff = 0;
						for (size_t i = 0; i < data.size(); i++)
							for (size_t j = 0; j < data[i].size(); j++)
								diff += min(data[i][j], m[i][j]);

						return max_min(0., 1., 1. - diff);
					}
				};

				double curvature;
				double cluster;
				double center_dist;
				double path_index;
				double track_index;
				double speed;
				double density;
				double neighbours;
				double H_proj, V_proj;
				int begin_flag, end_flag;
				double extrema;
				Context context;
				Vector list;
				std::string nb_descriptor;

				PointFeatures()
				{
					curvature = 0;
					cluster = 0;
					center_dist = 0;
					path_index = 0;
					track_index = 0;
					speed = 0;
					density = 0;
					neighbours = 0;
					H_proj = 0, V_proj = 0;
					begin_flag = 0, end_flag = 0;
					extrema = 0;
				}

				double dist(const Vector& l) const
				{
					double diff = 0;
					for (size_t j = 0; j < list.size(); j++)
						diff += calc_dist(list[j], l[j]);

					return diff;
				}
			};

			PointPtr * self;
			int index;
			Coord coord;
			PointFeatures features;

			Point()
			{
				index = 0;
				self = 0;
			}

			Point(double x, double y, int i = 0)
			{
				coord.X = x;
				coord.Y = y;
				index = i;
				self = 0;
			}

			double distance_to_features(const Point& P, int nfunc = 0) const
			{
				double sim_context = features.context.compare(P.features.context.data);
				double sim_features = features.dist(P.features.list);

				return 0.5 * sim_context + 0.05 * sim_features;
			}
		};

		struct Segment
		{
			Points points;
			Points points_reduced;

			Segment()
			{
			}
		};

		struct Curve
		{
			struct CurveFeatures
			{
				double slope_angle;
				double speed;
				double negative_X, negative_Y;
				double positive_X, positive_Y;
				double segments_dist;
				double width, height;
				double ratio;
				double left, right, top, bottom;
				int points_count;
				int points_reduced_count;
				int points_count_bm;
				double points_ratio;
				Vector H_proj, V_proj;
				Matrix points_features;

				CurveFeatures()
				{
					slope_angle = 0;
					speed = 0;
					negative_X = 0, negative_Y = 0;
					positive_X = 0, positive_Y = 0;
					left = 0, right = 0;
					top = 0, bottom = 0;
					segments_dist = 0;
					width = 0, height = 0;
					ratio = 0;
					points_count = 0;
					points_reduced_count = 0;
					points_count_bm = 0;
					points_ratio = 0;
				}
			};

			PointsPtr path;
			PointsPtr path_reduced;
			Segments segments;
			CurveFeatures features;
			Raster raster;

			FuzzyCMeans::FCM fuzzy;
			NaiveBayes::NBC nb;

			Curve()
			{
			}

			void synchronize_points(bool write_reduced = true)
			{
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					if (write_reduced) {
						*pPoint = **pPoint->self;
					}
					else
						**pPoint->self = *pPoint;
				}
			}

			void calc_bounds()
			{
				features.left = 0, features.right = 0;
				features.top = 0, features.bottom = 0;
				for (size_t i = 0; i < 1; i++) {
					Point * pPoint = path[i];
					features.left = pPoint->coord.X;
					features.right = pPoint->coord.X;
					features.bottom = pPoint->coord.Y;
					features.top = pPoint->coord.Y;
				}
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					if (pPoint->coord.X < features.left)
						features.left = pPoint->coord.X;
					if (pPoint->coord.X > features.right)
						features.right = pPoint->coord.X;
					if (pPoint->coord.Y < features.bottom)
						features.bottom = pPoint->coord.Y;
					if (pPoint->coord.Y > features.top)
						features.top = pPoint->coord.Y;
				}
				features.width = max(features.right - features.left, EPSILON);
				features.height = max(features.top - features.bottom, EPSILON);
				features.ratio = features.width / features.height;
			}

			void calc_bitmap(bool sparse_flag = false)
			{
				raster.init(sparse_flag);
				for (size_t i = 0; i < path.size() - 1; i++) {
					Point * pPoint_l = path[i];
					Point * pPoint_r = path[i + 1];
					Coord l = curve_to_bitmap_coord(pPoint_l->coord);
					Coord r = curve_to_bitmap_coord(pPoint_r->coord);
					bresenham_algorithm(l, r, &raster);
				}

				calc_bm_points();
			}

			void calc_path()
			{
				path.clear();
				for (size_t i = 0; i < segments.size(); i++)
					for (size_t j = 0; j < segments[i].points.size(); j++)
						path.push_back(&segments[i].points[j]);
				for (size_t i = 0; i < path.size(); i++)
					path[i]->self = &path[i];
				calc_bounds();
			}

			int calc_points(bool recalc = true)
			{
				if (!recalc)
					return features.points_count;

				features.points_count = 0;
				for (size_t i = 0; i < segments.size(); i++) {
					int nSize = segments[i].points.size();
					if (nSize > 0) {
						segments[i].points[0].features.begin_flag = 1;
						segments[i].points[nSize - 1].features.end_flag = 1;
						features.points_count += nSize;
					}
				}
				calc_path();

				return features.points_count;
			}

			void calc_reduced_path()
			{
				synchronize_points(false);
				path_reduced.clear();
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++)
						path_reduced.push_back(&segments[i].points_reduced[j]);
				}
			}

			int calc_reduced_points(bool recalc = true)
			{
				if (!recalc)
					return features.points_reduced_count;

				features.points_reduced_count = 0;
				for (size_t i = 0; i < segments.size(); i++)
					features.points_reduced_count += segments[i].points_reduced.size();
				calc_reduced_path();

				return features.points_reduced_count;
			}

			void norm_vec(Vector * v, bool oper_max = false) {
				double max_v = 0;
				double sum_v = 0;
				for (size_t i = 0; i < v->size(); i++) {
					if (oper_max)
						max_v = max(max_v, abs((*v)[i]));
					else
						sum_v += abs((*v)[i]);
				}
				max_v = max(max_v, EPSILON);
				sum_v = max(sum_v, EPSILON);
				for (size_t i = 0; i < v->size(); i++) {
					if (oper_max)
						(*v)[i] /= max_v;
					else
						(*v)[i] /= sum_v;
				}
			}

			int calc_bm_points(bool recalc = true)
			{
				if (!recalc)
					return features.points_count_bm;

				raster.fill_coords();
				features.points_count_bm = raster.coords.size();
				features.H_proj.clear();
				features.H_proj.resize(raster.size + 1, 0);
				features.V_proj.clear();
				features.V_proj.resize(raster.size + 1, 0);
				for (size_t i = 0; i < raster.coords.size(); i++) {
					Coord pt;
					if (raster.sparse_flag)
						pt = curve_to_bitmap_coord(raster.coords[i], true);
					else
						pt = raster.coords[i];
					features.H_proj[(int)pt.X]++;
					features.V_proj[(int)pt.Y]++;
				}
				norm_vec(&features.H_proj);
				norm_vec(&features.V_proj);

				return features.points_count_bm;
			}

			Coord calc_median_speed()
			{
				Coord med;
				Vector X_dist, Y_dist;
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size() - 1; j++) {
						Point * pPoint_l = &segments[i].points[j];
						Point * pPoint_r = &segments[i].points[j + 1];
						double dX = abs(pPoint_r->coord.X - pPoint_l->coord.X);
						double dY = abs(pPoint_r->coord.Y - pPoint_l->coord.Y);
						X_dist.push_back(dX);
						Y_dist.push_back(dY);
					}
				}
				std::sort(X_dist.begin(), X_dist.end());
				std::sort(Y_dist.begin(), Y_dist.end());
				int nSize = X_dist.size();
				if (nSize % 2 == 0) {
					med.X = (X_dist[nSize / 2] + X_dist[nSize / 2 + 1]) / 2.;
					med.Y = (Y_dist[nSize / 2] + Y_dist[nSize / 2 + 1]) / 2.;
				} else {
					med.X = X_dist[nSize / 2];
					med.Y = Y_dist[nSize / 2];
				}
				features.speed = med.distance_to_point(Coord());

				return med;
			}

			void calc_points_speed()
			{
				for (size_t i = 0; i < segments.size(); i++) {
					int nSize = segments[i].points.size();
					if (nSize > 0) {
						segments[i].points[0].features.extrema = 1.;
						segments[i].points[nSize - 1].features.extrema = 1.;
					}
					for (size_t j = 1; j < segments[i].points.size() - 1; j++) {
						Point * pPoint_l = &segments[i].points[j - 1];
						Point * pPoint = &segments[i].points[j];
						Point * pPoint_r = &segments[i].points[j + 1];
						double spd_from = pPoint->coord.distance_to_point(pPoint_l->coord);
						double spd_to = pPoint->coord.distance_to_point(pPoint_r->coord);
						pPoint->features.speed = (spd_to + spd_from) / 2.;
						double ex = (pPoint_r->coord.X - pPoint->coord.X) * (pPoint->coord.X - pPoint_l->coord.X);
						double ey = (pPoint_r->coord.Y - pPoint->coord.Y) * (pPoint->coord.Y - pPoint_l->coord.Y);
						if (ex < 0 && ey < 0)
							pPoint->features.extrema = 1.;
						else if (ex < 0 || ey < 0)
							pPoint->features.extrema = .5;
					}
				}
				synchronize_points();
			}

			void calc_path_index()
			{
				if (features.points_count == 0)
					return;

				int nCur = 0;
				double nSize = (double)features.points_count;
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					pPoint->features.path_index = (tanh(((double)nCur / nSize - 0.5) * exp(1)) + 1.) / 2.;
					nCur++;
				}
				synchronize_points();
			}

			void calc_track_index(int tail = 5)
			{
				if (tail < 1)
					return;
				double length = 0;
				for (size_t i = 0; i < path.size() - 1; i++) {
					Point * pPoint_l = path[i];
					Point * pPoint = path[i + 1];
					length += pPoint->coord.distance_to_point(pPoint_l->coord);
					pPoint->features.track_index = length;
				}
				double max_length = length;
				for (int i = path.size() - 1; i >= tail ; i--) {
					Point * pPoint = path[i];
					Point * pPoint_tail = path[i - tail];
					length = pPoint->features.track_index - pPoint_tail->features.track_index;
					pPoint->features.track_index = length;
				}
				max_length = max(max_length, EPSILON);
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					pPoint->features.track_index /= max_length;
				}
				synchronize_points();
			}

			void calc_density(double radius = .05)
			{
				double rad = sqrt(features.width * features.height) * radius;
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					for (size_t j = 0; j < path.size(); j++) {
						if (pPoint->coord.distance_to_point(path[j]->coord) < rad)
							pPoint->features.density++;
					}
					pPoint->features.density /= (double)(features.points_count + 1.);
				}
				synchronize_points();
			}

			int simplify_path(double speed_koef = .25)
			{
				double spd = features.speed * speed_koef;
				for (size_t i = 0; i < segments.size(); i++)
					rdp_algorithm(segments[i].points, segments[i].points_reduced, spd);
		
				return calc_reduced_points();
			}

			Coord calc_points_center(CoordsPtr pts)
			{
				Coord pt;
				for (size_t i = 0; i < pts.size(); i++) {
					pt.X += pts[i]->X;
					pt.Y += pts[i]->Y;
				}
				double nPoints = pts.size();
				if (nPoints > 0) {
					pt.X /= nPoints;
					pt.Y /= nPoints;
				}

				return pt;
			}

			void move_to(Coord pt)
			{
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					pPoint->coord.X -= pt.X;
					pPoint->coord.Y -= pt.Y;
				}
				synchronize_points();
				calc_bounds();
			}

			void rotate(double angle)
			{
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					double X_r = pPoint->coord.X * cos(angle) + pPoint->coord.Y * sin(angle);
					double Y_r = -pPoint->coord.X * sin(angle) + pPoint->coord.Y * cos(angle);
					pPoint->coord.X = X_r;
					pPoint->coord.Y = Y_r;
				}
				synchronize_points();
				calc_bounds();
			}

			void scale()
			{
				double scl_x = max(abs(features.left), abs(features.right));
				double scl_y = max(abs(features.top), abs(features.bottom));
				double scl = max(scl_x, scl_y);
				scl = max(scl, EPSILON);
				for (size_t i = 0; i < path.size(); i++) {
					Point * pPoint = path[i];
					pPoint->coord.X /= scl;
					pPoint->coord.Y /= scl;
					pPoint->features.speed /= scl;
				}
				features.speed /= scl;
				synchronize_points();
				calc_bounds();
			}

			void calc_projs()
			{
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					Coord bm_pt = curve_to_bitmap_coord(pPoint->coord, true);
					pPoint->features.H_proj = features.H_proj[(int)bm_pt.X];
					pPoint->features.V_proj = features.V_proj[(int)bm_pt.Y];
				}
				synchronize_points(false);
			}

			double calc_neighbours_count(Coord pt, int radius = 3)
			{
				int bm_size = BITMAP_SIZE;
				Coord bm_pt = curve_to_bitmap_coord(pt);
				int x = (int)bm_pt.X;
				int y = (int)bm_pt.Y;

				double nn = 0;
				for (int i = 0; i < 2 * radius + 1; i++)
				{
					int xi = x + i - radius;
					if (xi < 0 || xi > bm_size)
						continue;
					for (int j = 0; j < 2 * radius + 1; j++)
					{
						int yj = y + j - radius;
						if (yj < 0 || yj > bm_size)
							continue;
						if (raster.is_xy(xi, yj))
							nn++;
					}
				}
				nn /= (double)pow(2 * radius + 1, 2);

				return nn;
			}

			void calc_neighbours()
			{
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					pPoint->features.neighbours = calc_neighbours_count(pPoint->coord);
				}
				synchronize_points(false);
			}

			void calc_curvatures()
			{
				for (size_t i = 1; i < path_reduced.size() - 1; i++) {
					Point * pPoint_l = path_reduced[i - 1];
					Point * pPoint = path_reduced[i];
					Point * pPoint_r = path_reduced[i + 1];
					pPoint->features.curvature = calc_angle_abc(pPoint_l->coord, pPoint->coord, pPoint_r->coord);
				}
				synchronize_points(false);
			}

			void calc_center_dists()
			{
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					pPoint->features.center_dist = pPoint->coord.distance_to_point(Coord());
				}
				synchronize_points(false);
			}

			void calc_shape_context()
			{
				const int szAB = CONTEXT_ANGLE_BIN_SIZE;
				const int szDB = CONTEXT_DIST_BIN_SIZE;
				const int window = CONTEXT_WINDOW_SIZE;
				int nSize = path_reduced.size();
				for (int i = 0; i < nSize; i++) {
					Point * pPoint = path_reduced[i];
					pPoint->features.context.init(szAB, szDB);
					for (int j = 0; j < nSize; j++) {
						if (i == j || !(abs(i - j) < window || (nSize - j) < window))
							continue;
						Point * pPoint_to = path_reduced[j];
						PolarCoord pc = Coord(pPoint->coord, pPoint_to->coord).to_polar();
						int angle_bin = (int)max_min(0, szAB - 1, (int)(((double)szAB) * (pc.Theta / M_PI / 2. + 0.5)));
						int length_bin = (int)max_min(0, szDB - 1, (int)pc.R);
						pPoint->features.context.data[angle_bin][length_bin] += 1;
					}
				}
				for (int i = 0; i < nSize; i++) {
					Point * pPoint = path_reduced[i];
					pPoint->features.context.normalize();
				}
				synchronize_points(false);
			}

			void calc_neg_pos()
			{
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size() - 1; j++) {
						Point * pPoint_l = &segments[i].points[j];
						Point * pPoint_r = &segments[i].points[j + 1];
						double dX = pPoint_r->coord.X - pPoint_l->coord.X;
						double dY = pPoint_r->coord.Y - pPoint_l->coord.Y;
						if (dX > 0)
							features.positive_X += dX;
						else
							features.negative_X -= dX;
						if (dY > 0)
							features.positive_Y += dY;
						else
							features.negative_Y -= dY;
					}
				}
				double dX = features.positive_X + features.negative_X;
				dX = max(dX, EPSILON);
				features.positive_X /= dX;
				features.negative_X /= dX;
				double dY = features.positive_Y + features.negative_Y;
				dY = max(dY, EPSILON);
				features.positive_Y /= dY;
				features.negative_Y /= dY;
				features.points_ratio = max((double)features.points_count_bm, EPSILON) / max((double)features.points_count, EPSILON);
			}

			void calc_segments_dist()
			{
				for (size_t i = 0; i < segments.size() - 1; i++) {
					Point * pPoint_l = &segments[i].points[segments[i].points.size() - 1];
					Point * pPoint_r = &segments[i + 1].points[0];
					features.segments_dist += pPoint_l->coord.distance_to_point(pPoint_r->coord);
				}
				if (segments.size() > 0)
					features.segments_dist /= segments.size();
			}

			Coord bitmap_to_curve_coord(Coord bm_pt, bool force_scale = false)
			{
				if (raster.sparse_flag && !force_scale)
					return bm_pt;

				double bm_size = BITMAP_SIZE;
				double center = BITMAP_SIZE / 2.;
				Coord pt;
				double max_wh = max(features.width, features.height);
				pt.X = ((bm_pt.X - .5) / bm_size - .5) * max_wh + features.width / 2. + features.left;
				pt.Y = ((bm_pt.Y - .5) / bm_size - .5) * max_wh + features.height / 2. + features.bottom;

				return pt;
			}

			Coord curve_to_bitmap_coord(Coord pt, bool force_scale = false)
			{
				if (raster.sparse_flag && !force_scale)
					return pt;

				double bm_size = BITMAP_SIZE;
				Coord bm_pt;
				double max_wh = max(features.width, features.height);
				bm_pt.X = bm_size * ((pt.X - features.left - features.width / 2.) / max_wh + .5) + .5;
				bm_pt.Y = bm_size * ((pt.Y - features.bottom - features.height / 2.) / max_wh + .5) + .5;

				return bm_pt;
			}

			void calc_clusters()
			{
				Matrix points;
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					Vector v;
					v.push_back(pPoint->features.speed);
					points.push_back(v);
				}
				fuzzy.fcm(points, 4);
				for (size_t i = 0; i < path_reduced.size(); i++)
					path_reduced[i]->features.cluster = fuzzy.clusters[i];
				synchronize_points(false);
			}

			void norm_features()
			{
				if (features.points_features.size() == 0)
					return;

				int num_features = features.points_features[0].size();
				Vector f_min(num_features), f_max(num_features);
				for (size_t i = 0; i < features.points_features.size(); i++) {
					for (int j = 0; j < num_features; j++) {
						f_min[j] = min(f_min[j], features.points_features[i][j]);
						f_max[j] = max(f_max[j], features.points_features[i][j]);
					}
				}
				Vector f_dist(num_features);
				for (int j = 0; j < num_features; j++)
					f_dist[j] = max(f_max[j] - f_min[j], EPSILON);
				for (size_t i = 0; i < features.points_features.size(); i++) {
					for (int j = 0; j < num_features; j++) {
						features.points_features[i][j] -= f_min[j];
						features.points_features[i][j] /= f_dist[j];
					}
				}
			}

			void calc_features()
			{
				features.points_features.clear();
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					Vector v;
					// base features
					v.push_back(tan(pPoint->features.curvature));
					v.push_back(pPoint->features.extrema);
					v.push_back(pPoint->features.begin_flag);
					v.push_back(pPoint->features.end_flag);
					v.push_back(pPoint->features.center_dist);
					v.push_back(pPoint->features.path_index);
					v.push_back(pPoint->features.track_index);
					v.push_back(pPoint->features.density);
					v.push_back(pPoint->features.neighbours);
					v.push_back(pPoint->features.H_proj);
					v.push_back(pPoint->features.V_proj);
					v.push_back(pPoint->features.speed);
					// some synthetic features
					v.push_back(tan(pPoint->features.curvature) * pPoint->features.center_dist);
					v.push_back(tan(pPoint->features.curvature) * pPoint->features.speed);
					v.push_back(tan(pPoint->features.curvature) * pPoint->features.density);
					v.push_back(tan(pPoint->features.curvature) * pPoint->features.neighbours);
					v.push_back(pPoint->features.path_index * pPoint->features.speed);
					v.push_back(pPoint->features.path_index * pPoint->features.track_index);
					v.push_back(pPoint->features.density * pPoint->features.neighbours);
					v.push_back(pPoint->features.density * pPoint->features.speed);
					v.push_back(pPoint->features.H_proj * pPoint->features.V_proj);
					v.push_back(pPoint->features.H_proj * pPoint->features.speed);
					v.push_back(pPoint->features.V_proj * pPoint->features.speed);
					v.push_back(pPoint->features.track_index * pPoint->features.density);
					v.push_back(pPoint->features.track_index * pPoint->features.neighbours);

					pPoint->features.list = v;
					features.points_features.push_back(v);
				}
				norm_features();
				synchronize_points(false);
			}

			void calc_nb_data()
			{
				if (features.points_features.size() == 0)
					return;

				int num_features = features.points_features[0].size();
				for (size_t i = 0; i < path_reduced.size(); i++) {
					Point * pPoint = path_reduced[i];
					std::string s_c = int_str(pPoint->coord.get_quadrant(), "Q");
					for (int j = 0; j < num_features; j++) {
						double f = features.points_features[i][j];
						std::string s_w = int_str(j, "F") + int_str((int)(f * 10.), "R");
						pPoint->features.nb_descriptor += s_w;
						if (j < num_features - 1)
							pPoint->features.nb_descriptor += " ";
						nb.nb_data[s_c][s_w] += 1;
						nb.nb_C[s_c] += 1;
						nb.nb_W[s_w] += 1;
					}
				}
				synchronize_points(false);
			}

			void center_rotate()
			{
				calc_bitmap();
				CoordsPtr pts, pts_left, pts_right;
				for (size_t i = 0; i < raster.coords.size(); i++)
					pts.push_back(&raster.coords[i]);
				Coord pt_center = calc_points_center(pts);
				for (size_t i = 0; i < raster.coords.size(); i++) {
					if (raster.coords[i].X <= pt_center.X)
						pts_left.push_back(&raster.coords[i]);
					if (raster.coords[i].X >= pt_center.X)
						pts_right.push_back(&raster.coords[i]);
				}
				Coord pt_left = calc_points_center(pts_left);
				Coord pt_right = calc_points_center(pts_right);
				features.slope_angle = calc_angle_xy(pt_left, pt_right);
				rotate(features.slope_angle);
				move_to(bitmap_to_curve_coord(pt_center));
				calc_bitmap();
			}

			void normalize()
			{
				calc_points();
				center_rotate();
				calc_median_speed();
				simplify_path();
			}

			void calc()
			{
				normalize();
				calc_points_speed();
				calc_path_index();
				calc_density();
				calc_segments_dist();
				calc_neg_pos();
				calc_track_index();
				calc_curvatures();
				calc_neighbours();
				calc_center_dists();
				calc_shape_context();
				calc_projs();
				calc_features();
				calc_nb_data();
			}
		};

		Curve curve;
	};
}
