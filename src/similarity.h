/**
* Copyright © 2015-2016 Stanislav Semenov. All rights reserved.
* Contacts: <stas.semenov@gmail.com>
*
* This C++ library implements algorithm for computing similarity between curves.
*/

#pragma once

#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define EPSILON (1E-15)
#define INF (std::numeric_limits<double>::infinity())
#define CONTEXT_ANGLE_BIN_SIZE 4
#define CONTEXT_DIST_BIN_SIZE 6
#define CONTEXT_WINDOW_SIZE 32

namespace Similarity
{
	using std::min;
	using std::max;

	static double max_min(double left, double right, double value)
	{
		return max(left, min(right, value));
	}

	class SmartCurve
	{
	public:
		struct Point;
		struct Segment;
		struct Curve;

		typedef std::vector<double> Vector;
		typedef std::vector<Vector> Matrix;
		typedef std::vector<Point> Points;
		typedef std::vector<Segment> Segments;

		static double calc_angle_abc(const Point& a, const Point& b, const Point& c)
		{
			double abX = b.X - a.X;
			double abY = b.Y - a.Y;
			double cbX = b.X - c.X;
			double cbY = b.Y - c.Y;

			double dot = abX * cbX + abY * cbY;
			double cross = abX * cbY - abY * cbX;

			return atan2(cross, dot);
		}

		static double calc_angle_xy(const Point& X, const Point& Y)
		{
			double dX = Y.X - X.X;
			double dY = Y.Y - X.Y;

			return atan2(dY, dX);
		}

		static double distance_to_features(const Point& left, const Point& right)
		{
			return left.distance_to_features(right);
		}

		// Dynamic time warping algorithm
		static double dtw_simple_algorithm(const Points& left, const Points& right, int wnd = 0, double (*d_f)(const Point&, const Point&) = distance_to_features)
		{
			int N = left.size();
			int M = right.size();
			int window = (int)max_min(abs(N - M), max(N, M), wnd);
			double cost = INF;

			Matrix matrix(N, Vector(M, INF));
			matrix[0][0] = 0;

			for (int i = 0; i < N; i++) {
				for (int j = max(0, i - window); j < min(M, i + window + 1); j++) {
					cost = INF;
					if (i == 0 && j == 0)
						matrix[i][j] = (*d_f)(left[i], right[j]);
					else {
						if (i > 0)
							cost = matrix[i - 1][j];
						if (j > 0)
							cost = min(cost, matrix[i][j - 1]);
						if (i > 0 && j > 0)
							cost = min(cost, matrix[i - 1][j - 1]);
						matrix[i][j] = cost + (*d_f)(left[i], right[j]);
					}
				}
			}

			return 2 * matrix[N - 1][M - 1] / (double)(N + M + 1);
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
					double dist = input[i].distance_to_line(anchor, floater);
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

		struct Point
		{
			struct PointFeatures
			{
				struct Contexts
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
								sum += abs(data[i][j]);

						if (sum > 0) {
							for (size_t i = 0; i < data.size(); i++)
								for (size_t j = 0; j < data[i].size(); j++)
									data[i][j] /= sum;
						}
					}
				};

				double curvature;
				double center_dist;
				double path_index;
				double speed;
				double density;
				Contexts context;
			};

			double X;
			double Y;
			int index;
			PointFeatures features;

			Point()
			{
				X = 0;
				Y = 0;
				index = 0;
			}

			Point(double x, double y, int i = 0)
			{
				X = x;
				Y = y;
				index = i;
			}

			double distance_to_features(const Point& P) const
			{
				if (index == P.index) {
					if (index == -1)
						return 0.;
				} else {
					if (index == -1 || P.index == -1)
						return 1.;
				}

				const double a_part = M_PI / CONTEXT_ANGLE_BIN_SIZE;
				double sim_context = 0;
				for (size_t i = 0; i < features.context.data.size(); i++) {
					double sim_angle = abs(sin(2. * a_part * (double)i - M_PI + a_part) + EPSILON) / (1. + EPSILON);
					for (size_t j = 0; j < features.context.data[i].size(); j++) {
						sim_context += abs(features.context.data[i][j] - P.features.context.data[i][j]) * sim_angle;
					}
				}
				sim_context /= 2.;

				Vector f;
				f.push_back(sim_context);
				f.push_back(abs(sin(features.curvature) - sin(P.features.curvature)) / 2.);
				f.push_back(abs(features.center_dist - P.features.center_dist) / sqrt(2));
				f.push_back(min(exp(abs(features.path_index - P.features.path_index)) - 1., 1.));
				f.push_back(min((exp(abs(features.speed - P.features.speed)) - 1.) / 4., 1.));
				f.push_back(min(exp(abs(features.density - P.features.density)) - 1., 1.));
				
				double res = std::accumulate(f.begin(), f.end(), 0.0) / (double)f.size();
				return res;
			}

			double distance_to_point(const Point& P) const
			{
				return sqrt(pow((X - P.X), 2) + pow((Y - P.Y), 2));
			}

			double distance_to_line(const Point& a, const Point& b) const
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
				double angle;
				double speed;

				CurveFeatures()
				{
					angle = 0;
					speed = 0;
				}
			};

			Segments segments;
			Points full_path;
			CurveFeatures features;
			int nTotalPoints;
			int nReducedPoints;

			Curve()
			{
				nTotalPoints = 0;
				nReducedPoints = 0;
			}

			int calc_total_points(bool recalc = true)
			{
				if (!recalc)
					return nTotalPoints;

				nTotalPoints = 0;
				for (size_t i = 0; i < segments.size(); i++)
					nTotalPoints += segments[i].points.size();

				return nTotalPoints;
			}

			int calc_reduced_points(bool recalc = true)
			{
				if (!recalc)
					return nReducedPoints;

				nReducedPoints = 0;
				for (size_t i = 0; i < segments.size(); i++)
					nReducedPoints += segments[i].points_reduced.size();

				return nReducedPoints;
			}

			Point calc_median_speed()
			{
				Vector X_dist, Y_dist;
				double X_med = 0, Y_med = 0;
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size() - 1; j++) {
						Point * pPoint_l = &segments[i].points[j];
						Point * pPoint_r = &segments[i].points[j + 1];
						X_dist.push_back(abs(pPoint_r->X - pPoint_l->X));
						Y_dist.push_back(abs(pPoint_r->Y - pPoint_l->Y));
					}
				}
				std::sort(X_dist.begin(), X_dist.end());
				std::sort(Y_dist.begin(), Y_dist.end());
				int nSize = X_dist.size();
				if (nSize % 2 == 0)
					X_med = (X_dist[nSize / 2] + X_dist[nSize / 2 + 1]) / 2.;
				else
					X_med = X_dist[nSize / 2];
				if (nSize % 2 == 0)
					Y_med = (Y_dist[nSize / 2] + Y_dist[nSize / 2 + 1]) / 2.;
				else
					Y_med = Y_dist[nSize / 2];
				features.speed = sqrt(pow(X_med, 2) + pow(Y_med, 2));

				return Point(X_med, Y_med);
			}

			void calc_points_speed()
			{
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 1; j < segments[i].points.size() - 1; j++) {
						Point * pPoint_l = &segments[i].points[j - 1];
						Point * pPoint = &segments[i].points[j];
						Point * pPoint_r = &segments[i].points[j + 1];
						double spd_from = pPoint->distance_to_point(*pPoint_l);
						double spd_to = pPoint->distance_to_point(*pPoint_r);
						pPoint->features.speed = (spd_to + spd_from) / 2.;
					}
				}
			}

			void calc_path_index()
			{
				int nCur = 0;
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						pPoint->features.path_index = (double)(nCur + 1) / (double)(nTotalPoints + 1);
						nCur++;
					}
				}
			}

			void calc_density(double radius_koef = 8.)
			{
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						for (size_t ii = 0; ii < segments.size(); ii++) {
							for (size_t jj = 0; jj < segments[ii].points.size(); jj++) {
								if (pPoint->distance_to_point(segments[ii].points[jj]) < features.speed * radius_koef)
									pPoint->features.density++;
							}
						}
						pPoint->features.density /= (double)(nTotalPoints + 1);
					}
				}
			}

			int simplify_path(double speed_koef = 4.)
			{
				for (size_t i = 0; i < segments.size(); i++)
					rdp_algorithm(segments[i].points, segments[i].points_reduced, features.speed / speed_koef);
		
				return calc_reduced_points();
			}

			Point calc_reduced_center(int lr = 0)
			{
				double nPoints = 0;
				Point pt;
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++) {
						Point * pPoint = &segments[i].points_reduced[j];
						if (lr == -1 && pPoint->X > 0)
							continue;
						if (lr == 1 && pPoint->X < 0)
							continue;
						pt.X += pPoint->X;
						pt.Y += pPoint->Y;
						nPoints++;
					}
				}
				if (nPoints > 0) {
					pt.X /= nPoints;
					pt.Y /= nPoints;
				}

				return pt;
			}

			void move_to(Point pt)
			{
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++) {
						Point * pPoint = &segments[i].points_reduced[j];
						pPoint->X -= pt.X;
						pPoint->Y -= pt.Y;
					}
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						pPoint->X -= pt.X;
						pPoint->Y -= pt.Y;
					}
				}
			}

			void rotate(double angle)
			{
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						double X_r = pPoint->X * cos(angle) - pPoint->Y * sin(angle);
						double Y_r = pPoint->X * sin(angle) + pPoint->Y * cos(angle);
						pPoint->X = X_r;
						pPoint->Y = Y_r;
					}
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++) {
						Point * pPoint = &segments[i].points_reduced[j];
						double X_r = pPoint->X * cos(angle) - pPoint->Y * sin(angle);
						double Y_r = pPoint->X * sin(angle) + pPoint->Y * cos(angle);
						pPoint->X = X_r;
						pPoint->Y = Y_r;
					}
				}
			}

			void scale()
			{
				double X_Y_max = 0;
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						X_Y_max = max(X_Y_max, abs(pPoint->X));
						X_Y_max = max(X_Y_max, abs(pPoint->Y));
					}
				}
				for (size_t i = 0; i < segments.size(); i++) {
					for (size_t j = 0; j < segments[i].points.size(); j++) {
						Point * pPoint = &segments[i].points[j];
						if (X_Y_max > 0) {
							pPoint->X /= X_Y_max;
							pPoint->Y /= X_Y_max;
							pPoint->features.speed /= X_Y_max;
						}
					}
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++) {
						Point * pPoint = &segments[i].points_reduced[j];
						if (X_Y_max > 0) {
							pPoint->X /= X_Y_max;
							pPoint->Y /= X_Y_max;
							pPoint->features.speed /= X_Y_max;
						}
					}
				}
			}

			void calc_full_path()
			{
				for (size_t i = 0; i < segments.size(); i++) {
					Point pt;
					pt.index = -1;
					full_path.push_back(pt);
					for (size_t j = 0; j < segments[i].points_reduced.size(); j++)
						full_path.push_back(segments[i].points_reduced[j]);
				}
			}

			void calc_curvatures()
			{
				for (size_t i = 1; i < full_path.size() - 1; i++) {
					Point * pPoint_l = &full_path[i - 1];
					Point * pPoint = &full_path[i];
					Point * pPoint_r = &full_path[i + 1];
					pPoint->features.curvature = calc_angle_abc(*pPoint_l, *pPoint, *pPoint_r);
				}
			}

			void calc_center_dists()
			{
				for (size_t i = 1; i < full_path.size() - 1; i++) {
					Point * pPoint = &full_path[i];
					if (pPoint->index == -1)
						continue;
					pPoint->features.center_dist = sqrt(pPoint->X * pPoint->X + pPoint->Y * pPoint->Y);
				}
			}

			void calc_shape_context()
			{
				int a_bsz = CONTEXT_ANGLE_BIN_SIZE;
				int d_bsz = CONTEXT_DIST_BIN_SIZE;
				int window = CONTEXT_WINDOW_SIZE;
				int nFullPoints = full_path.size();
				for (int i = 0; i < nFullPoints; i++) {
					Point * pPoint = &full_path[i];
					pPoint->features.context.init(a_bsz, d_bsz);
					if (pPoint->index == -1)
						continue;
					for (int j = 0; j < nFullPoints; j++) {
						if (i == j || !(abs(i - j) < window || (nFullPoints - j) < window))
							continue;
						Point * pPoint_to = &full_path[j];
						int angle_bin = (int)max_min(0, a_bsz - 1, (int)floor(a_bsz * ((M_PI + calc_angle_abc(*pPoint, Point(), *pPoint_to)) / M_PI / 2.)));
						int length_bin = (int)max_min(0, d_bsz - 1, (int)floor(d_bsz * (pPoint->distance_to_point(*pPoint_to) / sqrt(2.) / 2.)));
						pPoint->features.context.data[angle_bin][length_bin] += 1;
					}
				}
				for (int i = 0; i < nFullPoints; i++) {
					Point * pPoint = &full_path[i];
					pPoint->features.context.normalize();
				}
			}

			void normalize()
			{
				calc_total_points();
				calc_median_speed();
				calc_path_index();
				simplify_path();
				move_to(calc_reduced_center());
				features.angle = calc_angle_xy(calc_reduced_center(-1), calc_reduced_center(1));
				rotate(features.angle);
				scale();
				calc_points_speed();
				calc_density();
			}

			void calc()
			{
				normalize();
				calc_full_path();
				calc_curvatures();
				calc_center_dists();
				calc_shape_context();
			}
		};

		Curve curve;
	};
}
