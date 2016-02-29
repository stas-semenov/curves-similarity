/*
*  Copyright © 2015-2016 Stanislav Semenov. All rights reserved.
*  Contacts: <stas.semenov@gmail.com>
*  https://github.com/stas-semenov
*
*  This C++ library implements Fuzzy c-means algorithm.
*
*  Source code is provided under the GNU General Public License v3.
*  To access the source code please visit the public code repository.
*/

#pragma once

#define NOMINMAX

#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

#define EPSILON (1E-15)
#define INF (std::numeric_limits<double>::infinity())

namespace FuzzyCMeans
{
	using std::min;
	using std::max;

	typedef std::vector<double> Vector;
	typedef std::vector<Vector> Matrix;
	typedef std::pair<double, int> ValueKey;
	typedef std::vector<ValueKey> ValueKeyVec;

	static double max_min(double left, double right, double value)
	{
		return max(left, min(right, value));
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

	class FCM
	{
	public:
		double fuzziness;
		int num_points;
		int num_clusters;
		int num_dim;

		Matrix u;
		Matrix centers;
		Matrix membership;
		std::vector<int> relations;
		std::vector<int> clusters;

		FCM()
		{
			std::srand(unsigned(std::time(0)));
			fuzziness = 2;
			num_points = 0;
			num_clusters = 2;
			num_dim = 1;
		}

		double norm(const Matrix& points, int i, int j)
		{
			double sum = 0;
			for (int k = 0; k < num_dim; k++)
				sum += calc_dist(points[i][k], centers[j][k]);

			return max(sum / (double)num_dim, EPSILON);
		}

		double calc_wij(const Matrix& points, int i, int j)
		{
			double sum = 0;
			double p = 2. / (fuzziness - 1.);
			for (int k = 0; k < num_clusters; k++) {
				double nj = norm(points, i, j);
				double nk = norm(points, i, k);
				sum += pow(nj / nk, p);
			}
			return 1.0 / sum;
		}

		double update_membership(const Matrix& points)
		{
			double max_diff = 0;
			for (int j = 0; j < num_clusters; j++) {
				for (int i = 0; i < num_points; i++) {
					double nms = calc_wij(points, i, j);
					double diff = nms - membership[i][j];
					membership[i][j] = nms;
					if (diff > max_diff)
						max_diff = diff;
				}
			}

			return max_diff;
		}

		int calc_centers(const Matrix& points)
		{
			for (int i = 0; i < num_points; i++) {
				for (int j = 0; j < num_clusters; j++) {
					u[i][j] = pow(membership[i][j], fuzziness);
				}
			}
			for (int j = 0; j < num_clusters; j++) {
				for (int k = 0; k < num_dim; k++) {
					double n = 0;
					double d = 0;
					for (int i = 0; i < num_points; i++) {
						n += u[i][j] * points[i][k];
						d += u[i][j];
					}
					centers[j][k] = n / max(d, EPSILON);
				}
			}

			return 0;
		}

		void fcm(const Matrix& points, int number_clusters = 2, double eps = .001, int max_iter = 30, double fuzz = 1.1)
		{
			num_points = points.size();
			if (num_points)
				num_dim = points[0].size();
			else
				return;

			if (number_clusters < 2 || fuzz < (1. + EPSILON))
				return;

			num_clusters = number_clusters;
			fuzziness = fuzz;
			u.resize(num_points, Vector(num_clusters));
			centers.resize(num_clusters, Vector(num_dim));
			membership.clear();
			membership.resize(num_points, Vector(num_clusters));
			for (int i = 0; i < num_points; i++) {
				double n = 0;
				for (int j = 0; j < num_clusters; j++)
					n += membership[i][j] = (rand() / (double)RAND_MAX);
				n = max(n, EPSILON);
				for (int j = 0; j < num_clusters; j++)
					membership[i][j] /= n;
			}
			double diff = 0;
			int count = 0;
			while (true) {
				calc_centers(points);
				diff = update_membership(points);
				count++;
				if (count > max_iter || diff < eps)
					break;
			}
			relations.clear();
			relations.resize(num_clusters);
			ValueKeyVec order;
			for (int j = 0; j < num_clusters; j++)
				order.push_back(ValueKey(centers[j][0], j));
			std::sort(order.begin(), order.end());
			for (int j = 0; j < num_clusters; j++)
				relations[order[j].second] = j;
			clusters.clear();
			for (int i = 0; i < num_points; i++) {
				int max_c = std::distance(membership[i].begin(), std::max_element(membership[i].begin(), membership[i].end()));
				clusters.push_back(relations[max_c]);
			}
			u.clear();
			centers.clear();
			membership.clear();
			relations.clear();
		}
	};
}
