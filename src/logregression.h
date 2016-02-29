/*
*  Copyright © 2015-2016 Stanislav Semenov. All rights reserved.
*  Contacts: <stas.semenov@gmail.com>
*  https://github.com/stas-semenov
*
*  This C++ library implements logistic regression with L1-regularization.
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
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace LogisticRegression
{
	using std::min;
	using std::max;

	static double max_min(double left, double right, double value)
	{
		return max(left, min(right, value));
	}

	typedef std::vector<std::string> Strings;
	typedef std::vector<double> Vector;
	typedef std::vector<Vector> Matrix;
	typedef std::map<int, double> IDMap;
	typedef std::vector<IDMap> IDMaps;

	static Strings split(const std::string& str, char delim) {
		Strings elems;
		std::stringstream stream(str);
		std::string item;
		while (getline(stream, item, delim))
			elems.push_back(item);

		return elems;
	}

	static double vec_norm(IDMap& left, IDMap& right)
	{
		double sum = 0;
		for (IDMap::iterator it = left.begin(); it != left.end(); it++)
			sum += pow(left[it->first] - right[it->first], 2);

		return sqrt(sum);
	}

	class LRegression
	{
	public:
		IDMap weights;
		IDMap total_l1;

		LRegression()
		{
			std::srand(unsigned(std::time(0)));
		}

		void load(const char * filename)
		{
			weights.clear();
			std::ifstream file;
			std::string line;
			file.open(filename);
			while (getline(file, line)) {
				if (line.length()) {
					Strings tokens = split(line, ' ');
					if (tokens.size() == 2)
						weights[atoi(tokens[0].c_str())] = atof(tokens[1].c_str());
				}
			}
			file.close();
		}

		void save(const char * filename)
		{
			std::ofstream file;
			file.open(filename);
			for (IDMap::iterator it = weights.begin(); it != weights.end(); it++)
				file << it->first << " " << it->second << std::endl;
			file.close();
		}

		void learn(IDMaps& examples, double eps = .005, int max_iter = 10000, double alpha = .001, double l1 = .0001)
		{
			weights.clear();
			total_l1.clear();

			for (size_t i = 0; i < examples.size(); i++) {
				for (IDMap::iterator it = examples[i].begin(); it != examples[i].end(); it++) {
					weights[it->first] = 0;
					total_l1[it->first] = 0;
				}
			}

			double mu = 0;
			double norm = 1;
			int count = 0;
			std::vector<int> index;
			for (size_t i = 0; i < examples.size(); i++)
				index.push_back(i);

			while (norm > eps) {
				IDMap old_weights(weights);
				std::random_shuffle(index.begin(), index.end());
				for (size_t i = 0; i < examples.size(); i++) {
					mu += l1 * alpha;
					double real = examples[index[i]][0];
					double predicted = classify(examples[index[i]]);
					for (IDMap::iterator it = examples[index[i]].begin(); it != examples[index[i]].end(); it++) {
						if (it->first != 0) {
							weights[it->first] += alpha * (real - predicted) * it->second;
							if (l1) {
								double z = weights[it->first];
								if (weights[it->first] > 0)
									weights[it->first] = max(0., (double)(weights[it->first] - (mu + total_l1[it->first])));
								else if (weights[it->first] < 0)
									weights[it->first] = min(0., (double)(weights[it->first] + (mu - total_l1[it->first])));
								total_l1[it->first] += weights[it->first] - z;
							}
						}
					}
				}
				norm = vec_norm(weights, old_weights);
				count++;
				if (count > max_iter)
					break;
			}
		}

		double classify(IDMap& features)
		{
			double logit = 0;
			for (IDMap::iterator it = features.begin(); it != features.end(); it++)
				if (it->first != 0)
					logit += it->second * weights[it->first];

			return 1. / (1. + exp(-max_min(-15, 15, logit)));
		}

	};
}
