/*
*  Copyright © 2015-2016 Stanislav Semenov. All rights reserved.
*  Contacts: <stas.semenov@gmail.com>
*  https://github.com/stas-semenov
*
*  This C++ library implements Naive Bayes Classifier.
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
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

#define EPSILON (1E-15)
#define INF (std::numeric_limits<double>::infinity())

namespace NaiveBayes
{
	using std::min;
	using std::max;

	typedef std::map<std::string, std::map <std::string, int> > SSIMap;
	typedef std::map<std::string, int> SIMap;
	typedef std::map<std::string, double> SDMap;
	typedef std::vector<std::string> Strings;

	static Strings explode(const std::string &str, const std::string &delimiter)
	{
		Strings arr;
		int str_len = str.length();
		int del_len = delimiter.length();
		if (del_len == 0)
			return arr;

		int i = 0, k = 0;
		while (i < str_len) {
			int j = 0;
			while (i + j < str_len && j < del_len && str[i + j] == delimiter[j])
				j++;
			if (j == del_len) {
				arr.push_back(str.substr(k, i - k));
				i += del_len;
				k = i;
			}
			else
				i++;
		}
		arr.push_back(str.substr(k, i - k));

		return arr;
	}

	class NBC
	{
	public:
		SSIMap nb_data;
		SIMap nb_C;
		SIMap nb_W;
		SDMap nb_classes;

		NBC()
		{
		}

		void calc(const Strings& words)
		{
			const double coef = 3.;
			int cc = nb_data.size();
			SDMap pwl;
			for (size_t i = 0; i < words.size(); i++) {
				std::string w = words[i];
				SDMap clw;
				SDMap pw;
				double sum = 0;
				int wc = 0;
				if (nb_W.count(w))
					wc = nb_W[w];
				for (SIMap::iterator ic = nb_C.begin(), ice = nb_C.end(); ic != ice; ++ic) {
					if (nb_data.count(ic->first) && nb_data[ic->first].count(w))
						clw[ic->first] = nb_data[ic->first][w];
					sum += clw[ic->first] / ic->second;
				}
				for (SIMap::iterator ic = nb_C.begin(), ice = nb_C.end(); ic != ice; ++ic) {
					if (wc == 0)
						pw[ic->first] = 1. / cc;
					else
						pw[ic->first] = (coef / cc + wc * clw[ic->first] / ic->second / sum) / (coef + wc);
					pwl[ic->first] += log(pw[ic->first]);
				}
			}
			double div = 0;
			for (SIMap::iterator ic = nb_C.begin(), ice = nb_C.end(); ic != ice; ++ic)
				div += exp(pwl[ic->first]);
			for (SIMap::iterator ic = nb_C.begin(), ice = nb_C.end(); ic != ice; ++ic)
				nb_classes[ic->first] = exp(pwl[ic->first]) / div;
		}

		void calc(const std::string& words)
		{
			calc(explode(words, " "));
		}
	};
}
