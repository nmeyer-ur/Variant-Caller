/*
    Allele-specific variant caller
    Copyright (C) 2019  Nils Meyer, University of Regensburg, Germany
    <nils.meyer@ur.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/  // end legal

#ifndef __DBSCAN_H__
#define __DBSCAN_H__

// Original code written by: Ibraheem Alhashim / ialhashim
// https://gist.github.com/ialhashim/b29e5455333aa6ae0071 ,
// which is adapted from https://github.com/propanoid/DBSCAN

#include <vector>
#include <algorithm>
#include "ecdefs.hpp"
#include "Eigen/Core"

template<typename Vector, typename Matrix, typename Epsilon>
class DBSCAN {

public:
	typedef Matrix DistanceMatrix;
	typedef std::vector<unsigned int> Neighbors;
	typedef std::vector<int> Labels;

private:
	Epsilon m_eps;
	Labels m_labels;

public:
	DBSCAN(Epsilon eps): m_eps( eps ) {

		clear();
	}

	const Labels & get_labels() const {

		return m_labels;
	}

	void clear() {
		m_labels.clear();
	}

	void reinit(Epsilon eps) {
		m_eps = eps;
        clear();
	}

	void fit_precomputed(const DistanceMatrix& D) {

        m_labels.resize(D.rows(), -1);
		dbscan( D );
	}


private:
	Neighbors find_neighbors(const DistanceMatrix& D, unsigned int pid) {

		Neighbors ne;

		for (unsigned int j = 0; j < D.rows(); ++j)
			if 	(D(pid, j) <= m_eps)
				ne.push_back(j);

		return ne;
	}

	void dbscan(const DistanceMatrix & dm) {

		std::vector<unsigned int> visited(dm.rows());
		unsigned int cluster_id = 0;

		for (unsigned int pid = 0; pid < dm.rows(); ++pid) {
			if ( !visited[pid] ) {
				visited[pid] = 1;
				Neighbors ne = find_neighbors(dm, pid);

				if (ne.size() >= 1) {
					m_labels[pid] = cluster_id;

					for (unsigned int i = 0; i < ne.size(); ++i) {
						unsigned int nPid = ne[i];

						if (!visited[nPid]) {
							visited[nPid] = 1;
							Neighbors ne1 = find_neighbors(dm, nPid);

							if (ne1.size() >= 1)
								for (const auto& n1 : ne1)
								   ne.push_back(n1);
						}

						if (m_labels[nPid] == -1)
							m_labels[nPid] = cluster_id;
					}
					++cluster_id;
				}
			}
		}
	}
};

#endif // __DBSCAN_H__
