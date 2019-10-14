#pragma once

#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;

class  SparseTensor
{
private:
	std::vector<SpMat> _val;
	std::map<int, int> _vid;//last dim, access index to inner storage id
	int _size;

public:
	SparseTensor() {}
	~SparseTensor() {}

	inline int update_minimal_siz() { _size = _val.size(); return _size; }
	inline void setSize(int siz) { _size = siz; }
	inline int getSize(){return _size;}
	inline void insert(SpMat &mat, int id) {
		//assert(id < _size);
		std::pair<std::map<int, int>::iterator, bool> ret;
		ret = _vid.insert(std::pair<int, int>(id, _vid.size()));
		if (ret.second == false) {
			//already exist, should add matirx to coresponding one
			int currectid = ret.first->second;
			_val[currectid] += mat;
		}
		else {
			_val.push_back(mat);

		}
	}
	inline SpMat* getElement(int i) { if (_vid.find(i) == _vid.end()) return NULL;  return &(_val[_vid[i]]); }//have not check wheter the idx is valid
	SparseTensor(std::vector<SpMat> &vmat, std::vector<int> &id, int size)//work as triplets
	{
		_size = size;

		for (int i = 0; i < id.size(); i++) {
			
			insert(vmat[i], id[i]);
		}
		
	}

	inline SparseTensor mul(SpMat &m) {
		std::vector<SpMat> mats;
		std::vector<int> ids;
		for (int k = 0; k<m.outerSize(); ++k)
			for (SpMat::InnerIterator it(m, k); it; ++it)
			{
				double vv=it.value();
				int m=it.row();   // row index
				int n=it.col();   // col index (here it is equal to k)
				//it.index(); // inner index, here it is equal to it.row()
				mats.push_back(_val[m] * vv);
				ids.push_back(n);
			}

		return SparseTensor(mats, ids, m.cols());
	}
};

