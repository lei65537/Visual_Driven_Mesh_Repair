#pragma once

#include "typeBinary.h"
#include <stdio.h>
#include "MRFEnergy.h"
#include <vector>
class MRF_solver {
private:
	MRFEnergy<TypeBinary>* mrf;
	MRFEnergy<TypeBinary>::NodeId* nodes;
	MRFEnergy<TypeBinary>::Options options;
	TypeBinary::REAL energy, lowerBound;

	int _nodeNum;
public:
	MRF_solver(int nodeNum){
		mrf = new MRFEnergy<TypeBinary>(TypeBinary::GlobalSize());
		nodes = new MRFEnergy<TypeBinary>::NodeId[nodeNum];
		_nodeNum = nodeNum;
	}

	void add_nodes(int idx, double val1, double val2) {
		nodes[idx] = mrf->AddNode(TypeBinary::LocalSize(), TypeBinary::NodeData(val1, val2)); // scores for label 0,1
	}

	void add_edge(int idv0, int idv1, double cost00, double cost01, double cost10, double cost11) {
		mrf->AddEdge(nodes[idv0], nodes[idv1], TypeBinary::EdgeData(cost00, cost01, cost10, cost11));
	}

	void solve(std::vector<int> &res) {
		mrf->SetAutomaticOrdering();
		options.m_iterMax = 100;

#ifndef  _OUTPUT_MRF_LOG
		options.m_printMinIter = 100;
#endif
		mrf->Minimize_TRW_S(options, lowerBound, energy);
		res.clear();
		for (int i = 0; i < _nodeNum; i++)
		{
			res.push_back(mrf->GetSolution(nodes[i]));
		}
	}

	~MRF_solver() {
		delete[] nodes;
		delete mrf;
	}
};
