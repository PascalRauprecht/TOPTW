#ifndef HSSA_H
#define HSSA_H

#include <vector>
#include <time.h>
#include <random>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "data.h"
#include "dynpro.h"

using namespace std;

//the sequence class is storing a sequence of nodes as well as the seqValue calculated in the dynamic programming
class Sequence{

public:
	Sequence(int Number):number(Number){};
	Sequence(vector<Node> Nodes, int Number, string ProblemType):nodes(Nodes),number(Number){determineSeqValue(ProblemType);};
	void setNodes(vector<Node> Nodes, string ProblemType);
	vector<Node> getNodes(){ return nodes;};
	double getSeqValue(){return seqValue;};

private:
	int number;
	vector<Node> nodes;
	double seqValue;
	void determineSeqValue(string ProblemType);
};

//the class SimA is the hybrid selective simulated annealing heuristic
class SimA{
public:
	//the class is initialized with the respective nodes of an instance
	SimA(Instance Inst):nodes(Inst.getNodes()), inst(Inst){};
	//this method is starting the heuristic
	double startSA(int M, int SAType, int Parameter, double T, double Alpha, int IIterB);
	double startSAOld(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue);
	//this method is returning the best found sequence in the heuristidc
	vector<Sequence> getSequences(){return sequences;};
	//this method is returning the computational time needed to execute the heuristic
	double getDurationOfSimAn(){return durationOfSimAn;};

private:
	Instance inst;
	vector<Node> nodes;
	int variant;
	int m, saType, tMax, nNonImproving, iIter, iSel, iIterB;
	double t, alpha, optimalValue;
	vector<int> stockedNodes;
	vector<int> selectedNodes;
	vector<int> xBest, stockedBest;
	double fBest;
	vector<Sequence> sequences;
	double durationOfSimAn;

	void insertZeroesInSequence();
	void checkVerticFeasibilityInSequence();

	void randomSelectionMethod();

	void swapSeq(int pos1, int pos2);
	void insertSeq(int pos1, int pos2);
	void reversionSeq(int pos1, int pos2);
	void exchangeSeq(int posSel, int posStock);

};

#endif