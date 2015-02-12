#ifndef MIP_H
#define MIP_H

#include <vector>
#include <ilocplex.h>

#include "data.h"

using namespace std;

//the MIP class is the class for the mixed integer programming approach to solve a given instance exactly
class MIP{

public:
	//the class is initialized with a certain instance
	MIP(Instance inst);
	//the solve method depends on the number of vehicles M
	double solve(int M);
	//here, the solution time of the MIP is returned
	double getSolutionTime(){return solutionTime;};

private:
	vector<Node> nodes;
	int m;
	double solutionTime;
};

#endif