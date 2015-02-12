#ifndef DYNPRO_H
#define DYNPRO_H

#include <string>
#include <vector>

#include "data.h"

//The class State is used in the dynamic programming class. Each state has a best next state, an optimal value and an optimal duration. The paper to this heuristic is providing further insights.
class State{
	
public:
	State(int T):t(T){};
	void setOptimalValue(double OptimalValue){optimalValue = OptimalValue;};
	double getOptimalValue(){return optimalValue;};

	void setOptimalStateInNextStage(double OptimalStateInNextStage){optimalStateInNextStage = OptimalStateInNextStage;};
	double getOptimalStateInNextStage(){return optimalStateInNextStage;};

	void setOptimalDuration(double OptimalDuration){optimalDuration = OptimalDuration;};
	double getOptimalDuration(){return optimalDuration;};

	int getT(){return t;};

private:
	//t is the time of the state,
	int t;
	int optimalDuration;
	double optimalValue;
	int optimalStateInNextStage;
};

//The Stage class is also used in the dynamic programming class. Each Stage has a number of States and is linked to exactly one location of the respective instance.
class Stage{

public:
	Stage(int Number, int MinState, int MaxState, Node StageNode);
	void setOptimalValue(int StateNumber, double OptimalValue){states[StateNumber].setOptimalValue(OptimalValue);};
	double getOptimalValue(int StateNumber){return states[StateNumber].getOptimalValue();};

	void setOptimalStateInNextStage(int StateNumber, double OptimalStateInNextStage){states[StateNumber].setOptimalStateInNextStage(OptimalStateInNextStage);};
	double getOptimalStateInNextStage(int StateNumber){return states[StateNumber].getOptimalStateInNextStage();};

	void setOptimalDuration(int StateNumber, double OptimalDuration){states[StateNumber].setOptimalDuration(OptimalDuration);};
	double getOptimalDuration(int StateNumber){return states[StateNumber].getOptimalDuration();};
	
	
	vector<State> getStates(){return states;};
	Node getStageNode(){return stageNode;};

private:
	//Number of the stage
	int number;
	//Range of states within the stage
	int minState, maxState;
	//all states of the stage
	vector<State> states;
	//the linked location of the stage
	Node stageNode;
};

//In the DP class, the dynamic programming is taking place to solve the scheduling problem described in the paper.
class DP{

public:
	//as soon as a dp class is created, the dp is solving the scheduling problem for a given sequence "Nodes"
	DP(vector<Node> Nodes):nodes(Nodes){solve();};
	//the sequence is returned including the updated starting times and service times for each node (=location)
	vector<Node> getNodes(){return nodes;};
	//the objective value of the Nodes sequence is returned
	double getBestValue(){return bestValue;};

private:
	//this is the internal method to solve a given DP-sequence
	void solve();
	//in this function, the profit is calculated given the base profit pi, the minimum service time li and the maximum service time Li as well as the service time (Duration)
	double ServiceTimeDependentProfitFunction(double pi, int li, int Li, int Duration);
	//all nodes for of the dp sequence are saved in nodes
	vector<Node> nodes;
	//the objective value of the dp sequence is saved in bestValue
	double bestValue;
	//all stages are kept in stages
	vector<Stage> stages;

};

#endif