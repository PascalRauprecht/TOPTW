#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <math.h>

using namespace std;


//The class Node is used to save all required information for each location (=Node). The variable names should be selfexplaining. The vector<int> distanceToNode gives the information about the distance to all other nodes in the specific instance. 

class Node{

public:
	Node(int NumberOfNode, int XCoord, int YCoord, int ReadyTime, int Demand, int DueDate, int MinServiceTime, int MaxServiceTime, double BaseProfit, string Name):
		numberOfNode(NumberOfNode), xCoord(XCoord), yCoord(YCoord), readyTime(ReadyTime), demand(Demand), dueDate(DueDate), minServiceTime(MinServiceTime), maxServiceTime(MaxServiceTime), baseProfit(BaseProfit), name(Name){};


	//getter and setter methods
	int getNumberOfNode(){return numberOfNode;};
	int getXCoord(){return xCoord;};
	int getYCoord(){return yCoord;};
	int getReadyTime(){return readyTime;};
	int getDemand(){return demand;};
	int getDueDate(){return dueDate;};
	int getMinServiceTime(){return minServiceTime;};
	int getMaxServiceTime(){return maxServiceTime;};
	double getBaseProfit(){return baseProfit;};
	void setBaseProfit(double BaseProfit){baseProfit = BaseProfit;};
	string getName(){return name;};

	int getStartingTime(){return startingTime;};
	void setStartingTime(int StartingTime){startingTime = StartingTime;};

	int getDuration(){return duration;};
	void setDuration(int Duration){duration = Duration;};

	int getDistanceToNode(int i){return distanceToNode[i];}; 
	void setDistanceToNode(vector<int> DistanceToNode){distanceToNode = DistanceToNode;};

	int getDelay(int FromNode, int ToNode){ return delayFromNodeToNode[FromNode][ToNode]; };
	void setDelay(vector<vector<int>> DelayFromNodeToNode){ delayFromNodeToNode = DelayFromNodeToNode; };

	int isFeasible(int FromNode, int ToNode){ return feasibleFromNodeToNode[FromNode][ToNode]; };
	void setFeasibility(vector<vector<bool>> FeasibleFromNodeToNode){ feasibleFromNodeToNode = FeasibleFromNodeToNode; };


private:
	int numberOfNode, xCoord, yCoord, readyTime, demand, dueDate, minServiceTime, maxServiceTime;
	double baseProfit;
	string name;
	vector<int> distanceToNode;
	vector<vector<int>> delayFromNodeToNode;
	vector<vector<bool>> feasibleFromNodeToNode;
	int startingTime, duration;
	
};


//The Instance class is storing all information gathered from the specific instance txt file from Vansteenwegen (see paper)
class Instance{

public:
	//constructor of class
	Instance(string Path, string ProblemType, int DecimalPlaces, int MinServiceTime, int MaxServiceTime);
	//a specific node at NodePosition is returned
	Node getNode(int NodePosition){return nodes[NodePosition];};
	//the nodes from the instance are returned here
	vector<Node> getNodes(){return nodes;};
	//the name of the instance is returned here
	string getName(){return name;};
	int getInstanceSize(){return nodes.size();};
	double getValue(int M, int SAType, int Parameter, double T, double Alpha, int IIter, int ISel, int Variant);
	string getProblemType(){return problemType;};

private:
	//all locations are saved in this class
	vector<Node> nodes;
	string name;
	string problemType;
	int numberOfSequences;
	int numberOfNodes;
};

#endif