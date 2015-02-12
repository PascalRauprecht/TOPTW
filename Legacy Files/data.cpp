#include "data.h"

Instance::Instance(string Path, string ProblemType, int DecimalPlaces, int MinServiceTime, int MaxServiceTime){

	//TOPTW-Problem
	problemType = ProblemType;

	string fileName = Path.substr( Path.find_last_of("\/") + 1 );
	fileName = fileName.substr(0, fileName.find_first_of("\."));
	name = fileName;

	string filename = Path;
	string line;
	double numberOfNode, xCoord, yCoord, readyTime, demand, dueDate, baseProfit;

	int decimalPlaces = pow((double) 10,DecimalPlaces);


	//here the txt file is read and transmitted to Node objects
	ifstream infile(filename.c_str());
	int tmp;

	getline(infile, line);
	stringstream(line) >> tmp >> numberOfSequences >> numberOfNodes;

	for(int i = 0; i < 2;i++) getline(infile, line);
	
	
	//this is the first line, which belongs to the depot
	stringstream(line) >> numberOfNode >> xCoord >> yCoord >> demand >> baseProfit >> tmp >> tmp >> readyTime >> dueDate;

	Node depotNode(numberOfNode, xCoord, yCoord, readyTime*decimalPlaces, demand*decimalPlaces, dueDate*decimalPlaces, 0, 0, 0, "Depot");
	nodes.push_back(depotNode);

	for(int i = 0; i < numberOfNodes;i++)
	{
		getline(infile, line);
		stringstream(line) >> numberOfNode >> xCoord >> yCoord >> demand >> baseProfit >> tmp >> tmp >> tmp >> readyTime >> dueDate;
		numberOfNode = i + 1;

		if(problemType == "TOPTW")
		{
			Node tmpNode(numberOfNode, xCoord, yCoord, readyTime*decimalPlaces, demand*decimalPlaces, dueDate*decimalPlaces, demand*decimalPlaces, demand*decimalPlaces, baseProfit, "Node" + (i + 1));
			nodes.push_back(tmpNode);
		}
		else if (problemType == "TOPTWSTDP")
		{
			Node tmpNode(numberOfNode, xCoord, yCoord, readyTime*decimalPlaces, demand*decimalPlaces, dueDate*decimalPlaces, MinServiceTime*decimalPlaces, MaxServiceTime*decimalPlaces, baseProfit, "Node" + (i + 1));
			nodes.push_back(tmpNode);
		}

		
	}

	infile.close();

	//in this loop, the euclidean distances are calculated
	for(int i = 0; i < nodes.size();i++)
	{
		vector<int> distanceToNode;
		for(int j = 0; j < nodes.size();j++)
		{
			if(i == j) distanceToNode.push_back(0);
			else 
			{
				int xI = nodes[i].getXCoord();
				int yI = nodes[i].getYCoord();
				int xJ = nodes[j].getXCoord();
				int yJ = nodes[j].getYCoord();
				
				int difX = abs(xI-xJ);
				int difY = abs(yI-yJ);
				
				int arcIJ = floor(sqrt(pow((double)difX,2) + pow((double)difY,2))*decimalPlaces);
				distanceToNode.push_back(arcIJ);
			}
		}
		nodes[i].setDistanceToNode(distanceToNode);
	}

	//calculation of delays for every scenario
	//1. Distance from predecessor to successor
	//2. New distance via new node
	for (int i = 0; i < nodes.size(); i++)
	{
		vector<vector<int>> delayFromNodeJToNodeK;
		for (int j = 0; j < nodes.size(); j++)
		{
			vector<int> tmpDelayFromJ;
			for (int k = 0; k < nodes.size(); k++)
			{
				int distanceFromJToK = nodes[j].getDistanceToNode(k);
				int distanceFromJToI = nodes[j].getDistanceToNode(i);
				int distanceFromIToK = nodes[i].getDistanceToNode(k);
				int alternativeTime = distanceFromJToI + distanceFromIToK + nodes[i].getMinServiceTime();
				int delay = alternativeTime - distanceFromJToK;
				tmpDelayFromJ.push_back(delay);

				if (delay > nodes[k].getDueDate() - nodes[k].getReadyTime())
				{

				}
				else
				{
				
				}
			}
			delayFromNodeJToNodeK.push_back(tmpDelayFromJ);
		}
		nodes[i].setDelay(delayFromNodeJToNodeK);
	}

	//if delay is bigger than closing time window - opening time window -> infeasible



}
