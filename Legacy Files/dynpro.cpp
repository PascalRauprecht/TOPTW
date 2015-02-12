#include "dynpro.h"

Stage::Stage(int Number, int MinState, int MaxState, Node StageNode):number(Number), minState(MinState), maxState(MaxState), stageNode(StageNode)
{
	//initialization of all states of the stage
	for(int i = minState; i <= maxState;i++)
	{
		State tmpState(i);
		states.push_back(tmpState);
	}
}


double DP::ServiceTimeDependentProfitFunction(double pi, int li, int Li, int Duration)
{
	int t = Duration - li;
	//this is the objective function of the TOPTWSTDP, which is also illustrated and explained in the paper
	double fctValue = li * pi + 0.5 * (pi * t + (pi * pow((double) t, 2)) / (Li - li));

	return fctValue;
}

void DP::solve()
{
	//only if the sequence is not only including the depot, an evaluation is taking part. otherwise the objective value is set to 0
	if(nodes.size() > 2)
	{
		vector<int> minStates;
		vector<int> maxStates;


		//The process of the evaluation of the soonest Arriving time SAT and latest arriving time for each stage is also explained in the paper
		for(int i = 0; i < nodes.size(); i++)
		{
			int soonestArrivingTime = 0;
			int latestArrivingTime = 0;
		
			for(int j = 0; j < i; j++)
			{
				soonestArrivingTime = soonestArrivingTime + nodes[j].getMinServiceTime() + nodes[j].getDistanceToNode(nodes[j+1].getNumberOfNode());
				if(soonestArrivingTime < nodes[j+1].getReadyTime()) soonestArrivingTime = nodes[j+1].getReadyTime();

				latestArrivingTime = latestArrivingTime + nodes[j].getMaxServiceTime() + nodes[j].getDistanceToNode(nodes[j+1].getNumberOfNode());
				if(latestArrivingTime < nodes[j+1].getReadyTime()) latestArrivingTime = nodes[j+1].getReadyTime();
				if(latestArrivingTime > nodes[j+1].getDueDate()) latestArrivingTime = nodes[j+1].getDueDate();
			}
			//latestArrivingTime = nodes[i].getDueDate();
			maxStates.push_back(latestArrivingTime);
			minStates.push_back(soonestArrivingTime);
		}

		int latestArrivingTime = 0;
		for(int i = nodes.size()-1; i >=0; i--)
		{
			if(i == nodes.size()-1) latestArrivingTime = nodes[i].getDueDate();
			else
			{
				latestArrivingTime = latestArrivingTime - nodes[i].getDistanceToNode(nodes[i+1].getNumberOfNode())-nodes[i].getMinServiceTime();
				if(latestArrivingTime > nodes[i].getDueDate()) latestArrivingTime = nodes[i].getDueDate();

				if(latestArrivingTime < maxStates[i]) maxStates[i] = latestArrivingTime;
			}

		}

		/*cout << "SAT:  \t";
		for(int i = 0; i < minStates.size(); i++) cout  << i << ": " << minStates[i] << "  \t";
		cout << endl;
		cout << "LAT:  \t";
		for(int i = 0; i < maxStates.size(); i++) cout  << i << ": " << maxStates[i] << "  \t";
		cout << endl << endl;
		for(int i = 0; i < nodes.size(); i++) cout << nodes[i].getNumberOfNode() << " <- Node Number" << endl;*/
		//system("pause");

		//all neccassary stages are initialized
		for(int i = 0; i < nodes.size();i++)
		{
			Stage tmpStage(i+1,minStates[i],maxStates[i], nodes[i]);
			stages.push_back(tmpStage);
		}

		//in this loop, the minimum and maximum stay time for every stage is evaluated
		for(int i = stages.size()-1; i >= 0; i--)
		{
			if(i == stages.size()-1) 
			{
				for(int j = 0; j < stages[i].getStates().size();j++) 
				{
					stages[i].setOptimalValue(j,0);
					stages[i].setOptimalDuration(j,0);
				}
			}
			else if(i == 0)
			{
				for(int j = 0; j < stages[i].getStates().size();j++) stages[i].setOptimalDuration(j,0);
			}
			else 
			{
				for(int j = 0; j < stages[i].getStates().size(); j++)
				{
					int minDurationOfStayInNodeI;
					int minDurationOfStayInNodeIStateInIPlusOne;
					int maxDurationOfStayInNodeI;
					int maxDurationOfStayInNodeIStateInIPlusOne;

					for(int k = 0; k < stages[i+1].getStates().size();k++)
					{
						minDurationOfStayInNodeI = stages[i+1].getStates()[k].getT()-
													stages[i].getStageNode().getDistanceToNode(stages[i+1].getStageNode().getNumberOfNode())- 
													stages[i].getStates()[j].getT();

						minDurationOfStayInNodeIStateInIPlusOne = k;

						if(minDurationOfStayInNodeI > stages[i].getStageNode().getMaxServiceTime()) 
						{
							minDurationOfStayInNodeI = stages[i].getStageNode().getMaxServiceTime();
							break;
						}
						
						if(minDurationOfStayInNodeI >= stages[i].getStageNode().getMinServiceTime()) break;
					}




					for(int k = stages[i+1].getStates().size()-1; k >=0;k--)
					{
						maxDurationOfStayInNodeI = stages[i+1].getStates()[k].getT()-
													stages[i].getStageNode().getDistanceToNode(stages[i+1].getStageNode().getNumberOfNode())- 
													stages[i].getStates()[j].getT();

						maxDurationOfStayInNodeIStateInIPlusOne = k;

						if(maxDurationOfStayInNodeI >= stages[i].getStageNode().getMinServiceTime()
							&& maxDurationOfStayInNodeI <= stages[i].getStageNode().getMaxServiceTime()) break;
					}

					/*cout << "Stage " << i << " / State " << stages[i].getStates()[j].getT() << ": " << minDurationOfStayInNodeI << " and " << stages[i+1].getStates()[minDurationOfStayInNodeIStateInIPlusOne].getT() << " - " << maxDurationOfStayInNodeI << " and " << stages[i+1].getStates()[maxDurationOfStayInNodeIStateInIPlusOne].getT()  << endl;

					cout << minDurationOfStayInNodeIStateInIPlusOne << " <- minDurationOfStayInNodeIStateInIPlusOne; " << maxDurationOfStayInNodeIStateInIPlusOne << " <- maxDurationOfStayInNodeIStateInIPlusOne" << endl;
					cout << stages[i].getStates()[j].getT() << " <- State in Stage " << i << "; " << stages[i+1].getStates()[minDurationOfStayInNodeIStateInIPlusOne].getT() << " <- MinState; " << stages[i+1].getStates()[maxDurationOfStayInNodeIStateInIPlusOne].getT() << " <- MaxState; " << stages[i].getStageNode().getDistanceToNode(stages[i+1].getStageNode().getNumberOfNode()) << " Distance to Node; " << stages[i+1].getStageNode().getDueDate() << " <- Due Date Node + 1" << endl;*/

					double optimalValue = 0;
					int optimalStateInIPlusOne = 0;
					int optimalDuration = 0;

					//In this loop, the optimal value and duration for each stage are computed. Every state of each stage is looped through.

					for(int k = minDurationOfStayInNodeIStateInIPlusOne; k <= maxDurationOfStayInNodeIStateInIPlusOne;k++)
					{
						int durationOfStayInNodeI = stages[i+1].getStates()[k].getT()-
													stages[i].getStageNode().getDistanceToNode(stages[i+1].getStageNode().getNumberOfNode())- 
													stages[i].getStates()[j].getT();
						if(durationOfStayInNodeI > stages[i].getStageNode().getMaxServiceTime()) 
						durationOfStayInNodeI = stages[i].getStageNode().getMaxServiceTime();

						double tmpValue = ServiceTimeDependentProfitFunction(stages[i].getStageNode().getBaseProfit(), stages[i].getStageNode().getMinServiceTime(),stages[i].getStageNode().getMaxServiceTime(),durationOfStayInNodeI) + stages[i+1].getOptimalValue(k);	

						if(tmpValue >= optimalValue)
						{
							optimalValue = tmpValue;
							optimalStateInIPlusOne = k;
							optimalDuration = durationOfStayInNodeI;
						}
						
					}

					stages[i].setOptimalValue(j,optimalValue);
					stages[i].setOptimalStateInNextStage(j,optimalStateInIPlusOne);
					stages[i].setOptimalDuration(j,optimalDuration);

				}
			}
		}

		//the optimal value and optimal state in the next stage for stage 0 are set individually
		stages[0].setOptimalValue(0,stages[1].getOptimalValue(0));
		stages[0].setOptimalStateInNextStage(0,0);

		//the bestValue of the first stage equals the objective value of the whole dp-sequence
		bestValue = stages[1].getOptimalValue(0);

		//starting time of location 1 is set to 0, as well as the duration of the depot
		nodes[0].setStartingTime(0);
		nodes[0].setDuration(0);
		nodes[nodes.size()-1].setDuration(0);

		int optimalStateInThisStage = 0;
		int optimalStateInNextStage = 0;
		//the beforehand calculated starting times and durations are forwarded to the nodes (=locations)
		for(int i = 1; i < nodes.size(); i++)
		{
			nodes[i].setStartingTime(stages[i].getStates()[stages[i-1].getOptimalStateInNextStage(optimalStateInThisStage)].getT());
			nodes[i].setDuration(stages[i].getOptimalDuration(stages[i-1].getOptimalStateInNextStage(optimalStateInThisStage)));
		}

	}
	else 
	{
		bestValue = 0;
		nodes[0].setStartingTime(0);
		nodes[0].setDuration(0);
		nodes[nodes.size()-1].setDuration(0);
		nodes[nodes.size()-1].setStartingTime(0);
	}

	
	//system("pause");
}