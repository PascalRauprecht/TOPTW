#include "hssa.h"
#include <list>
#include <iterator>	

void Sequence::setNodes(vector<Node> Nodes, string ProblemType)
{
	nodes = Nodes;
	determineSeqValue(ProblemType);
}

void Sequence::determineSeqValue(string ProblemType)
{
	//Here, the objective value of the nodes sequence is calculated and stored in seqValue
	seqValue = 0;
	if(ProblemType == "TOPTW")
	{
		for(int i = 0; i < nodes.size(); i++) seqValue += nodes[i].getBaseProfit();
	}
	else if(ProblemType == "TOPTWSTDP")
	{
		DP dynPro(nodes);
		nodes = dynPro.getNodes();
		seqValue = dynPro.getBestValue();
	}
}

double SimA::startSA(int M, int SAType, int Parameter, double T, double Alpha, int IIterB)
{
	m = M;
	saType = SAType;
	if (saType == 1)	tMax = Parameter;
	else if (saType == 2) nNonImproving = Parameter;

	t = T;
	alpha = Alpha;
	iIter = IIterB * (nodes.size() + m - 1);
	iIterB = IIterB;
	
	/*
		1. Array mit Pointer auf Nodes
		2. Array mit Startzeiten
		3. Array mit Ankunftszeiten

		4. dynamic bitset mit information, welcher Node nicht in Sequenz ist
		5. Liste mit Route

		6. insertion move
			6.1. Zufällig node welcher nicht eingefügt wurde auswählen
			6.2. Zufällige position in tour auswählen, in die node eingefügt wird

	*/
	return 1;
}

//M = number of vehicles, SAType = 1 (FHSSA), 2 (SHSSA), Parameter= tMax (FHSSA), nNonImproving (SHSSA), Variant = always 1 (potentially extandable for other modifications of the heuristic), T = temperature, OptimalValue = optional, if optimal value of MIP exists
double SimA::startSAOld(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
{
	m = M;
	saType = SAType;
	optimalValue = OptimalValue;
	if(saType == 1)	tMax = Parameter;
	else if(saType == 2) nNonImproving = Parameter;

	variant = Variant;
	t = T;
	alpha = Alpha;
	iIter = IIterB * (nodes.size() + m - 1);
	iIterB = IIterB;
	iSel = ISel;

	for(int i = 1; i <= nodes.size()-1; i++) stockedNodes.push_back(i);

	for(int i = 0; i < m; i++) 
	{
			Sequence tmpSequence(i+1);
			sequences.push_back(tmpSequence);
	}

	//the time is taken to evaluate the computational time of the heuristic
	clock_t sT, eT;
	sT = clock();

	if(variant == 1) randomSelectionMethod();

	eT = clock();
	durationOfSimAn = (double)((double)(eT-sT) / CLOCKS_PER_SEC);

	//the best found objective value is returned.
	return fBest;
}

void SimA::randomSelectionMethod()
{

	//here, the first decision is made, how much nodes are going to be included in the sequence
	int numberOfDraws = rand() % (stockedNodes.size()-m) + m + 1;

	//cout << numberOfDraws << " <- numberOfDraws" << endl;

	//-------------------Initial Sequence--------------------------
	//initial nodes are drawn randomly from the set of locations
	for(int i = 0; i < numberOfDraws; i++)
	{
		int draw = rand() % stockedNodes.size();
		int drawedNumber = stockedNodes[draw];
		/*cout << draw << "<- draw" << endl;
		cout << drawedNumber << "<- drawedNumber" << endl;
		cout << stockedNodes.size() << " <- nodeNumbers.size() " << endl;*/
		//the drawn nodes are erased from the stocked nodes
		stockedNodes.erase(stockedNodes.begin()+draw, stockedNodes.begin()+draw+1);
		/*cout << stockedNodes.size() << " <- nodeNumbers.size() " << endl;
		cout << stockedNodes[draw] << " <- nodeNumbers[draw]" << endl;
		cout << stockedNodes[draw-1] << " <- nodeNumbers[draw-1]" << endl;*/
		selectedNodes.push_back(drawedNumber);
		//system("pause");
	}

	//depending on the amount of vehicles m, m-1 zeroes are inserted (see paper)
	insertZeroesInSequence();

	//in this method, the sequence is checked, nodes which can be included in the tour are included in the dp-sequence sequences[i]
	checkVerticFeasibilityInSequence();

	//the best sequence xBest is the current initial sequence
	xBest = selectedNodes;
	fBest = 0;
	//the current objective value of all dp-sequences are computed using DP
	for(int i = 0; i < sequences.size(); i++) fBest += sequences[i].getSeqValue();

	//cout << fBest << endl;


	//-------------------actual SA Heuristic--------------------------

	int exitCriteria;

	int exitIter = 0;
	clock_t sT, eT;
	sT = clock();
	
	//depending on the heuristic type, slow or fast, the exitCriteria is initialized accordingly
	if(saType == 1) exitCriteria = tMax;
	else if(saType == 2) exitCriteria = nNonImproving;


	int loopIter = 1;
	int selIter = 1;
	int noImprovementIter = 0;
	double fLastTemperatureDecrease = fBest;
	vector<int> xIterBest = selectedNodes;
	vector<int> stockedIterBest = stockedNodes;
	double fIterBest = fBest;
	double initialT = t;

	//iIter is calculated using iIterB
	iIter = (selectedNodes.size()+m)*iIterB;


	//the heuristic loops as long as either the exitCriteria is reached or the optimalValue, if existing, from the MIP
	while(exitIter < exitCriteria && fBest < optimalValue && fIterBest < optimalValue)
	{
		vector<int> xOld = selectedNodes;
		vector<int> stockedOld = stockedNodes;


		//Swap, Insert, Reversion
		double randomNumber = rand() % 1000;
		randomNumber = randomNumber /  1000;
		
		vector<int> possiblePositions;
		for(int i = 0; i < selectedNodes.size(); i++) possiblePositions.push_back(i);
		int pos1 = possiblePositions[rand() % possiblePositions.size()];
		possiblePositions.erase(possiblePositions.begin()+pos1, possiblePositions.begin()+pos1+1);
		int pos2 = possiblePositions[rand() % possiblePositions.size()];

		//pos1 and pos2 of the sequence are either swapped, inserted, reversed or exchanged (see paper)
		if(randomNumber <= 0.3) swapSeq(pos1, pos2);
		else if(randomNumber > 0.3 && randomNumber <= 0.6) insertSeq(pos1, pos2);
		else if(randomNumber > 0.6 && randomNumber <= 0.9) reversionSeq(pos1, pos2);
		else
		{
			if(stockedNodes.size()>0 && selectedNodes[pos1] > 0)
			{
				int posInStockedNodes = 0;
				double bestValue = 0;
				int predecessor = 0;
				int successor = 0;

				//only if the node in the sequence is not zero, an exchange is taking place
				if(pos1 != 0) 
					if(selectedNodes[pos1-1] != 0)
						predecessor = selectedNodes[pos1-1];

				if(pos1 != selectedNodes.size() - 1)
					if(selectedNodes[pos1+1] != 0)
						successor = selectedNodes[pos1+1];
				
				//here, the best stocked node to be exchanged with the selected nodes is computed (see paper for more details)
				for(int i = 0; i < stockedNodes.size();i++)
				{
					double valueStockedNodes = nodes[stockedNodes[i]].getBaseProfit() / (nodes[predecessor].getDistanceToNode(nodes[stockedNodes[i]].getNumberOfNode()) + nodes[stockedNodes[i]].getMinServiceTime() + nodes[stockedNodes[i]].getDistanceToNode(nodes[successor].getNumberOfNode()));
					//cout << i << " <- i; " << valueStockedNodes << " <- VSN" << endl;
					if(valueStockedNodes > bestValue)
					{
						bestValue = valueStockedNodes;
						posInStockedNodes = i;
					}
				}

				//thats the actual exchange move
				exchangeSeq(pos1, posInStockedNodes);

			}
		}

		checkVerticFeasibilityInSequence();

		//Best Value Check
		double tmpSeqValue = 0;
		for(int i = 0; i < sequences.size(); i++) tmpSeqValue += sequences[i].getSeqValue();


		randomNumber = rand() % 1000;
		randomNumber = randomNumber /  1000;

		//if the value of current sequence is better than the best found value so far, the solution is transfered to be the best one, if it is not better, the current sequence is reset to a certain probability
		if(tmpSeqValue > fIterBest)
		{
			fIterBest = tmpSeqValue;
			xIterBest = selectedNodes;
			stockedIterBest = stockedNodes;
		}
		else if(randomNumber > exp(abs(fIterBest-tmpSeqValue)/t))
		{
			selectedNodes = xOld;
			stockedNodes = stockedOld;
		}

		//the iteration counter is increased by one
		loopIter++;

		//iIter-Check
		//if the iteration counter equals the iteration interval iIter, a local search is performed
		if(loopIter % iIter == 0)
		{
			//the temperature is decreased by alpha
			t = t * alpha;
			loopIter = 0;
			//local search

			//the current sequence is set to be the best found sequence so far
			selectedNodes = xIterBest;
			stockedNodes = stockedIterBest;
			xOld = selectedNodes;
			stockedOld = stockedNodes;

			//all possible swap moves are applied
			for(int i = 0; i < selectedNodes.size()-1; i++) 
			{
				for(int j = i+1; j < selectedNodes.size(); j++) 
				{
					swapSeq(i,j);
					checkVerticFeasibilityInSequence();

					tmpSeqValue = 0;
					for(int i = 0; i < sequences.size(); i++) tmpSeqValue += sequences[i].getSeqValue();

					if(tmpSeqValue > fIterBest)
					{
						fIterBest = tmpSeqValue;
						xIterBest = selectedNodes;
						stockedIterBest = stockedNodes;
					}

					selectedNodes = xOld;
					stockedNodes = stockedOld;
				}
			}

			xOld = xIterBest;
			stockedOld = stockedIterBest;
			selectedNodes = xIterBest;
			stockedNodes = stockedIterBest;
			//all possible insertion moves are applied
			for(int i = 0; i < selectedNodes.size()-1; i++) 
			{
				for(int j = i+1; j < selectedNodes.size(); j++) 
				{
					insertSeq(i,j);
					checkVerticFeasibilityInSequence();

					tmpSeqValue = 0;
					for(int i = 0; i < sequences.size(); i++) tmpSeqValue += sequences[i].getSeqValue();

					if(tmpSeqValue > fIterBest)
					{
						fIterBest = tmpSeqValue;
						xIterBest = selectedNodes;
						stockedIterBest = stockedNodes;
					}

					selectedNodes = xOld;
					stockedNodes = stockedOld;
				}
			}

			selectedNodes = xIterBest;
			stockedNodes = stockedIterBest;
			//iSel-Check
			//if the iteration interval iSel is reached, a random number of the currently selected nodes are exchanged with those not selected yet
			if(selIter % iSel == 0)
			{

				if(fIterBest > fBest) 
				{
					fBest = fIterBest;
					xBest = xIterBest;
					stockedBest = stockedIterBest;
				}

				vector<int> nodesFromSelectedNodes;
				//Move of Nodes into stocked Nodes
				
				//only if the current sequence is larger than 1, nodes are potentially moved to the stocked nodes
				if(selectedNodes.size() > 1)
				{
					randomNumber = rand() % (selectedNodes.size()-1);

					for(int i = 0; i < randomNumber; i++)
					{
						int tmpPosInSelectedNodes = rand() % selectedNodes.size();
						int tmpNumber = selectedNodes[tmpPosInSelectedNodes];
						if(tmpNumber != 0)
						{
							selectedNodes.erase(selectedNodes.begin()+tmpPosInSelectedNodes);
							nodesFromSelectedNodes.push_back(tmpNumber);
						}
					}
				}
				
				//only if stocked nodes exist, a random number of them is chosen to be inserted in the current sequence
				if(stockedNodes.size() > 0)
				{
					randomNumber = rand() % (stockedNodes.size()+1);

					//cout << randomNumber << " <- Random Number" << endl;
					for(int i = 0; i < randomNumber; i++)
					{
						int tmpPosInStockedNodes = rand() % stockedNodes.size();
						int tmpNumber = stockedNodes[tmpPosInStockedNodes];
					
						stockedNodes.erase(stockedNodes.begin()+tmpPosInStockedNodes);

						int tmpPosInSelectedNodes = 0;
						if(selectedNodes.size() > 0) tmpPosInSelectedNodes = rand() % selectedNodes.size();

						selectedNodes.insert(selectedNodes.begin()+tmpPosInSelectedNodes,tmpNumber);
					
					}
				}

				//if selected nodes were chosen to be deleted, they are inserted in the stocked nodes
				if(nodesFromSelectedNodes.size() > 0) stockedNodes.insert(stockedNodes.begin(),nodesFromSelectedNodes.begin(), nodesFromSelectedNodes.end());

				checkVerticFeasibilityInSequence();

				tmpSeqValue = 0;
				for(int i = 0; i < sequences.size(); i++) tmpSeqValue += sequences[i].getSeqValue();

				if(tmpSeqValue > fBest)
				{
					fBest = tmpSeqValue;
					xBest = selectedNodes;
					stockedBest = stockedNodes;
				}

				xIterBest = selectedNodes;
				stockedIterBest = stockedNodes;
				fIterBest = tmpSeqValue;

				/*cout << endl;
				cout << "After: " << endl;
				cout << "Selected Nodes: " << endl;
				copy (selectedNodes.begin(), selectedNodes.end(),
					  ostream_iterator<int>(cout," "));
				cout << endl;

				cout << "Stocked Nodes: " << endl;
				copy (stockedNodes.begin(), stockedNodes.end(),
					  ostream_iterator<int>(cout," "));
				cout << endl;
				cout << endl;
				cout << fBest << " <- FBest " << endl;
				cout << endl;*/
				//cout << endl;
				//system("pause");

				selIter = 0;
				//the temperature is reset to be the initial temperature
				t = initialT;
				iIter = (selectedNodes.size()+m)*iIterB;

				//if the objective value has not changed since the last iteration, noImprovementIter is increased by one
				if(fBest  == fLastTemperatureDecrease) noImprovementIter++;
				else 
				{
					//cout << fBest << endl;
					noImprovementIter = 0;
				}

				fLastTemperatureDecrease = fBest;
			}

			
			selIter++;
		}
		
		//depending on the heuristic variant (fast or slow), the exitIterator is updated
		if(saType == 2) exitIter = noImprovementIter;
		else if(saType == 1) 
		{
			eT = clock(); 
			exitIter = (double)((double)(eT-sT) / CLOCKS_PER_SEC);
		}
	}

	//If the objective value of the last iteration interval fIterBest is bigger than fBest, fBest is taking those values
	if(fIterBest > fBest) 
	{
		fBest = fIterBest;
		xBest = xIterBest;
		stockedBest = stockedIterBest;
	}

	selectedNodes = xBest;
	stockedNodes = stockedBest;
	checkVerticFeasibilityInSequence();
}

void SimA::insertZeroesInSequence()
{
	for(int i = 0; i < m-1; i++)
	{
		int positionInSelectedNodes = rand() % selectedNodes.size();
		selectedNodes.insert(selectedNodes.begin()+positionInSelectedNodes,0);
	}
}

void SimA::swapSeq(int pos1, int pos2)
{
	swap(selectedNodes[pos1],selectedNodes[pos2]);
}

void SimA::insertSeq(int pos1, int pos2)
{
	int tmpValue = selectedNodes[pos1];
	selectedNodes.erase(selectedNodes.begin()+pos1,selectedNodes.begin()+pos1+1);
	if(pos1 < pos2) pos2 = pos2 - 1;
	selectedNodes.insert(selectedNodes.begin()+pos2,tmpValue);

	//cout << pos1 << " <- pos1, " << pos2 << " <- pos2" << endl;
	//cout << "After: " << endl;
	//copy (selectedNodes.begin(), selectedNodes.end(),
 //         ostream_iterator<int>(cout," "));
 //   cout << endl;
}

void SimA::reversionSeq(int pos1, int pos2)
{

	int tmpPos = pos1;
	if(pos1 > pos2) 
	{
		pos1 = pos2;
		pos2 = tmpPos;
	}
	reverse(selectedNodes.begin()+pos1,selectedNodes.begin()+pos2+1);

}

void SimA::exchangeSeq(int posSel, int posStock)
{

	int tmpValue = selectedNodes[posSel];
	selectedNodes.insert(selectedNodes.begin()+posSel,stockedNodes[posStock]);
	stockedNodes.push_back(tmpValue);
	stockedNodes.erase(stockedNodes.begin()+posStock);
	selectedNodes.erase(selectedNodes.begin()+posSel+1);

}

void SimA::checkVerticFeasibilityInSequence()
{
	int iterSequence = 0;
	int travelTime = 0;
	vector<Node> tmpNodes;

	int tmpNodeIter = 0;
	tmpNodes.push_back(nodes[0]);
	//in this loop, a check is performed whether the next node can still be inserted while not hurting the depot closing time window condition
	for(int i = 0; i < selectedNodes.size(); i++)
	{
		
		if(selectedNodes[i] == 0)
		{
			tmpNodes.push_back(nodes[0]);
			sequences[iterSequence].setNodes(tmpNodes, inst.getProblemType());
			travelTime = 0;
			iterSequence++;
			tmpNodes.clear();
			tmpNodeIter = 0;
			tmpNodes.push_back(nodes[0]);
		}
		else
		{
			//cout << "SN: " << selectedNodes[i] << "; TT: " << travelTime << "; DTN: " << tmpNodes[tmpNodeIter].getDistanceToNode(selectedNodes[i]) << "; MinServiceTime: " << nodes[selectedNodes[i]].getMinServiceTime() << "; DTN0: " << nodes[selectedNodes[i]].getDistanceToNode(0) << "; N0DD: " << nodes[0].getDueDate() << "; SNDT: " << nodes[selectedNodes[i]].getDueDate() << "; RT: " << nodes[selectedNodes[i]].getReadyTime() << endl;

			int distanceFromLastFeasibleNodeToCurrentNode = tmpNodes[tmpNodeIter].getDistanceToNode(selectedNodes[i]);
			int minServiceTimeCurrentNode = nodes[selectedNodes[i]].getMinServiceTime();
			int distanceFromCurrentNodeToDepot = nodes[selectedNodes[i]].getDistanceToNode(0);
			int dueDateOfDepot = nodes[0].getDueDate();
			int dueDateOfCurrentNode = nodes[selectedNodes[i]].getDueDate();
			int readyTimeOfCurrentNode = nodes[selectedNodes[i]].getReadyTime();


			if(travelTime + distanceFromLastFeasibleNodeToCurrentNode + minServiceTimeCurrentNode + 
				distanceFromCurrentNodeToDepot <= dueDateOfDepot
				&& travelTime + distanceFromLastFeasibleNodeToCurrentNode <= dueDateOfCurrentNode
				&& readyTimeOfCurrentNode + minServiceTimeCurrentNode + distanceFromCurrentNodeToDepot <= dueDateOfDepot)
			{
				/*cout << travelTime << " <- TravelTime" << endl;
				cout << distanceFromLastFeasibleNodeToCurrentNode << " <- Distance from Last to Current" << endl;
				cout << minServiceTimeCurrentNode << " <- minService Time current node" << endl;*/
				if(travelTime + distanceFromLastFeasibleNodeToCurrentNode >= readyTimeOfCurrentNode)
				travelTime += distanceFromLastFeasibleNodeToCurrentNode + minServiceTimeCurrentNode;
				else travelTime = readyTimeOfCurrentNode + minServiceTimeCurrentNode;
				//cout << travelTime << " <- new TravelTime" << endl;
				//system("pause");
				tmpNodes.push_back(nodes[selectedNodes[i]]);
				tmpNodeIter++;
			}
			
		}
	}

	tmpNodes.push_back(nodes[0]);
	sequences[iterSequence].setNodes(tmpNodes, inst.getProblemType());
}