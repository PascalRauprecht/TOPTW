#include "mip.h"

MIP::MIP(Instance inst)
{
	//aditionally to the existing nodes in the instance, the depot is added a second time as last node
	nodes = inst.getNodes();
	nodes.push_back(nodes[0]);
}

double MIP::solve(int M)
{
	m = M;
	double value = 0;
	IloEnv env;
	try 
	{
		IloModel model(env);

		typedef IloArray<IloNumVarArray> NumVarMatrix;
		typedef IloArray<NumVarMatrix>   NumVar3Matrix;
		typedef IloArray<NumVar3Matrix>   NumVar4Matrix;

		IloInt numberOfClients = nodes.size();
		int bigM = 10000;

		
		//here, all decision variables are initialized
		vector<vector <int> > tForEachNode;
		for(int i = 0; i < nodes.size(); i++)
		{
			vector<int> tmpVector;
			for(int j = 0; j <= nodes[i].getMaxServiceTime() - nodes[i].getMinServiceTime(); j++) tmpVector.push_back(j);
			tForEachNode.push_back(tmpVector);
			//cout << tForEachNode[i].size() << " Size TForEachNode" << endl;
			//system("pause");
		}

		NumVar3Matrix leavingBool(env, numberOfClients);
		 /* initialize this matrix */
		for(int i=0; i < numberOfClients; i++) 
		{

			leavingBool[i] = NumVarMatrix(env, tForEachNode[i].size());
			for(int o = 0; o < tForEachNode[i].size(); o++)
			{
				leavingBool[i][o] = IloNumVarArray(env, m);
				for(int d=0; d < m; d++) 
				{
					leavingBool[i][o][d] = IloNumVar(env, 0.0, 1.0, ILOBOOL);
				}
			}
		}

		NumVarMatrix durationOfStay(env, numberOfClients);
		 /* initialize this matrix */
		for(int i=0; i < numberOfClients; i++) 
		{
			durationOfStay[i] = IloNumVarArray(env, m);
			for(int d=0; d < m; d++) 
			{
				durationOfStay[i][d] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
			}
		}

		NumVarMatrix startingTime(env, numberOfClients);
		 /* initialize this matrix */
		for(int i=0; i < numberOfClients; i++) 
		{
			startingTime[i] = IloNumVarArray(env, m);
			for(int d=0; d < m; d++) 
			{
				startingTime[i][d] = IloNumVar(env, 0.0, IloInfinity, ILOINT);
			}
		}

		NumVarMatrix isVisited(env, numberOfClients);
		 /* initialize this matrix */
		for(int i=0; i < numberOfClients; i++) 
		{
			isVisited[i] = IloNumVarArray(env, m);
			for(int d=0; d < m; d++) 
			{
				isVisited[i][d] = IloNumVar(env, 0.0, 1.0, ILOBOOL);
			}
		}

		NumVar3Matrix arcUsed(env, numberOfClients);
		 /* initialize this matrix */
		for(int i=0; i < numberOfClients; i++) 
		{

			arcUsed[i] = NumVarMatrix(env, numberOfClients);
			for(int j = 0; j < numberOfClients; j++)
			{
				arcUsed[i][j] = IloNumVarArray(env, m);
				for(int d=0; d < m; d++) 
				{
					arcUsed[i][j][d] = IloNumVar(env, 0.0, 1.0, ILOBOOL);
				}
			}
		}

		//IloNumVarArray var(env);
		//for(int i = 0; i < NbProducts;i++)
		//{
		//	var.add(IloNumVar(env, 0.0, 1.0, ILOINT));
		//};


		//this is the objective function facilitated using cplex and concert technology (see paper for more details about the objective function of the mathematical model of the TOPTWSTDP)
		IloExpr exprOF(env);

		for(int i = 1; i < numberOfClients-1;i++)
		{
			for(int d = 0; d < m; d++)
			{
				exprOF += nodes[i].getMinServiceTime() * nodes[i].getBaseProfit() * isVisited[i][d];

				IloExpr exprLeavingTime(env);
				for(int o = 0; o < tForEachNode[i].size(); o++)
				{
					exprLeavingTime += tForEachNode[i][o] * leavingBool[i][o][d];
				}

				exprOF += 0.5 * nodes[i].getBaseProfit() * exprLeavingTime;
				exprOF += (nodes[i].getBaseProfit() * exprLeavingTime * exprLeavingTime) / (2 * (nodes[i].getMaxServiceTime() - nodes[i].getMinServiceTime()));
				
				//exprOF += ((nodes[i].getBaseProfit() * durationOfStay[i][d]) - ((nodes[i].getBaseProfit() * (durationOfStay[i][d] - nodes[i].getMinServiceTime())) / 2) + ((nodes[i].getBaseProfit() * (durationOfStay[i][d] - nodes[i].getMinServiceTime()) * (durationOfStay[i][d] - nodes[i].getMinServiceTime())) / (2 * (nodes[i].getMaxServiceTime()-nodes[i].getMinServiceTime())))) ;//* isVisited[i][d];
			}
		};

		//the objective function is added to the model
		model.add(IloMaximize(env, exprOF));
		
		
		IloRangeArray minTimeCon(env);
		IloRangeArray maxTimeCon(env);
		IloRangeArray maxNotVisitedTimeCon(env);
		IloRangeArray leaveOnceCon(env);

		//in this constraint, the connection between the leavingBool decision variable and the duration of stay is made
		for(int i = 0; i < numberOfClients;i++)
		{
			for(int d = 0; d < m; d++)
			{
				IloExpr exprMinMaxTimeCon(env);
				IloExpr exprLeaveOnceCon(env);
				for(int o = 0; o < tForEachNode[i].size(); o++)
				{
					exprLeaveOnceCon += leavingBool[i][o][d];
					exprMinMaxTimeCon += tForEachNode[i][o] * leavingBool[i][o][d];
				}

				IloExpr exprDurationTimeCon(env);
				exprDurationTimeCon += durationOfStay[i][d] - nodes[i].getMinServiceTime() - exprMinMaxTimeCon;

				minTimeCon.add(exprDurationTimeCon  - isVisited[i][d] * bigM >= (-1 * bigM));
				maxTimeCon.add(exprDurationTimeCon  + isVisited[i][d] * bigM <= bigM);
				maxNotVisitedTimeCon.add((durationOfStay[i][d] - isVisited[i][d] * nodes[i].getMaxServiceTime()) <= 0);
				leaveOnceCon.add(exprLeaveOnceCon - isVisited[i][d] == 0);
			}
		}

		model.add(minTimeCon);
		model.add(maxTimeCon);
		model.add(maxNotVisitedTimeCon);
		model.add(leaveOnceCon);


		//in this constraint, it is stated that the depot is departed and arrived as often as they are vehicles in the problem
		IloRangeArray outCon(env);
		IloRangeArray inCon(env);

		IloExpr exprOutCon(env);
		IloExpr exprInCon(env);

		for(int i = 0; i < numberOfClients;i++)
		{
			for(int d = 0; d < m; d++)
			{
				exprOutCon += arcUsed[0][i][d];
				exprInCon += arcUsed[i][numberOfClients-1][d];
			}
		}

		outCon.add(exprOutCon == m);
		inCon.add(exprInCon == m);

		model.add(outCon);
		model.add(inCon);

		//this constraint says, that every location can be visited at most once and has to be visited and left to and from one location if it is visited
		IloRangeArray outVarCon(env);
		IloRangeArray inVarCon(env);
		for(int k = 1; k < numberOfClients-1;k++)
		{
			for(int d = 0; d < m; d++)
			{
				IloExpr exprOutVarCon(env);
				IloExpr exprInVarCon(env);
				for(int i = 0; i < numberOfClients - 1; i++)
				{
					exprOutVarCon += arcUsed[i][k][d];
				}
				for(int j = 1; j < numberOfClients; j++)
				{
					exprInVarCon += arcUsed[k][j][d];
				}

				outVarCon.add(exprOutVarCon - isVisited[k][d] == 0);
				inVarCon.add(exprInVarCon - isVisited[k][d] == 0);

			}
		}

		model.add(outVarCon);
		model.add(inVarCon);

		//this constraints is forbidding a tour from location k to location k
		IloRangeArray noCirclesCon(env);
		for(int k = 0; k < numberOfClients;k++)
		{
			for(int d = 0; d < m; d++)
			{
				noCirclesCon.add(arcUsed[k][k][d] == 0);
			}
		}
		model.add(noCirclesCon);

		//the leavingBool can in sum only be 0 or 1, if a location is visited the sum is 1, if not the sum is 0
		IloRangeArray boolLeavingTimeCon(env);
		for(int i = 0; i < numberOfClients;i++)
		{
			for(int d = 0; d < m; d++)
			{
				IloExpr exprBoolLeavingTimeCon(env);
				for(int o = 0; o < tForEachNode[i].size(); o++)
				{
					exprBoolLeavingTimeCon += leavingBool[i][o][d];
				}
				boolLeavingTimeCon.add(exprBoolLeavingTimeCon - isVisited[i][d] == 0);
			}
		}

		model.add(boolLeavingTimeCon);

		//this constraint is ensuring the logical time connection between two cities (see paper)
		IloRangeArray timeCon(env);
		for(int i = 0; i < numberOfClients;i++)
		{
			for(int j = 0; j < numberOfClients;j++)
			{
				for(int d = 0; d < m; d++)
				{
					IloExpr exprTimeCon(env);
					int numberOfJ;
					if(j < numberOfClients - 1) numberOfJ = j;
					else numberOfJ = 0;
					exprTimeCon = startingTime[i][d] + durationOfStay[i][d] + nodes[i].getDistanceToNode(numberOfJ) - startingTime[j][d] + bigM * arcUsed[i][j][d];
					timeCon.add(exprTimeCon <= bigM);

				}
			}
		}

		model.add(timeCon);


		//this constraint sets the starting time of the depot to 0
		IloRangeArray depotStartingTimeCon(env);
		for(int d = 0; d < m; d++)
		{
			depotStartingTimeCon.add(startingTime[0][d] == 0);
		}
		model.add(depotStartingTimeCon);

		//the depot is set to be visited by all locations when leaving and arriving
		IloRangeArray depotIsVisitedCon(env);
		for(int d = 0; d < m; d++)
		{
			depotIsVisitedCon.add(isVisited[0][d] == 1);
			depotIsVisitedCon.add(isVisited[numberOfClients-1][d] == 1);
		}
		model.add(depotIsVisitedCon);

		//a location can be at maximum only be visited by one vehicle
		IloRangeArray visitRestCon(env);
		for(int k = 1; k < numberOfClients-1;k++)
		{
			IloExpr exprVisitRestCon(env);
			for(int d = 0; d < m; d++)
			{
				exprVisitRestCon += isVisited[k][d];
			}
			visitRestCon.add(exprVisitRestCon <= 1);
		}

		model.add(visitRestCon);

		//constraints which ensures the depot closing time windows constraint is satisfied (see paper)
		IloRangeArray maxTotalTimeCon(env);
		for(int d = 0; d < m; d++)
		{
			IloExpr exprMaxTotalTimeCon(env);
			for(int i = 0; i < numberOfClients - 1; i++)
			{
				IloExpr exprMaxTotalTimeSTCon(env);
				for(int j = 1; j < numberOfClients; j++)
				{
					int numberOfNodeJ;
					if(j < numberOfClients - 1) numberOfNodeJ = j;
					else numberOfNodeJ = 0;

					exprMaxTotalTimeSTCon += nodes[i].getDistanceToNode(numberOfNodeJ) * arcUsed[i][j][d];
				}
				exprMaxTotalTimeCon += durationOfStay[i][d] + exprMaxTotalTimeSTCon;
			}
			maxTotalTimeCon.add(exprMaxTotalTimeCon <= nodes[0].getDueDate());
		}

		model.add(maxTotalTimeCon);

		//this constraint ensures, that the starting time is only between the opening time window of a location and the closing time window of the location
		IloRangeArray readyTimeCon(env);
		IloRangeArray dueDateCon(env);
		IloRangeArray noVisitStartingTimeCon(env);

		for(int d = 0; d < m; d++)
		{
			for(int i = 0; i < numberOfClients; i++)
			{
				readyTimeCon.add(startingTime[i][d] + bigM * (1 - isVisited[i][d]) >= nodes[i].getReadyTime());
				dueDateCon.add(startingTime[i][d] - bigM * (1 - isVisited[i][d]) <= nodes[i].getDueDate());
				noVisitStartingTimeCon.add(startingTime[i][d] - bigM * isVisited[i][d] <= 0);
			}
		}

		model.add(readyTimeCon);
		model.add(dueDateCon);
		model.add(noVisitStartingTimeCon);



		IloCplex cplex(model);

		clock_t sT, eT;
		sT = clock();

		//with this method, the mip is started
		cplex.solve();

		eT = clock();

		solutionTime = (double)((double)(eT-sT) / CLOCKS_PER_SEC);

		//the objective value is retrieved
		value = cplex.getObjValue();

		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;

		IloNumArray vals(env);

		vector<vector<vector<bool> > > usedArcs;

		for(int i = 0; i < numberOfClients; i++)
		{
			vector<vector<bool> > tmpVecVec;
			for(int j = 0; j < numberOfClients; j++)
			{
				vector<bool> tmpVec;
				cplex.getValues(vals, arcUsed[i][j]);
				for(int d = 0; d < m; d++)
				{
					tmpVec.push_back(vals[d]);
				}
				tmpVecVec.push_back(tmpVec);
			}
			usedArcs.push_back(tmpVecVec);
		}

		/*for(int i = 0; i < numberOfClients; i++)
		{
			for(int j = 0; j < numberOfClients; j++)
			{
				cplex.getValues(vals, arcUsed[i][j]);
				env.out() << "From Node " << i << " to Node " << j << ": " << vals << endl;
			}
		}*/
		//system("pause");

		//the optimal solution is printed on the terminal
		for(int d = 0; d < m; d++)
		{
			cout << "Vehicle " << d + 1 << ": " << endl;
			int j = 0;
			cout << "Depot -> ";
			int breakOutOfInfinity = 0;
			while(j < numberOfClients - 1 && breakOutOfInfinity < 500)
			{
				for(int i = 0; i < numberOfClients; i++)
				{
					if(usedArcs[j][i][d])
					{
						if(i < numberOfClients - 1) cout << i << " -> ";
						j = i;
					}
				}
				breakOutOfInfinity++;
			}
			cout << "Depot" << endl;

			for(int i = 0; i < numberOfClients; i++)
			{
				cplex.getValues(vals, startingTime[i]);
				if(i == 0 || i == numberOfClients-1) cout << "Depot --> \t";
				else cout << "Vertic " << i << " --> \t";
				env.out() << "ST Values        = " << vals[d] << "\t";
				cplex.getValues(vals, durationOfStay[i]);
				env.out() << "DS Values        = " << vals[d] << "\t";
				cplex.getValues(vals, isVisited[i]);
				env.out() << "IV Values        = " << vals[d] << endl;
			}
			cout << "------------------" << endl;
		}
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();

	return value;
}