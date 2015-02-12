#ifndef BENCHRUN_H
#define BENCHRUN_H

//if the variant entered in the command prompt equals 1, benchmark runs are performed, if the variant is 2, preliminary checks are performed (see paper for the meaning of them)
/*	if(programNoInt == 1)
	{
		vector<string> instNames;
		vector<double> mipValues;
		vector<double> mipDurations;
		vector<double> fhssaValues;
		vector<double> fhssaDurations;
		vector<double> fhssaGaps;
		vector<double> shssaValues;
		vector<double> shssaDurations;
		vector<double> shssaGaps;
		vector<double> fhssaBestValue;
		vector<double> shssaBestValue;

		for(int instNumber = 7; instNumber < argc; instNumber++)
		{
			string filePath = argv[1];
			filePath.append("\/");
			filePath.append(argv[instNumber]);

			//partitioning of instance
			int amountOfInstances = (int) floor((double) (100 / instanceSizeInt));

			//iteration through all instance parts of size instanceSizeInt
			for(int instMultiFirstNode = 0; instMultiFirstNode < amountOfInstances; instMultiFirstNode++)
			{
				Instance inst(filePath,1 + (instMultiFirstNode*instanceSizeInt),instanceSizeInt,5,15);

				cout << endl << endl;
				cout << "------------------------------------------------------------------" << endl;
				cout << endl << "Instance: " << inst.getName() << endl << endl;
				cout << "Node ->  ";
				for(int i = 0; i < inst.getInstanceSize(); i++)
				{
					cout << i << " ";
					if(i < 10) cout << " ";
				}

				cout << endl;
				cout << endl;
				for(int i = 0; i < inst.getInstanceSize(); i++)
				{
					cout << "Node " << i << ": ";
					if(i < 10) cout << " ";
					for(int j = 0; j < inst.getInstanceSize(); j++)
					{
						cout << inst.getNode(i).getDistanceToNode(j);
						cout << " ";
						if(inst.getNode(i).getDistanceToNode(j) < 10) cout << " ";
					}
					cout << endl;
				}

				cout << endl;
				cout << endl;


				//here, the MIP is initialized and executed, the results are later on printed on the command prompt
				cout << "MIP: " << endl;

				MIP mip(inst);
				double mipValue = 100000;
				double mipDuration = 30;

				mipValue = mip.solve(numberOfVehiclesInt);
				mipDuration = mip.getSolutionTime();


				cout << "MIP Objective Value: " << mipValue << endl;


				//here, 5 runs for each instance are performed for the SHSSA and the FHSSA
				int runs = 5;
				double fhssaValue = 0;
				double durationOfFHSSA = 0;
				double shssaValue = 0;
			    double durationOfSHSSA = 0;
				double bestFHSSAValue = 0;
				double bestSHSSAValue = 0;

				//iteration through runs
				for(int numberOfRun = 0; numberOfRun < runs; numberOfRun++)
				{
					double tmpFHSSAValue = 0;
					double tmpSHSSAValue = 0;

					SimA fhssa(inst);
					vector<Sequence> sequences;

					//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
					int runTimeFHSSA = numberOfVehiclesInt * instanceSizeInt;

					//the FHSSA heuristic is started, the best found solution of the FHSSA is stored in tmpFHSSAValue, afterwards, the results are printed on the command prompt
					tmpFHSSAValue = fhssa.startSA(numberOfVehiclesInt,1,runTimeFHSSA,1,0.1,0.999,iIterBInt,selIterInt,mipValue);
					if(tmpFHSSAValue > bestFHSSAValue) bestFHSSAValue = tmpFHSSAValue;
					fhssaValue += tmpFHSSAValue;
					sequences = fhssa.getSequences();
					durationOfFHSSA += fhssa.getDurationOfSimAn();

					cout << "+++++++++++++++++++++++" << endl << endl;
					cout << "Fast Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
					cout << "Duration of SA-Heuristic: " << durationOfFHSSA << endl;
					cout << "Sequences:" << endl << endl;
					for(int i = 0; i < sequences.size(); i++)
					{
						cout << "Sequence " << i + 1 << ":" << endl;
						cout << "Number\tStarting Time\tDuration" << endl;
						for(int j = 0; j < sequences[i].getNodes().size(); j++)
						{
							cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
						}
						cout << endl;
					}

					cout << "Best Value: " << fhssaValue << endl;


					//here, the same happens again for the SHSSA version of the heuristic
					SimA shssa(inst);

					//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
					tmpSHSSAValue = shssa.startSA(numberOfVehiclesInt,2,30,1,0.3,0.99,iIterBInt,selIterInt,mipValue);
					if(tmpSHSSAValue > bestSHSSAValue) bestSHSSAValue = tmpSHSSAValue;
					shssaValue += tmpSHSSAValue;
					sequences = shssa.getSequences();
					durationOfSHSSA  += shssa.getDurationOfSimAn();

					cout << "+++++++++++++++++++++++" << endl << endl;
					cout << "Slow Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
					cout << "Duration of SA-Heuristic: " << durationOfSHSSA << endl;
					cout << "Sequences:" << endl << endl;
					for(int i = 0; i < sequences.size(); i++)
					{
						cout << "Sequence " << i + 1 << ":" << endl;
						cout << "Number\tStarting Time\tDuration" << endl;
						for(int j = 0; j < sequences[i].getNodes().size(); j++)
						{
							cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
						}
						cout << endl;
					}

					cout << "Best Value: " << shssaValue << endl;
				}

				shssaValue = shssaValue / runs;
				durationOfSHSSA = durationOfSHSSA / runs;
				fhssaValue = fhssaValue / runs;
				durationOfFHSSA = durationOfFHSSA / runs;

				//all values are stored in the respective vectors, which later on are printed into a file
				instNames.push_back(inst.getName());
				mipValues.push_back(mipValue);
				mipDurations.push_back(mipDuration);
				fhssaValues.push_back(fhssaValue);
				fhssaDurations.push_back(durationOfFHSSA);
				fhssaGaps.push_back(1 - (fhssaValue/mipValue));
				shssaGaps.push_back(1 - (shssaValue/mipValue));
				shssaValues.push_back(shssaValue);
				shssaDurations.push_back(durationOfSHSSA);
				fhssaBestValue.push_back(bestFHSSAValue);
				shssaBestValue.push_back(bestSHSSAValue);

				//periodically after each instance, an update file is written, since the cluster takes a while to solve all instances
				ofstream latexFile;
				stringstream instanceSizeIntS;
				instanceSizeIntS << instanceSizeInt;
				stringstream numberOfVehiclesIntS;
				numberOfVehiclesIntS << numberOfVehiclesInt;
				string instanceName = inst.getName();
				stringstream instanceNameS;
				instanceNameS << instanceName;

				string latexFileName = "Update_Benchmark_" + instanceSizeIntS.str() + "_" + numberOfVehiclesIntS.str() + "_" + ".txt";
				latexFile.open(latexFileName);
				for(int i = 0; i < instNames.size(); i++)
				{
					latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << mipValues[i] <<  "&" << mipDurations[i] << "&&" << fhssaBestValue[i] << "&" << fhssaValues[i] << "&" << fhssaGaps[i] << "&" << fhssaDurations[i] << "&&" << shssaBestValue[i] << "&" <<shssaValues[i] << "&" << shssaGaps[i]  << "&" << shssaDurations[i] << "\\\\" << endl;
				}

				double averageFHSSAGap = 0;
				double averageSHSSAGap = 0;
				double maxFHSSAGap = 0;
				double maxSHSSAGap = 0;

				for(int i = 0; i < fhssaGaps.size(); i++)
				{
					if(fhssaGaps[i] > maxFHSSAGap) maxFHSSAGap = fhssaGaps[i];
					if(shssaGaps[i] > maxSHSSAGap) maxSHSSAGap = shssaGaps[i];
					averageFHSSAGap += fhssaGaps[i];
					averageSHSSAGap += shssaGaps[i];
				}
				averageFHSSAGap = averageFHSSAGap / fhssaGaps.size();
				averageSHSSAGap = averageSHSSAGap / shssaGaps.size();

				latexFile << "\\hline" << endl;

				latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average &&&&&&" << averageFHSSAGap << "&&&&&" <<  averageSHSSAGap << "\\\\" << endl;

				latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Max &&&&&&" << maxFHSSAGap << "&&&&&" << maxSHSSAGap << "\\\\" << endl;

				latexFile.close();



			}
		}

		cout << endl;
		//Table for LaTex

		//this file writes all results into a txt file, which can later on be directly copied into a latex table
		ofstream latexFile;
		stringstream instanceSizeIntS;
		instanceSizeIntS << instanceSizeInt;
		stringstream numberOfVehiclesIntS;
		numberOfVehiclesIntS << numberOfVehiclesInt;
		string latexFileName = "LaTexTable_Benchmark_" + instanceSizeIntS.str() + "_" + numberOfVehiclesIntS.str() + ".txt";
		latexFile.open(latexFileName);
		for(int i = 0; i < instNames.size(); i++)
		{
			latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << mipValues[i] <<  "&" << mipDurations[i] << "&&" << fhssaBestValue[i] << "&" << fhssaValues[i] << "&" << fhssaGaps[i] << "&" << fhssaDurations[i] << "&&" << shssaBestValue[i] << "&" <<shssaValues[i] << "&" << shssaGaps[i]  << "&" << shssaDurations[i] << "\\\\" << endl;
		}

		double averageFHSSAGap = 0;
		double averageSHSSAGap = 0;
		double maxFHSSAGap = 0;
		double maxSHSSAGap = 0;

		for(int i = 0; i < fhssaGaps.size(); i++)
		{
			if(fhssaGaps[i] > maxFHSSAGap) maxFHSSAGap = fhssaGaps[i];
			if(shssaGaps[i] > maxSHSSAGap) maxSHSSAGap = shssaGaps[i];
			averageFHSSAGap += fhssaGaps[i];
			averageSHSSAGap += shssaGaps[i];
		}
		averageFHSSAGap = averageFHSSAGap / fhssaGaps.size();
		averageSHSSAGap = averageSHSSAGap / shssaGaps.size();

		latexFile << "\\hline" << endl;

		latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average &&&&&&" << averageFHSSAGap << "&&&&&" <<  averageSHSSAGap << "\\\\" << endl;

		latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Max &&&&&&" << maxFHSSAGap << "&&&&&" << maxSHSSAGap << "\\\\" << endl;

		latexFile.close();
	}
	else if(programNoInt == 2)
	{
		//this is the preliminary check run
		vector<string> instNames;
		vector<double> fhssaBestValue;
		vector<double> shssaBestValue;
		vector<double> fhssaValues;
		vector<double> shssaValues;
		vector<double> fhssaDurations;
		vector<double> shssaDurations;
		vector<int> numberOfVehicles;
		//from vehicle numbers 1 till 4, all scenarios are played through
		for(int m = 1; m < 5; m++)
		{
			for(int instNumber = 7; instNumber < argc; instNumber++)
			{
				string filePath = argv[1];
				filePath.append("\/");
				filePath.append(argv[instNumber]);

				//"C:\\Users\\PascalRauprecht\\Dropbox\\Universität\\Viadrina Master\\Management Science\\Seminar\\Programm\\Instances\\Vansteenwegen2009\\r101.txt"

				//Instance(string Path, int FirstNode, int NumberOfNodes, int MinServiceTime, int MaxServiceTime);
				vector<int> instanceSizes;
				instanceSizes.push_back(10);
				instanceSizes.push_back(20);
				instanceSizes.push_back(30);
				instanceSizes.push_back(50);
				instanceSizes.push_back(100);
				for(int i = 0; i < instanceSizes.size(); i++)
				{
					instanceSizeInt = instanceSizes[i];
					int amountOfInstances = (int) floor((double) (100 / instanceSizeInt));
					//if the partitioned amount of instances is larger than 2, it is limited to two, due to computational time reasons
					if(amountOfInstances > 2) amountOfInstances = 2;
					for(int instMultiFirstNode = 0; instMultiFirstNode < amountOfInstances; instMultiFirstNode++)
					{
						Instance inst(filePath,1 + (instMultiFirstNode*instanceSizeInt),instanceSizeInt,5,15);

						
						double bestFHSSAValue = 0;
						double bestSHSSAValue = 0;
						double fhssaValue = 0;
						double shssaValue = 0;
						double averageFHSSAValue = 0;
						double averageSHSSAValue = 0;
						double durationOfFHSSA = 0;
						double durationOfSHSSA = 0;
						double averageDurationOfFHSSA = 0;
						double averageDurationOfSHSSA = 0;

						//in the pre check, 3 runs per instance were performed
						int runs = 3;
						for(int j = 0; j < runs; j++)
						{
							SimA fhssa(inst);
							
							int runTimeFHSSA = m * instanceSizeInt;
							fhssaValue = fhssa.startSA(m,1,runTimeFHSSA,1,0.1,0.999,iIterBInt,selIterInt,100000);
							durationOfFHSSA = fhssa.getDurationOfSimAn();
							averageDurationOfFHSSA += durationOfFHSSA;
							averageFHSSAValue += fhssaValue;

							if(fhssaValue > bestFHSSAValue) bestFHSSAValue = fhssaValue;

							//startSA(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
							SimA shssa(inst);
							
							shssaValue = shssa.startSA(m,2,30,1,0.3,0.99,iIterBInt,selIterInt,100000);
							durationOfSHSSA  = shssa.getDurationOfSimAn();
							averageDurationOfSHSSA += durationOfSHSSA;
							averageSHSSAValue += shssaValue;

							if(shssaValue > bestSHSSAValue) bestSHSSAValue = shssaValue;
						}

						averageDurationOfSHSSA = averageDurationOfSHSSA / runs;
						averageSHSSAValue = averageSHSSAValue / runs;
						averageDurationOfFHSSA = averageDurationOfFHSSA / runs;
						averageFHSSAValue = averageFHSSAValue / runs;

						cout << "----------------------------" << endl;
						cout << "Instance: " << inst.getName() << endl;
						cout << "Number of vehicles: " << m << endl;
						cout << "Number of vertices: " << instanceSizeInt << endl;
						cout << "averageFHSSAValue: " << averageFHSSAValue << endl;
						cout << "averageSHSSAValue: " << averageSHSSAValue << endl;
						cout << "averageDurationOfFHSSA: " << averageDurationOfFHSSA << endl;
						cout << "averageDurationOfSHSSA: " << averageDurationOfSHSSA << endl;
						cout << "bestFHSSAValue: " << bestFHSSAValue << endl;
						cout << "bestSHSSAValue: " << bestSHSSAValue << endl << endl;

						//the results are stored in the according vectors
						instNames.push_back(inst.getName());
						fhssaValues.push_back(averageFHSSAValue);
						shssaValues.push_back(averageSHSSAValue);
						fhssaDurations.push_back(averageDurationOfFHSSA);
						shssaDurations.push_back(averageDurationOfSHSSA);
						fhssaBestValue.push_back(bestFHSSAValue);
						shssaBestValue.push_back(bestSHSSAValue);
						numberOfVehicles.push_back(m);

						ofstream latexFile;
						stringstream instanceSizeIntS;
						instanceSizeIntS << instanceSizeInt;
				
						stringstream numberOfVehiclesIntS;
						numberOfVehiclesIntS << numberOfVehiclesInt;

						stringstream selIterIntS;
						selIterIntS << selIterInt;

						stringstream iIterBIntS;
						iIterBIntS << iIterBInt;

						int sizeNames = instNames.size(); 
						stringstream sizeNamesS;
						sizeNamesS << sizeNames;

						//the Update file serves as check, how far the cluster already computed the instances
						string latexFileName = "Update_PreCheck_" + iIterBIntS.str() + "_" + selIterIntS.str() + "_" +  sizeNamesS.str() + ".txt";
						latexFile.open(latexFileName);

						double averageTimeSHSSA = 0;
						double averageTimeFHSSA = 0;
						double averageValueSHSSA = 0;
						double averageValueFHSSA = 0;
						latexFile << "iIterB = " << iIterBInt << endl << endl;
						latexFile << "sIter = " << selIterInt << endl << endl;
						for(int i = 0; i < instNames.size(); i++)
						{
							averageTimeSHSSA += shssaDurations[i];
							averageTimeFHSSA += fhssaDurations[i];
							averageValueSHSSA += shssaValues[i];
							averageValueFHSSA += fhssaValues[i];

							latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << numberOfVehicles[i] << "&" << fhssaValues[i] << "&" << fhssaDurations[i] << "&&" << shssaValues[i] << "&" << shssaDurations[i] << "\\\\" << endl;
						}

						averageTimeSHSSA = averageTimeSHSSA / instNames.size();
						averageTimeFHSSA = averageTimeFHSSA / instNames.size();
						averageValueSHSSA = averageValueSHSSA / instNames.size();
						averageValueFHSSA = averageValueFHSSA / instNames.size();

						latexFile << "\\hline" << endl;

						latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average" << "&&" << averageValueFHSSA  << "&" << averageTimeFHSSA << "&&" << averageValueSHSSA << "&" << averageTimeSHSSA << "\\\\" << endl;

						latexFile << "\\hline" << endl;
						latexFile.close();

					}
				}

				//cout << endl;
				//Table for LaTex

				
			}
		}

		ofstream latexFile;
		stringstream instanceSizeIntS;
		instanceSizeIntS << instanceSizeInt;
				
		stringstream numberOfVehiclesIntS;
		numberOfVehiclesIntS << numberOfVehiclesInt;

		stringstream selIterIntS;
		selIterIntS << selIterInt;

		stringstream iIterBIntS;
		iIterBIntS << iIterBInt;

		//in the final file, all results are written as a latex table to later on make it handy to copy paste the table into the latex paper (see paper)
		string latexFileName = "LaTexTable_PreCheck_" + iIterBIntS.str() + "_" + selIterIntS.str() + ".txt";
		latexFile.open(latexFileName);

		double averageTimeSHSSA = 0;
		double averageTimeFHSSA = 0;
		double averageValueSHSSA = 0;
		double averageValueFHSSA = 0;
		latexFile << "iIterB = " << iIterBIntS << endl << endl;
		latexFile << "sIter = " << selIterInt << endl << endl;
		for(int i = 0; i < instNames.size(); i++)
		{
			averageTimeSHSSA += shssaDurations[i];
			averageTimeFHSSA += fhssaDurations[i];
			averageValueSHSSA += shssaValues[i];
			averageValueFHSSA += fhssaValues[i];

			latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << numberOfVehicles[i] << "&" << fhssaValues[i] << "&" << fhssaDurations[i] << "&&" << shssaValues[i] << "&" << shssaDurations[i] << "\\\\" << endl;
		}

		averageTimeSHSSA = averageTimeSHSSA / instNames.size();
		averageTimeFHSSA = averageTimeFHSSA / instNames.size();
		averageValueSHSSA = averageValueSHSSA / instNames.size();
		averageValueFHSSA = averageValueFHSSA / instNames.size();

		latexFile << "\\hline" << endl;

		latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average" << "&&" << averageValueFHSSA  << "&" << averageTimeFHSSA << "&&" << averageValueSHSSA << "&" << averageTimeSHSSA << "\\\\" << endl;

		latexFile << "\\hline" << endl;

		latexFile.close();
	}
*/
#endif