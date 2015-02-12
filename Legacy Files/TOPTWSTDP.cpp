#include "data.h"
#include "dynpro.h"
#include "hssa.h"
//#include "mip.h"

#include <boost/lambda/lambda.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <iterator>	
#include <algorithm>
#include <string>
#include <list>

using namespace std;
using boost::property_tree::ptree;

int main(int argc, char** argv)
{
	//to start the programm, you have to use the following command (example data): TOPTWSTDP "Path to the instance folder" 10 1 1 10 100 "r101.txt" "r102.txt" ... (other instance names)

	srand(time(NULL));

	

	string filePathConfig = argv[1];
	filePathConfig.append("\/");
	filePathConfig.append("Config.xml");

	ptree pt;
	read_xml(filePathConfig, pt);

	int instanceSizeInt = pt.get<int>("config.instanceParameters.instanceSize");
	int numberOfVehiclesInt = pt.get<int>("config.instanceParameters.vehicles");
	int programNoInt = pt.get<int>("config.appParameters.program");
	int selIterInt = pt.get<int>("config.saParameters.selIter");
	int iIterBInt = pt.get<int>("config.saParameters.iIterB");
	string problemType = pt.get<string>("config.instanceParameters.problemType");
	int decimalPlaces = pt.get<int>("config.instanceParameters.decimalPlaces");
	int minServiceTime = pt.get<int>("config.instanceParameters.minServiceTime");
	int maxServiceTime = pt.get<int>("config.instanceParameters.maxServiceTime");

	list<string> instancePaths;
	string pathToInstances = argv[2];
	pathToInstances.append("\/");

	BOOST_FOREACH(ptree::value_type &v, pt.get_child("config.instances")) instancePaths.push_back(pathToInstances + v.second.data());
	copy(instancePaths.begin(), instancePaths.end(),
		ostream_iterator<string>(cout," "));
	cout << endl << "Size:" << instancePaths.size() << endl;
	

	for (list<string>::iterator instancePath = instancePaths.begin(); instancePath != instancePaths.end(); instancePath++)
	{
		Instance inst(*instancePath, problemType, decimalPlaces, minServiceTime, maxServiceTime);
		cout << "------------------------------------------------------------------" << endl;
		cout << endl << "Instance: " << inst.getName() << endl << endl;
		cout << "Node ->  ";
		for (int i = 0; i < inst.getInstanceSize(); i++)
		{
			cout << i << " ";
			if (i < 10) cout << " ";
		}

		cout << endl;
		cout << endl;
		for (int i = 0; i < inst.getInstanceSize(); i++)
		{
			cout << "Node " << i << ": ";
			if (i < 10) cout << " ";
			for (int j = 0; j < inst.getInstanceSize(); j++)
			{
				cout << inst.getNode(i).getDistanceToNode(j);
				cout << " ";
				if (inst.getNode(i).getDistanceToNode(j) < 10) cout << " ";
			}
			cout << endl;
		}

		cout << endl;
		cout << endl;

		SimA fhssa(inst);
		vector<Sequence> sequences;

		//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
		int runTimeFHSSA = numberOfVehiclesInt * 100;

		//the FHSSA heuristic is started, the best found solution of the FHSSA is stored in tmpFHSSAValue, afterwards, the results are printed on the command prompt
		double FHSSAValue = fhssa.startSAOld(numberOfVehiclesInt, 1, runTimeFHSSA, 1, 0.1, 0.999, iIterBInt, selIterInt, 100000);
		sequences = fhssa.getSequences();
		double durationOfFHSSA = fhssa.getDurationOfSimAn();

		cout << "+++++++++++++++++++++++" << endl << endl;
		cout << "Fast Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
		cout << "Duration of SA-Heuristic: " << durationOfFHSSA << endl;
		cout << "Sequences:" << endl << endl;
		for (int i = 0; i < sequences.size(); i++)
		{
			cout << "Sequence " << i + 1 << ":" << endl;
			cout << "Number\tStarting Time\tDuration" << endl;
			for (int j = 0; j < sequences[i].getNodes().size(); j++)
			{
				cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
			}
			cout << endl;
		}

		cout << "Best Value: " << FHSSAValue << endl;

	}
		
	//system("pause");

	//-----------------------------------------------------------------------------------------------------------------------------
	//if(programNoInt == 4)
	//{
	//	for(int instNumber = 7; instNumber < argc; instNumber++)
	//	{
	//		string filePath = argv[1];
	//		filePath.append("\/");
	//		filePath.append(argv[instNumber]);
	//		Instance inst(filePath, 1);
	//		cout << "------------------------------------------------------------------" << endl;
	//			cout << endl << "Instance: " << inst.getName() << endl << endl;
	//			cout << "Node ->  ";
	//			for(int i = 0; i < inst.getInstanceSize(); i++)
	//			{
	//				cout << i << " ";
	//				if(i < 10) cout << " ";
	//			}

	//			cout << endl;
	//			cout << endl;
	//			for(int i = 0; i < inst.getInstanceSize(); i++)
	//			{
	//				cout << "Node " << i << ": ";
	//				if(i < 10) cout << " ";
	//				for(int j = 0; j < inst.getInstanceSize(); j++)
	//				{
	//					cout << inst.getNode(i).getDistanceToNode(j);
	//					cout << " ";
	//					if(inst.getNode(i).getDistanceToNode(j) < 10) cout << " ";
	//				}
	//				cout << endl;
	//			}

	//			cout << endl;
	//			cout << endl;

	//			SimA fhssa(inst);
	//			vector<Sequence> sequences;

	//			//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
	//			int runTimeFHSSA = numberOfVehiclesInt * 100;

	//			//the FHSSA heuristic is started, the best found solution of the FHSSA is stored in tmpFHSSAValue, afterwards, the results are printed on the command prompt
	//			double FHSSAValue = fhssa.startSA(numberOfVehiclesInt,1,runTimeFHSSA,1,0.1,0.999,iIterBInt,selIterInt,100000);
	//			sequences = fhssa.getSequences();
	//			double durationOfFHSSA = fhssa.getDurationOfSimAn();

	//			cout << "+++++++++++++++++++++++" << endl << endl;
	//			cout << "Fast Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
	//			cout << "Duration of SA-Heuristic: " << durationOfFHSSA << endl;
	//			cout << "Sequences:" << endl << endl;
	//			for(int i = 0; i < sequences.size(); i++)
	//			{
	//				cout << "Sequence " << i + 1 << ":" << endl;
	//				cout << "Number\tStarting Time\tDuration" << endl;
	//				for(int j = 0; j < sequences[i].getNodes().size(); j++)
	//				{
	//					cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
	//				}
	//				cout << endl;
	//			}

	//			cout << "Best Value: " << FHSSAValue << endl;
	//			
	//	}
	//}
	//else if(programNoInt == 7)
	//{
	//	vector<string> instNames;
	//	vector<double> mipValues;
	//	vector<double> mipDurations;
	//	vector<double> fhssaValues;
	//	vector<double> fhssaDurations;
	//	vector<double> fhssaGaps;
	//	vector<double> shssaValues;
	//	vector<double> shssaDurations;
	//	vector<double> shssaGaps;
	//	vector<double> fhssaBestValue;
	//	vector<double> shssaBestValue;

	//	for(int instNumber = 7; instNumber < argc; instNumber++)
	//	{
	//		string filePath = argv[1];
	//		filePath.append("\/");
	//		filePath.append(argv[instNumber]);

	//		//partitioning of instance
	//		int amountOfInstances = (int) floor((double) (100 / instanceSizeInt));

	//		//iteration through all instance parts of size instanceSizeInt
	//		for(int instMultiFirstNode = 0; instMultiFirstNode < amountOfInstances; instMultiFirstNode++)
	//		{
	//			Instance inst(filePath,1 + (instMultiFirstNode*instanceSizeInt),instanceSizeInt,5,15);

	//			cout << endl << endl;
	//			cout << "------------------------------------------------------------------" << endl;
	//			cout << endl << "Instance: " << inst.getName() << endl << endl;
	//			cout << "Node ->  ";
	//			for(int i = 0; i < inst.getInstanceSize(); i++)
	//			{
	//				cout << i << " ";
	//				if(i < 10) cout << " ";
	//			}

	//			cout << endl;
	//			cout << endl;
	//			for(int i = 0; i < inst.getInstanceSize(); i++)
	//			{
	//				cout << "Node " << i << ": ";
	//				if(i < 10) cout << " ";
	//				for(int j = 0; j < inst.getInstanceSize(); j++)
	//				{
	//					cout << inst.getNode(i).getDistanceToNode(j);
	//					cout << " ";
	//					if(inst.getNode(i).getDistanceToNode(j) < 10) cout << " ";
	//				}
	//				cout << endl;
	//			}

	//			cout << endl;
	//			cout << endl;


	//			//here, the MIP is initialized and executed, the results are later on printed on the command prompt
	//			cout << "MIP: " << endl;

	//			MIP mip(inst);
	//			double mipValue = 100000;
	//			double mipDuration = 30;

	//			mipValue = mip.solve(numberOfVehiclesInt);
	//			mipDuration = mip.getSolutionTime();


	//			cout << "MIP Objective Value: " << mipValue << endl;


	//			//here, 5 runs for each instance are performed for the SHSSA and the FHSSA
	//			int runs = 5;
	//			double fhssaValue = 0;
	//			double durationOfFHSSA = 0;
	//			double shssaValue = 0;
	//		    double durationOfSHSSA = 0;
	//			double bestFHSSAValue = 0;
	//			double bestSHSSAValue = 0;

	//			//iteration through runs
	//			for(int numberOfRun = 0; numberOfRun < runs; numberOfRun++)
	//			{
	//				double tmpFHSSAValue = 0;
	//				double tmpSHSSAValue = 0;

	//				SimA fhssa(inst);
	//				vector<Sequence> sequences;

	//				//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
	//				int runTimeFHSSA = numberOfVehiclesInt * instanceSizeInt;

	//				//the FHSSA heuristic is started, the best found solution of the FHSSA is stored in tmpFHSSAValue, afterwards, the results are printed on the command prompt
	//				tmpFHSSAValue = fhssa.startSA(numberOfVehiclesInt,1,runTimeFHSSA,1,0.1,0.999,iIterBInt,selIterInt,mipValue);
	//				if(tmpFHSSAValue > bestFHSSAValue) bestFHSSAValue = tmpFHSSAValue;
	//				fhssaValue += tmpFHSSAValue;
	//				sequences = fhssa.getSequences();
	//				durationOfFHSSA += fhssa.getDurationOfSimAn();

	//				cout << "+++++++++++++++++++++++" << endl << endl;
	//				cout << "Fast Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
	//				cout << "Duration of SA-Heuristic: " << durationOfFHSSA << endl;
	//				cout << "Sequences:" << endl << endl;
	//				for(int i = 0; i < sequences.size(); i++)
	//				{
	//					cout << "Sequence " << i + 1 << ":" << endl;
	//					cout << "Number\tStarting Time\tDuration" << endl;
	//					for(int j = 0; j < sequences[i].getNodes().size(); j++)
	//					{
	//						cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
	//					}
	//					cout << endl;
	//				}

	//				cout << "Best Value: " << fhssaValue << endl;


	//				//here, the same happens again for the SHSSA version of the heuristic
	//				SimA shssa(inst);

	//				//(int M, int SAType, int Parameter, int Variant, double T, double Alpha, int IIterB, int ISel, double OptimalValue)
	//				tmpSHSSAValue = shssa.startSA(numberOfVehiclesInt,2,30,1,0.3,0.99,iIterBInt,selIterInt,mipValue);
	//				if(tmpSHSSAValue > bestSHSSAValue) bestSHSSAValue = tmpSHSSAValue;
	//				shssaValue += tmpSHSSAValue;
	//				sequences = shssa.getSequences();
	//				durationOfSHSSA  += shssa.getDurationOfSimAn();

	//				cout << "+++++++++++++++++++++++" << endl << endl;
	//				cout << "Slow Hybrid Selective Simulated Annealing Heuristic: " << endl << endl;
	//				cout << "Duration of SA-Heuristic: " << durationOfSHSSA << endl;
	//				cout << "Sequences:" << endl << endl;
	//				for(int i = 0; i < sequences.size(); i++)
	//				{
	//					cout << "Sequence " << i + 1 << ":" << endl;
	//					cout << "Number\tStarting Time\tDuration" << endl;
	//					for(int j = 0; j < sequences[i].getNodes().size(); j++)
	//					{
	//						cout << sequences[i].getNodes()[j].getNumberOfNode() << "\t" << sequences[i].getNodes()[j].getStartingTime() << "\t\t" << sequences[i].getNodes()[j].getDuration() << endl;
	//					}
	//					cout << endl;
	//				}

	//				cout << "Best Value: " << shssaValue << endl;
	//			}

	//			shssaValue = shssaValue / runs;
	//			durationOfSHSSA = durationOfSHSSA / runs;
	//			fhssaValue = fhssaValue / runs;
	//			durationOfFHSSA = durationOfFHSSA / runs;

	//			//all values are stored in the respective vectors, which later on are printed into a file
	//			instNames.push_back(inst.getName());
	//			mipValues.push_back(mipValue);
	//			mipDurations.push_back(mipDuration);
	//			fhssaValues.push_back(fhssaValue);
	//			fhssaDurations.push_back(durationOfFHSSA);
	//			fhssaGaps.push_back(1 - (fhssaValue/mipValue));
	//			shssaGaps.push_back(1 - (shssaValue/mipValue));
	//			shssaValues.push_back(shssaValue);
	//			shssaDurations.push_back(durationOfSHSSA);
	//			fhssaBestValue.push_back(bestFHSSAValue);
	//			shssaBestValue.push_back(bestSHSSAValue);

	//			//periodically after each instance, an update file is written, since the cluster takes a while to solve all instances
	//			ofstream latexFile;
	//			stringstream instanceSizeIntS;
	//			instanceSizeIntS << instanceSizeInt;
	//			stringstream numberOfVehiclesIntS;
	//			numberOfVehiclesIntS << numberOfVehiclesInt;
	//			string instanceName = inst.getName();
	//			stringstream instanceNameS;
	//			instanceNameS << instanceName;

	//			string latexFileName = "Update_Benchmark_" + instanceSizeIntS.str() + "_" + numberOfVehiclesIntS.str() + "_" + ".txt";
	//			latexFile.open(latexFileName);
	//			for(int i = 0; i < instNames.size(); i++)
	//			{
	//				latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << mipValues[i] <<  "&" << mipDurations[i] << "&&" << fhssaBestValue[i] << "&" << fhssaValues[i] << "&" << fhssaGaps[i] << "&" << fhssaDurations[i] << "&&" << shssaBestValue[i] << "&" <<shssaValues[i] << "&" << shssaGaps[i]  << "&" << shssaDurations[i] << "\\\\" << endl;
	//			}

	//			double averageFHSSAGap = 0;
	//			double averageSHSSAGap = 0;
	//			double maxFHSSAGap = 0;
	//			double maxSHSSAGap = 0;

	//			for(int i = 0; i < fhssaGaps.size(); i++)
	//			{
	//				if(fhssaGaps[i] > maxFHSSAGap) maxFHSSAGap = fhssaGaps[i];
	//				if(shssaGaps[i] > maxSHSSAGap) maxSHSSAGap = shssaGaps[i];
	//				averageFHSSAGap += fhssaGaps[i];
	//				averageSHSSAGap += shssaGaps[i];
	//			}
	//			averageFHSSAGap = averageFHSSAGap / fhssaGaps.size();
	//			averageSHSSAGap = averageSHSSAGap / shssaGaps.size();

	//			latexFile << "\\hline" << endl;

	//			latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average &&&&&&" << averageFHSSAGap << "&&&&&" <<  averageSHSSAGap << "\\\\" << endl;

	//			latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Max &&&&&&" << maxFHSSAGap << "&&&&&" << maxSHSSAGap << "\\\\" << endl;

	//			latexFile.close();



	//		}
	//	}

	//	cout << endl;
	//	//Table for LaTex

	//	//this file writes all results into a txt file, which can later on be directly copied into a latex table
	//	ofstream latexFile;
	//	stringstream instanceSizeIntS;
	//	instanceSizeIntS << instanceSizeInt;
	//	stringstream numberOfVehiclesIntS;
	//	numberOfVehiclesIntS << numberOfVehiclesInt;
	//	string latexFileName = "LaTexTable_Benchmark_" + instanceSizeIntS.str() + "_" + numberOfVehiclesIntS.str() + ".txt";
	//	latexFile.open(latexFileName);
	//	for(int i = 0; i < instNames.size(); i++)
	//	{
	//		latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << instNames[i] << "&" << mipValues[i] <<  "&" << mipDurations[i] << "&&" << fhssaBestValue[i] << "&" << fhssaValues[i] << "&" << fhssaGaps[i] << "&" << fhssaDurations[i] << "&&" << shssaBestValue[i] << "&" <<shssaValues[i] << "&" << shssaGaps[i]  << "&" << shssaDurations[i] << "\\\\" << endl;
	//	}

	//	double averageFHSSAGap = 0;
	//	double averageSHSSAGap = 0;
	//	double maxFHSSAGap = 0;
	//	double maxSHSSAGap = 0;

	//	for(int i = 0; i < fhssaGaps.size(); i++)
	//	{
	//		if(fhssaGaps[i] > maxFHSSAGap) maxFHSSAGap = fhssaGaps[i];
	//		if(shssaGaps[i] > maxSHSSAGap) maxSHSSAGap = shssaGaps[i];
	//		averageFHSSAGap += fhssaGaps[i];
	//		averageSHSSAGap += shssaGaps[i];
	//	}
	//	averageFHSSAGap = averageFHSSAGap / fhssaGaps.size();
	//	averageSHSSAGap = averageSHSSAGap / shssaGaps.size();

	//	latexFile << "\\hline" << endl;

	//	latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Average &&&&&&" << averageFHSSAGap << "&&&&&" <<  averageSHSSAGap << "\\\\" << endl;

	//	latexFile << setw(2) << setfill('0') << fixed << setprecision(2) << "Max &&&&&&" << maxFHSSAGap << "&&&&&" << maxSHSSAGap << "\\\\" << endl;

	//	latexFile.close();
	//}

	//
	//system("pause");
	return 0;
}