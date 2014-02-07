#include "flowstarDll.h"
using namespace std;

DLLEXPORT ContinuousSystem* CDECL CreateContinuousSystem (int numVars)
{
	ContinuousSystem* system = new ContinuousSystem();
	for (int i = 0; i < numVars; ++i)
	{
		system->uncertainties.push_back(Interval(0));
		system->uncertainty_centers.push_back(Interval(0));
	}
	return system;
}
DLLEXPORT void CDECL DumpContinuousSystem(ContinuousSystem* system, int numVars, const char** varNames)
{
	vector<string> odeVarNames;
	vector<string> tmVarNames;
	vector<string> tmVarNamesWithout0;
	vector<string> localVarNames;
	tmVarNames.push_back("local_t");
	tmVarNamesWithout0.push_back("local_t");
	for (int i = 0; i < numVars; ++i)
	{
		localVarNames.push_back(varNames[i]);
		tmVarNamesWithout0.push_back(varNames[i]);		
		tmVarNames.push_back((string)(varNames[i]) + "0");
		odeVarNames.push_back((string)(varNames[i]) + "'");
	}
	cout << "Initial set:" << endl;
	system->initialSet.dump(stdout, localVarNames, tmVarNames);
	cout << endl << "ODE:" << endl;
	system->tmvOde.dump_constant(stdout, odeVarNames, tmVarNamesWithout0);
	cout << endl;
}
DLLEXPORT void CDECL SetODEContinuousSystem(ContinuousSystem* system, TaylorModelVec* tmvOde)
{
	system->tmvOde = TaylorModelVec(*tmvOde);
	
	int rangeDim = tmvOde->tms.size();
	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm hf;
		tmvOde->tms[i].expansion.toHornerForm(hf);
		system->hfOde.push_back(hf);
	}
}
DLLEXPORT void CDECL SetStrODEContinuousSystem(ContinuousSystem* system, int numOde, const char** strOde)
{
	system->strOde.clear();
	for (int i = 0; i < numOde; ++i)
		system->strOde.push_back(strOde[i]);
}
DLLEXPORT void CDECL SetInitialSetContinuousSystem(ContinuousSystem* system, Flowpipe* initialSet)
{
	system->initialSet = *initialSet;
}
DLLEXPORT Flowpipe* CDECL ReachLowDegreeContinuousSystem(ContinuousSystem* system, double step, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** varNames)
{
	vector<Interval> est;	
	vector<string> localVarNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		localVarNames.push_back(varNames[i]);
	}
	list<Flowpipe> flowpipes;
	bool nextEst = true;	
	if (!system->reach_low_degree(flowpipes, step, time, order, precondition, est, bPrint, localVarNames))
		return (Flowpipe*)0;
	/*while (!system->reach_low_degree(flowpipes, step, time, order, precondition, est, bPrint, localVarNames))
	{
		if (nextEst)
		{
			for (int i = 0; i < numVars; ++i)
				est[i] *= 2;
			if (est[0].sup() > 100000)
				return (Flowpipe*)0;
			nextEst = false;
			string str;
			est[0].toString(str);		
			cout << "Increasing error estimate to " << str << endl;
		}
		else
		{
			step /= 2;
			nextEst = true;
			cout << "Decreasing step to " << step << endl;
		}
	}*/
	
	// keep the last estimations
	for (int i = 0; i < numVars; ++i)
		estimation[i] = est[i].sup();
	
	return new Flowpipe(flowpipes.back());
}
DLLEXPORT Flowpipe* CDECL ReachLowDegreeAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames)
{
	vector<Interval> est;
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> flowpipes;
	if (!system->reach_low_degree(flowpipes, step, *miniStep, time, order, precondition, est, bPrint, varNames))
		return (Flowpipe*)0;
	
	for (int i = 0; i < numVars; ++i)
		estimation[i] = est[i].sup();
	
	return new Flowpipe(flowpipes.back());
}
DLLEXPORT Flowpipe* CDECL ReachLowDegreeAdaptiveOrderContinuousSystem(ContinuousSystem* system, double step, double time, int order, int maxOrder, int numVars, double* estimation, const char** stateVarNames)
{
	vector<Interval> est;	
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> results;	
	//system->reach_low_degree(results, step, time, order, maxOrder, ID_PRE, est, false, varNames);
	return new Flowpipe(results.back());
}
DLLEXPORT Flowpipe* CDECL ReachHighDegreeContinuousSystem(ContinuousSystem* system, double step, double time, int order, int numVars, double* estimation, const char** stateVarNames)
{
	vector<Interval> est;	
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> results;	
	system->reach_high_degree(results, step, time, order, ID_PRE, est, false, varNames);
	return new Flowpipe(results.back());
}
DLLEXPORT Flowpipe* CDECL ReachHighDegreeAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int numVars, double* estimation, const char** stateVarNames)
{
	vector<Interval> est;	
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> results;	
	system->reach_high_degree(results, step, *miniStep, time, order, ID_PRE, est, true, varNames);
	return new Flowpipe(results.back());
}
DLLEXPORT Flowpipe* CDECL ReachHighDegreeAdaptiveOrderContinuousSystem(ContinuousSystem* system, double step, double time, int order, int maxOrder, int numVars, double* estimation, const char** stateVarNames)
{
	vector<Interval> est;	
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> results;	
	system->reach_high_degree(results, step, time, order, maxOrder, ID_PRE, est, false, varNames);
	return new Flowpipe(results.back());
}
DLLEXPORT Flowpipe* CDECL ReachNonPolynomialContinuousSystem(ContinuousSystem* system, double step, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames)
{
	vector<Interval> est;
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> flowpipes;	
	bool nextEst = true;
	if (!system->reach_non_polynomial(flowpipes, step, time, order, precondition, est, bPrint, varNames))
		return (Flowpipe*)0;
	
	for (int i = 0; i < numVars; ++i)
		estimation[i] = est[i].sup();
	
	return new Flowpipe(flowpipes.back());	
}
DLLEXPORT Flowpipe* CDECL ReachNonPolynomialAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames)
{
	vector<Interval> est;
	vector<string> varNames;
	for (int i = 0; i < numVars; ++i)
	{
		est.push_back(Interval(-estimation[i],estimation[i]));
		varNames.push_back(stateVarNames[i]);
	}
	list<Flowpipe> flowpipes;	
	bool nextEst = true;
	if (!system->reach_non_polynomial(flowpipes, step, *miniStep, time, order, precondition, est, bPrint, varNames))
		return (Flowpipe*)0;
	
	for (int i = 0; i < numVars; ++i)
		estimation[i] = est[i].sup();
	
	return new Flowpipe(flowpipes.back());	
}
DLLEXPORT void CDECL DeleteContinuousSystem (ContinuousSystem* system)
{
	delete system;
}