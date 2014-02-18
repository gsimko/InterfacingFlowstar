#include "flowstarDll.h"
using namespace std;

DLLEXPORT Flowpipe* CDECL CreateFlowpipe (int numVars, Interval** domain, Interval* time)
{
	vector<Interval> box;
	for (int i = 0; i < numVars; ++i)
		box.push_back(*domain[i]);
	return new Flowpipe(box, *time);
	
	/*vector<string> statevarnames;
	statevarnames.push_back("t");
	statevarnames.push_back("x1");
	statevarnames.push_back("x2");	
	vector<string> dotvarnames;
	dotvarnames.push_back("x1'");
	dotvarnames.push_back("x2'");	
	cout << endl;
	((Flowpipe*)*result)->dump(stdout, dotvarnames, statevarnames);
	cout << endl;*/
}

DLLEXPORT void CDECL DeleteFlowpipe (Flowpipe* flowpipe)
{
	delete flowpipe;
}

DLLEXPORT void CDECL EvalFlowpipe (Flowpipe* flowpipe, Interval** result)
{
	vector<Interval> res;
	Flowpipe fp (flowpipe->tmvPre, flowpipe->tmv, flowpipe->domain);
	fp.domain[0] = Interval(flowpipe->domain[0].sup());
	fp.intEval(res);
	for (int i = 0; i < res.size(); ++i)
	{
		result[i] = new Interval(res[i]);
	}
}

DLLEXPORT bool CDECL AdvanceLowDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec *tmvOde, vector<HornerForm>* hfOde, double step, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
		
	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde->tms.size(); ++i)
	{
		polyODE.push_back(tmvOde->tms[i].expansion);
	}
	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order); 
	
	bool bvalid = flowpipe->advance_low_degree(**result, *hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, order, est, uncertainties);
	
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveStepLowDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec *tmvOde, vector<HornerForm>* hfOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<string> varNames;
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	
	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde->tms.size(); ++i)
	{
		polyODE.push_back(tmvOde->tms[i].expansion);
	}
	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order); 
	
	bool bvalid = flowpipe->advance_low_degree(**result, *hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, step, miniStep, order, est, uncertainties);
		
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveOrderLowDegreeFlowpipe(Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec *tmvOde, vector<HornerForm>* hfOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder); 
		
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	int modOrder = order;
		
	vector<Polynomial> polyODE;
	for(int i=0; i<tmvOde->tms.size(); ++i)
	{
		polyODE.push_back(tmvOde->tms[i].expansion);
	}
	vector<HornerForm> taylorExpansion;
	computeTaylorExpansion(taylorExpansion, polyODE, order); 
	
	bool bvalid = flowpipe->advance_low_degree(**result, *hfOde, taylorExpansion, precondition, step_exp_table, step_end_exp_table, modOrder, maxOrder, est, uncertainties);
		
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
		
	bool bvalid = flowpipe->advance_high_degree(**result, *hfOde, precondition, step_exp_table, step_end_exp_table, order, est, uncertainties);
	
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveStepHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	bool bvalid = flowpipe->advance_high_degree(**result, *hfOde, precondition, step_exp_table, step_end_exp_table, step, miniStep, order, est, uncertainties);	
	
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveOrderHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder); 
		
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	int modorder = order;
	bool bvalid = flowpipe->advance_high_degree(**result, *hfOde, precondition, step_exp_table, step_end_exp_table, modorder, maxOrder, est, uncertainties);
		
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<string> ode;
	for (int i = 0; i < numOde; ++i)
		ode.push_back(strOde[i]);
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
		
	bool bvalid = flowpipe->advance_non_polynomial_taylor(**result, ode, precondition, step_exp_table, step_end_exp_table, order, est, uncertainties, uncertainties);
	
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveStepNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*order); 
		
	vector<string> ode;
	for (int i = 0; i < numOde; ++i)
		ode.push_back(strOde[i]);
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	bool bvalid = flowpipe->advance_non_polynomial_taylor(**result, ode, precondition, step_exp_table, step_end_exp_table, step, miniStep, order, est, uncertainties, uncertainties);	
	
	return bvalid;
}

DLLEXPORT bool CDECL AdvanceAdaptiveOrderNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation)
{
	*result = new Flowpipe();
	vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, step, 2*maxOrder); 
		
	vector<string> ode;
	for (int i = 0; i < numOde; ++i)
		ode.push_back(strOde[i]);
	vector<Interval> est;
	for (int i = 0; i < numVars; ++i)
		est.push_back(Interval(-estimation[i],estimation[i]));
	vector<Interval> uncertainties;
	for (int i = 0; i < numVars; ++i)
		uncertainties.push_back(Interval(0));
	int modorder = order;
	bool bvalid = flowpipe->advance_non_polynomial_taylor(**result, ode, precondition, step_exp_table, step_end_exp_table, modorder, maxOrder, est, uncertainties, uncertainties);
		
	return bvalid;
}

DLLEXPORT void CDECL DumpFlowpipe (Flowpipe* A, int numStateVars, const char** stateVarNames, int numTmVars, const char** tmVarNames)
{
	vector<string> stateVars;
	for (int i = 0; i < numStateVars; ++i)
		stateVars.push_back(stateVarNames[i]);
	vector<string> tmVars;
	for (int i = 0; i < numTmVars; ++i)
		tmVars.push_back(tmVarNames[i]);
	A->dump(stdout, stateVars, tmVars);
}
