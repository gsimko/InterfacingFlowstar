#include "flowstarDll.h"
using namespace std;

DLLEXPORT Monomial* CDECL CreateMonomial (Interval* I, int numVars, int* degs)
{
	vector<int> deg;
	for (int i = 0; i < numVars; ++i)
		deg.push_back(degs[i]);
	return new Monomial(*I, deg);
}
DLLEXPORT Monomial* CDECL CreateConstantMonomial (Interval* I, int numVars)
{
	return new Monomial(*I, numVars);
}
DLLEXPORT void CDECL DeleteMonomial (Monomial* monomial)
{
	delete monomial;
}
DLLEXPORT Monomial* CDECL AddMonomial (Monomial* A, Monomial* B)
{
	Monomial* result = new Monomial(*A);
	*result += *B;
	return result;
}
DLLEXPORT Monomial* CDECL MulMonomial (Monomial* A, Monomial* B)
{
	Monomial* result = new Monomial(*A);
	*result *= *B;
	return result;
}
DLLEXPORT void CDECL DumpMonomial (Monomial* A, int numVars, const char** varNames, bool dumpIntervals)
{
	vector<string> vars;
	for (int i = 0; i < numVars; ++i)
		vars.push_back(varNames[i]);
	if (dumpIntervals)
		A->dump_interval(stdout, vars);
	else
		A->dump_constant(stdout, vars);
}

DLLEXPORT Polynomial* CDECL CreateEmptyPolynomial ()
{
	return new Polynomial();
}
DLLEXPORT Polynomial* CDECL CreatePolynomial (Monomial *monomial)
{
	return new Polynomial(*monomial);
}
DLLEXPORT Polynomial* CDECL CreateConstantPolynomial (Interval *constant, int numVars)
{
	return new Polynomial(*constant, numVars);
}
DLLEXPORT void CDECL DeletePolynomial (Polynomial* polynomial)
{
	delete polynomial;
}
DLLEXPORT void CDECL AddAssignPolynomial (Polynomial* polynomial, Monomial* monomial)
{	
	polynomial->add_assign(*monomial);
}

DLLEXPORT Polynomial* CDECL AddPolynomial(Polynomial* A, Polynomial* B)
{
	Polynomial* result = new Polynomial(*A);
	*result += *B;
	return result;
}
DLLEXPORT Polynomial* CDECL SubPolynomial(Polynomial* A, Polynomial* B)
{
	Polynomial* result = new Polynomial(*A);
	*result -= *B;
	return result;
}
DLLEXPORT Polynomial* CDECL MulPolynomial(Polynomial* A, Polynomial* B)
{
	Polynomial* result = new Polynomial(*A);
	*result *= *B;
	return result;
}
DLLEXPORT Polynomial* CDECL NegPolynomial(Polynomial* A)
{
	Polynomial* result = new Polynomial();
	A->inv(*result);
	return result;
}
DLLEXPORT Polynomial* CDECL SinPolynomial(Polynomial* A, int numVars, int order)
{
	Polynomial* result = new Polynomial();
	A->sin_taylor(*result, numVars, order);
	return result;
}
DLLEXPORT Polynomial* CDECL CosPolynomial(Polynomial* A, int numVars, int order)
{
	Polynomial* result = new Polynomial();
	A->cos_taylor(*result, numVars, order);
	return result;
}
DLLEXPORT Polynomial* CDECL ExpPolynomial(Polynomial* A, int numVars, int order)
{
	Polynomial* result = new Polynomial();
	A->exp_taylor(*result, numVars, order);
	return result;
}
DLLEXPORT Polynomial* CDECL RecPolynomial(Polynomial* A, int numVars, int order)
{
	Polynomial* result = new Polynomial();
	A->rec_taylor(*result, numVars, order);
	return result;
}
DLLEXPORT void CDECL CutoffPolynomial(Polynomial* A)
{
	A->cutoff();
}
DLLEXPORT void CDECL DumpPolynomial (Polynomial* A, int numVars, const char** varNames, bool dumpIntervals)
{
	vector<string> vars;
	for (int i = 0; i < numVars; ++i)
		vars.push_back(varNames[i]);
	if (dumpIntervals)
		A->dump_interval(stdout, vars);
	else
		A->dump_constant(stdout, vars);
}

DLLEXPORT TaylorModel* CDECL CreateTaylorModel (Polynomial* polynomial, Interval* ival)
{
	return new TaylorModel(*polynomial, *ival);
}
DLLEXPORT void CDECL DeleteTaylorModel (TaylorModel* tm)
{
	delete tm;
}

DLLEXPORT void CDECL DumpTaylorModel (TaylorModel* A, int numStateVars, const char** stateVarNames, bool dumpIntervals)
{
	vector<string> stateVars;
	for (int i = 0; i < numStateVars; ++i)
		stateVars.push_back(stateVarNames[i]);
	if (dumpIntervals)
		A->dump_interval(stdout, stateVars);
	else
		A->dump_constant(stdout, stateVars);
}

DLLEXPORT TaylorModelVec* CDECL CreateTaylorModelVec (TaylorModel** tmlist, int numTaylorModels)
{
	vector<TaylorModel> list;
	for (int i = 0; i < numTaylorModels; ++i)
		list.push_back(*tmlist[i]);	
	return new TaylorModelVec(list);
}
DLLEXPORT void CDECL DeleteTaylorModelVec (TaylorModelVec* tmv)
{
	delete tmv;
}
DLLEXPORT void CDECL DumpTaylorModelVec (TaylorModelVec* A, int numStateVars, const char** stateVarNames, int numTmVars, const char** tmVarNames, bool dumpIntervals)
{
	vector<string> stateVars;
	for (int i = 0; i < numStateVars; ++i)
		stateVars.push_back(stateVarNames[i]);
	vector<string> tmVars;
	for (int i = 0; i < numTmVars; ++i)
		tmVars.push_back(tmVarNames[i]);
	if (dumpIntervals)
		A->dump_interval(stdout, stateVars, tmVars);
	else
		A->dump_constant(stdout, stateVars, tmVars);
}

DLLEXPORT vector<HornerForm>* CDECL CreateODE (TaylorModelVec* ode_input)
{
	vector<HornerForm> *hf_output = new vector<HornerForm>();	
	
	int rangeDim = ode_input->tms.size();	
	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm hf;
		Interval I;
		ode_input->tms[i].toHornerForm(hf,I);
		hf_output->push_back(hf);
	} 
	
	/*cout << endl;
	vector<string> varnames;
	varnames.push_back("t");
	varnames.push_back("x1");
	varnames.push_back("x2");
	for (int i = 0; i < hf_output->size(); ++i)
	{
		(*hf_output)[i].dump(stdout, varnames);
		cout << endl;
	}*/
	return hf_output;
}

DLLEXPORT void CDECL DeleteODE (vector<HornerForm>* hfOde)
{
	delete hfOde;
}
DLLEXPORT void CDECL DumpODE (vector<HornerForm>* A, int numVars, const char** varNames)
{
	vector<string> vars;
	for (int i = 0; i < numVars; ++i)
		vars.push_back(varNames[i]);
	for (int i = 0; i < A->size(); ++i)
		(*A)[i].dump(stdout, vars);		
}

extern ContinuousReachability continuousProblem;
DLLEXPORT void CDECL DeclareStateVar(const char* varName)
{
	continuousProblem.declareStateVar(varName);
}

extern double cutoff_threshold;
DLLEXPORT void CDECL SetCutoffThreshold(double threshold)
{
	cutoff_threshold = threshold;
}
