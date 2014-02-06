#ifndef FLOWSTAR_DLL
#define FLOWSTAR_DLL

#include "Continuous.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define BUILD_DLL
#ifdef BUILD_DLL
#define DLLEXPORT __declspec(dllexport)
#define CDECL __cdecl
#else
#define DLLEXPORT
#define CDECL
#endif

Interval* CreateInterval(double left, double right);
Interval* CreateIntervalCopy(Interval* A);
void DeleteInterval(Interval* A);
double SupInterval(Interval* A);
double InfInterval(Interval* A);
double MidpointInterval(Interval* A);
double WidthInterval(Interval* A);
bool SubseteqInterval(Interval* A, Interval *B);
Interval* AddInterval(Interval* A, Interval *B);
Interval* SubInterval(Interval* A, Interval *B);
Interval* MulInterval(Interval* A, Interval *B);
Interval* DivInterval(Interval* A, Interval *B);
void AddAssignInterval(Interval* A, double c);
void SubAssignInterval(Interval* A, double c);
void MulAssignInterval(Interval* A, double c);
void DivAssignInterval(Interval* A, double c);
Interval* SqrtInterval(Interval* A);
Interval* InvInterval(Interval* A);
Interval* RecInterval(Interval* A);
Interval* SinInterval(Interval* A);
Interval* CosInterval(Interval* A);
Interval* ExpInterval(Interval* A);
void ToStringInterval(Interval* A, int length, char *buffer);

Flowpipe* CreateFlowpipe (int numVars, Interval** domain, Interval* time);
void DeleteFlowpipe (Flowpipe* flowpipe);
void EvalFlowpipe (Flowpipe* flowpipe, Interval** result);
bool AdvanceLowDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec* tmvOde, vector<HornerForm>* hfOde, double step, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveStepLowDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec* tmvOde,  vector<HornerForm>* hfOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveOrderLowDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, TaylorModelVec* tmvOde,  vector<HornerForm>* hfOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation);
bool AdvanceHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveStepHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveOrderHighDegreeFlowpipe (Flowpipe **result, Flowpipe* flowpipe, vector<HornerForm>* hfOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation);
bool AdvanceNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveStepNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, double miniStep, int order, int precondition, int numVars, double *estimation);
bool AdvanceAdaptiveOrderNonPolynomialFlowpipe (Flowpipe **result, Flowpipe* flowpipe, int numOde, const char** strOde, double step, int order, int maxOrder, int precondition, int numVars, double *estimation);
void DumpFlowpipe (Flowpipe* A, int numStateVars, const char** stateVarNames, int numTmVars, const char** tmVarNames);

ContinuousSystem* CreateContinuousSystem (int numVars);
void DumpContinuousSystem(ContinuousSystem* system, int numVars, const char** varNames);
void SetODEContinuousSystem(ContinuousSystem* system, TaylorModelVec* tmvOde);
void SetStrODEContinuousSystem(ContinuousSystem* system, int numOde, const char** strOde);
void SetInitialSetContinuousSystem(ContinuousSystem* system, Flowpipe* initialSet);
Flowpipe* ReachLowDegreeContinuousSystem(ContinuousSystem* system, double step, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** varNames);
Flowpipe* ReachLowDegreeAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames);
Flowpipe* ReachLowDegreeAdaptiveOrderContinuousSystem(ContinuousSystem* system, double step, double time, int order, int maxOrder, int numVars, double* estimation, const char** stateVarNames);
Flowpipe* ReachHighDegreeContinuousSystem(ContinuousSystem* system, double step, double time, int order, int numVars, double* estimation, const char** stateVarNames);
Flowpipe* ReachHighDegreeAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int numVars, double* estimation, const char** stateVarNames);
Flowpipe* ReachHighDegreeAdaptiveOrderContinuousSystem(ContinuousSystem* system, double step, double time, int order, int maxOrder, int numVars, double* estimation, const char** stateVarNames);
Flowpipe* ReachNonPolynomialContinuousSystem(ContinuousSystem* system, double step, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames);
Flowpipe* ReachNonPolynomialAdaptiveStepContinuousSystem(ContinuousSystem* system, double step, double* miniStep, double time, int order, int precondition, int numVars, double* estimation, bool bPrint, const char** stateVarNames);
void DeleteContinuousSystem (ContinuousSystem* system);

Monomial* CreateMonomial (Interval* I, int numVars, int* degs);
Monomial* CreateConstantMonomial (Interval* I, int numVars);
void DeleteMonomial (Monomial* monomial);
Monomial* AddMonomial (Monomial* A, Monomial* B);
Monomial* MulMonomial (Monomial* A, Monomial* B);
void DumpMonomial (Monomial* A, int numVars, const char** varNames, bool dumpIntervals);

Polynomial* CreateEmptyPolynomial ();
Polynomial* CreatePolynomial (Monomial *monomial);
Polynomial* CreateConstantPolynomial (Interval *constant, int numVars);
void DeletePolynomial (Polynomial* polynomial);
void AddAssignPolynomial (Polynomial* polynomial, Monomial* monomial);
Polynomial* AddPolynomial(Polynomial* A, Polynomial* B);
Polynomial* SubPolynomial(Polynomial* A, Polynomial* B);
Polynomial* MulPolynomial(Polynomial* A, Polynomial* B);
Polynomial* NegPolynomial(Polynomial* A);
Polynomial* SinPolynomial(Polynomial* A, int numVars, int order);
Polynomial* CosPolynomial(Polynomial* A, int numVars, int order);
Polynomial* ExpPolynomial(Polynomial* A, int numVars, int order);
Polynomial* RecPolynomial(Polynomial* A, int numVars, int order);
void CutoffPolynomial(Polynomial* A);
void DumpPolynomial (Polynomial* A, int numVars, const char** varNames, bool dumpIntervals);

TaylorModel* CreateTaylorModel (Polynomial* polynomial, Interval* ival);
void DeleteTaylorModel (TaylorModel* tm);
void DumpTaylorModel (TaylorModel* A, int numVars, const char** varNames, bool dumpIntervals);

TaylorModelVec* CreateTaylorModelVec (TaylorModel** tmlist, int numTaylorModels);
void DeleteTaylorModelVec (TaylorModelVec* tmv);
void DumpTaylorModelVec (TaylorModelVec* A, int numStateVars, const char** stateVarNames, int numTmVars, const char** tmVarNames, bool dumpIntervals);

vector<HornerForm>* CreateODE (TaylorModelVec* tmv);
void DeleteODE (vector<HornerForm>* hfOde);
void DumpODE (vector<HornerForm>* A, int numVars, const char** varNames);

void DeclareStateVar(const char* varName);
void SetCutoffThreshold(double threshold);

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif