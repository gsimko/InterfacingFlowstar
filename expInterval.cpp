#include "flowstarDll.h"
using namespace std;

DLLEXPORT Interval* CDECL CreateInterval (double left, double right)
{	
	return new Interval(left,right);
}

DLLEXPORT Interval* CDECL CreateIntervalCopy (Interval* A)
{	
	return new Interval(*A);
}

DLLEXPORT void CDECL DeleteInterval (Interval* interval)
{	
	delete interval;
}

DLLEXPORT double CDECL SupInterval (Interval* interval)
{
	return interval->sup();
}
DLLEXPORT double CDECL InfInterval (Interval* interval)
{
	return interval->inf();
}
DLLEXPORT double CDECL MidpointInterval (Interval* interval)
{
	return interval->midpoint();
}
DLLEXPORT double CDECL WidthInterval (Interval* interval)
{
	return interval->width();
}
DLLEXPORT bool CDECL SubseteqInterval (Interval* A, Interval* B)
{
	return A->subseteq(*B);
}
DLLEXPORT Interval* CDECL AddInterval (Interval* A, Interval* B)
{
	Interval *res = new Interval(*A);
	*res += *B;
	return res;
}
DLLEXPORT Interval* CDECL SubInterval (Interval* A, Interval* B)
{
	Interval *res = new Interval(*A);
	*res -= *B;
	return res;
}
DLLEXPORT Interval* CDECL MulInterval (Interval* A, Interval* B)
{
	Interval *res = new Interval(*A);
	*res *= *B;
	return res;
}
DLLEXPORT Interval* CDECL DivInterval (Interval* A, Interval* B)
{
	Interval *res = new Interval(*A);
	*res /= *B;
	return res;
}
DLLEXPORT void CDECL AddAssignInterval (Interval* A, double c)
{
	A->add_assign(c);
}
DLLEXPORT void CDECL SubAssignInterval (Interval* A, double c)
{
	A->sub_assign(c);
}
DLLEXPORT void CDECL MulAssignInterval (Interval* A, double c)
{
	A->mul_assign(c);
}
DLLEXPORT void CDECL DivAssignInterval (Interval* A, double c)
{
	A->div_assign(c);
}
DLLEXPORT Interval* CDECL SqrtInterval (Interval* A)
{
	Interval* result = new Interval();
	A->sqrt(*result);
	return result;
}
DLLEXPORT Interval* CDECL InvInterval (Interval* A)
{
	Interval* result = new Interval();
	A->inv(*result);
	return result;
}
DLLEXPORT Interval* CDECL RecInterval (Interval* A)
{
	Interval* result = new Interval();
	A->rec(*result);
	return result;
}
DLLEXPORT Interval* CDECL SinInterval (Interval* A)
{
	Interval* result = new Interval(A->sin());
	return result;
}
DLLEXPORT Interval* CDECL CosInterval (Interval* A)
{
	Interval* result = new Interval(A->cos());
	return result;
}
DLLEXPORT Interval* CDECL ExpInterval (Interval* A)
{
	Interval* result = new Interval(A->exp());
	return result;
}
DLLEXPORT void CDECL ToStringInterval(Interval* A, int length, char *buffer)
{
	string s;
	A->toString(s);
	strncpy(buffer, s.c_str(), length);
	buffer[length-1] = '\0';
}