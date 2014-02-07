/*
 * main.cpp
 *
 *  Created on: Feb 6, 2013
 *      Author: chen
 */

#include "Continuous.h"

int main()
{
/*
	// initial set [1.1,1.4] x [2.25,2.35], Van-der-Pol system
	Interval x(1.25,1.55), y(2.35,2.45), intZero, intOne(1), intMOne(-1), estimation(-1e-4,1e-4);
	vector<Interval> initialBox;
	initialBox.push_back(x);
	initialBox.push_back(y);

	Flowpipe initialSet(initialBox, intZero);

	vector<int> degrees;
	list<Monomial> monomials;

	degrees.push_back(0);
	degrees.push_back(0);
	degrees.push_back(1);
	Monomial m1(intOne, degrees);
	Polynomial p1(m1);
	TaylorModel dxdt(p1, intZero);

	monomials.push_back(m1);

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(1);
	degrees.push_back(0);
	Monomial m2(intMOne, degrees);
	monomials.push_back(m2);

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(2);
	degrees.push_back(1);
	Monomial m3(intMOne, degrees);
	monomials.push_back(m3);
	Polynomial p2(monomials);
	TaylorModel dydt(p2, intZero);

	vector<TaylorModel> tms;
	tms.push_back(dxdt);
	tms.push_back(dydt);
*/




	// initial set [0.8,1] x [0,0.2], Brusselator
	Interval x(0.99,1), y(0,0.01), intZero, intOne(1), intMOne(-1), estimation(-1e-4,1e-4), I1(-3.5), I2(2.5);
	vector<Interval> initialBox;
	initialBox.push_back(x);
	initialBox.push_back(y);

	Flowpipe initialSet(initialBox, intZero);

	vector<int> degrees;
	list<Monomial> monomials;

	degrees.push_back(0);
	degrees.push_back(0);
	degrees.push_back(0);
	Monomial m1(intOne, degrees);
	monomials.push_back(m1);

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(2);
	degrees.push_back(1);
	Monomial m2(intOne, degrees);
	monomials.push_back(m2);

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(1);
	degrees.push_back(0);
	Monomial m3(I1, degrees);
	monomials.push_back(m3);

	Polynomial p1(monomials);
	TaylorModel dxdt(p1, intZero);

	monomials.clear();

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(1);
	degrees.push_back(0);
	Monomial m4(I2, degrees);
	monomials.push_back(m4);

	degrees.clear();
	degrees.push_back(0);
	degrees.push_back(2);
	degrees.push_back(1);
	Monomial m5(intMOne, degrees);
	monomials.push_back(m5);

	Polynomial p2(monomials);
	TaylorModel dydt(p2, intZero);

	vector<TaylorModel> tms;
	tms.push_back(dxdt);
	tms.push_back(dydt);








	TaylorModelVec tmvOde(tms);
	ContinuousSystem system(tmvOde, initialSet);

	int d = 2;
	int v1 = 1, v2 = 2;
	Matrix A;
	Matrix tmpM(8,d);
	A = tmpM;

	A.set(1, 0, v1-1);

	A.set(1, 1, v2-1);

	A.set(-1, 2, v1-1);

	A.set(-1, 3, v2-1);

	A.set(1, 4, v1-1);
	A.set(1, 4, v2-1);

	A.set(1, 5, v1-1);
	A.set(-1, 5, v2-1);

	A.set(-1, 6, v1-1);
	A.set(1, 6, v2-1);

	A.set(-1, 7, v1-1);
	A.set(-1, 7, v2-1);







	vector<int> outputAxes;
	outputAxes.push_back(v1);
	outputAxes.push_back(v2);

	vector<int> orders, maxOrders;
	orders.push_back(4);
	orders.push_back(4);
	maxOrders.push_back(9);
	maxOrders.push_back(9);

	int globalMaxOrder = 9;

	ContinuousReachability problem(system, 0.03, 6, QR_PRE, outputAxes, true, A, UNIFORM, false, false, estimation, 0.001, orders, maxOrders, globalMaxOrder, true);

	string var1("t");
	string var2("x");
	string var3("y");

	problem.declareVariable(var1);
	problem.declareVariable(var2);
	problem.declareVariable(var3);
	sprintf(problem.outputFileName, "program_test");

	char filename[NAME_SIZE+10];
	sprintf(filename, "%s.plt", problem.outputFileName);

	char dumpname[NAME_SIZE+10];
	sprintf(dumpname, "%s.txt", problem.outputFileName);

	clock_t begin, end;
	begin = clock();
	problem.run();
	end = clock();
	printf("time cost for flowpipe construction: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	problem.composition();

	FILE *fp = fopen(filename, "w");
	FILE *fpDump = fopen(dumpname, "w");

	problem.output_2D_GNUPLOT(fp);

	problem.dump(fpDump);

	return 0;
}

