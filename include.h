/*---
  Flow*: A Taylor Model Based Flowpipe analyzer.
  Authors: Xin Chen, Erika Abraham and Sriram Sankaranarayanan.
  Email: Xin Chen <xin.chen@cs.rwth-aachen.de> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL). Please consult the file LICENSE.txt for
  further information.
---*/

#ifndef INCLUDE_H_
#define INCLUDE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <mpfr.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <cassert>
#include <map>
#include <time.h>
#include <algorithm>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <glpk.h>

const int intervalNumPrecision = 53;

#define THRESHOLD_HIGH			1e-12
#define THRESHOLD_LOW			1e-20

#define STOP_RATIO				0.99

#define PN 						10 			// the number of digits printed
#define INVALID 				-1e8

#define DC_THRESHOLD_SEARCH 	1e-5
#define DC_THRESHOLD_IMPROV		0.9

#define ID_PRE			0
#define QR_PRE			1

#define UNIFORM			0
#define MULTI			1

#define LAMBDA			0.5
#define NAME_SIZE		50

#define UNSAFE			-1
#define SAFE			1
#define UNKNOWN			0

#define PLOT_INTERVAL	0
#define PLOT_OCTAGON	1
#define PLOT_GRID		2

#define INTERVAL_AGGREG	0
#define PARA_AGGREG		1
#define PCA_AGGREG		2

#define MSG_SIZE		100
#define NUM_LENGTH		50

#define LOW_DEGREE		1
#define HIGH_DEGREE		2
#define NON_POLY		3

#define PI_UP			3.141592654
#define PI_LO			3.141592653

const char outputDir[] = "./outputs/";
const char imageDir[] = "./images/";
const char local_var_name[] = "local_x_";

const char str_prefix_taylor_picard[] = "taylor picard { ";
const char str_prefix_taylor_remainder[] = "taylor remainder { ";
const char str_prefix_taylor_polynomial[] = "taylor polynomial { ";
const char str_prefix_replace[] = "replace { ";

const char str_suffix[] = " }";

extern int lineNum;

using namespace std;

extern void parseODE();

#endif /* INCLUDE_H_ */
