/* A Bison parser, made by GNU Bison 2.4.2.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2006, 2009-2010 Free Software
   Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     IDENT = 259,
     STATEVAR = 260,
     TMVAR = 261,
     TM = 262,
     EQ = 263,
     GEQ = 264,
     LEQ = 265,
     ASSIGN = 266,
     END = 267,
     MODE = 268,
     INIT = 269,
     BELONGSTO = 270,
     POLYODE1 = 271,
     POLYODE2 = 272,
     VISUALIZE = 273,
     PARAAGGREG = 274,
     INTAGGREG = 275,
     TMAGGREG = 276,
     OUTPUT = 277,
     CONTINUOUS = 278,
     HYBRID = 279,
     SETTING = 280,
     FIXEDST = 281,
     FIXEDORD = 282,
     ADAPTIVEST = 283,
     ADAPTIVEORD = 284,
     MIN = 285,
     MAX = 286,
     REMEST = 287,
     INTERVAL = 288,
     OCTAGON = 289,
     GRID = 290,
     QRPRECOND = 291,
     IDPRECOND = 292,
     TIME = 293,
     MODES = 294,
     JUMPS = 295,
     INV = 296,
     GUARD = 297,
     RESET = 298,
     START = 299,
     MAXJMPS = 300,
     PRINTON = 301,
     PRINTOFF = 302,
     UNSAFESET = 303,
     CONTINUOUSFLOW = 304,
     HYBRIDFLOW = 305,
     TAYLOR_PICARD = 306,
     TAYLOR_REMAINDER = 307,
     TAYLOR_POLYNOMIAL = 308,
     EXP = 309,
     SIN = 310,
     COS = 311,
     LOG = 312,
     SQRT = 313,
     NPODE_TAYLOR = 314,
     CUTOFF = 315,
     PRECISION = 316,
     GNUPLOT = 317,
     MATLAB = 318,
     COMPUTATIONPATHS = 319,
     uminus = 320
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1685 of yacc.c  */
#line 22 "Flowstar/modelParser.y"

	double dblVal;
	string *identifier;
	vector<Interval> *intVec;
	vector<int> *iVec;
	vector<double> *dVec;
	vector<Monomial> *monoVec;
	vector<Polynomial> *polyVec;
	Monomial *mono;
	Polynomial *poly;
	TaylorModelVec *tmVec;
	Matrix *mat;
	vector<vector<double> > *dVecVec;
	vector<PolynomialConstraint> *vecConstraints;
	ResetMap *resetMap;
	Flowpipe *pFlowpipe;
	TaylorModel *ptm;
	Interval *pint;
	vector<string> *strVec;
	TreeNode *pNode;



/* Line 1685 of yacc.c  */
#line 140 "Flowstar/modelParser.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


