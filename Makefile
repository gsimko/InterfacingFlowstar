CXX = g++
include makefile.local
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -g -O3 -Wl,--stack,16777216
LINK_FLAGS = -g -L$(GMP_LIB_HOME) -Wl,--stack,16777216
OBJS = Interval.o Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o Geometry.o Constraints.o Hybrid.o

all: flowstar flowstarDLL

flowstar: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

flowstarDLL: $(OBJS) lex.yy.o modelParser.tab.o modelParser.o flowstarDll.o expInterval.o expFlowpipe.o expContinuousSystem.o
	g++ -O3 -w $(LINK_FLAGS) -shared -o flowstar.dll $^ $(LIBS) 
	
%.o: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

modelParser.tab.c: modelParser.y
	bison -d -v modelParser.y

lex.yy.c: modelLexer.l modelParser.tab.c
	flex modelLexer.l

clean: 
	rm -f flowstar *.o *~ modelParser.tab.c modelParser.tab.h modelParser.output lex.yy.c
