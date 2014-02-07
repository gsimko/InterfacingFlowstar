vpath %.cpp flowstar
vpath %.c flowstar

CXX = g++
GMP_HOME= /usr/local/include
GMP_LIB_HOME = /usr/local/lib 
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -I $(GSL_HOME) -g -O3 -Wl,--stack,16777216
LINK_FLAGS = -g -L$(GMP_LIB_HOME) -Wl,--stack,16777216
OBJS = $(patsubst %.o,../obj/%.o,Interval.o Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o Geometry.o Constraints.o Hybrid.o modelParser.o flowstarDll.o expInterval.o expFlowpipe.o expContinuousSystem.o)
GENOBJS = ../obj/lex.yy.o ../obj/modelParser.tab.o

all: ../bin/flowstar.dll

../bin/flowstar.dll: $(OBJS) $(GENOBJS)
	g++ -O3 -w $(LINK_FLAGS) -shared -o $@ $^ $(LIBS)

$(OBJS): ../obj/%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
	
../obj/modelParser.tab.o: flowstar/modelParser.tab.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

../obj/lex.yy.o: flowstar/lex.yy.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

flowstar/modelParser.tab.c: flowstar/modelParser.y
	bison -o flowstar/modelParser.tab.c -d -v flowstar/modelParser.y

flowstar/lex.yy.c: flowstar/modelLexer.l flowstar/modelParser.tab.c
	flex -o flowstar/lex.yy.c flowstar/modelLexer.l
	
clean: 
	rm -f $(OBJS) $(GENOBJS) flowstar/modelParser.tab.c flowstar/modelParser.tab.h flowstar/modelParser.output flowstar/lex.yy.c bin/flowstar.dll
