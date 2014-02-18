FLOWSTAR = Flowstar

vpath %.cpp $(FLOWSTAR)
vpath %.c $(FLOWSTAR)

CXX = g++
GMP_HOME= /usr/local/include
GMP_LIB_HOME = /usr/local/lib 
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -I $(GSL_HOME) -g -Wl,--stack,16777216
LINK_FLAGS = -g -L$(GMP_LIB_HOME) -Wl,--stack,16777216 -static -static-libgcc -static-libstdc++
OBJS = $(patsubst %.o,../obj/%.o,Interval.o Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o Geometry.o Constraints.o Hybrid.o modelParser.o flowstarDll.o expInterval.o expFlowpipe.o expContinuousSystem.o)
GENOBJS = ../obj/lex.yy.o ../obj/modelParser.tab.o

all: ../bin/flowstar.dll

../bin/flowstar.dll: $(OBJS) $(GENOBJS)
	g++ -w $(LINK_FLAGS) -shared -o $@ $^ $(LIBS)
	cp $@ c:\Work\if\sahvy\Src\More\Sahvy\Examples

$(OBJS): ../obj/%.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<
	
../obj/modelParser.tab.o: $(FLOWSTAR)/modelParser.tab.c
	$(CXX) -c $(CFLAGS) -o $@ $<

../obj/lex.yy.o: $(FLOWSTAR)/lex.yy.c
	$(CXX) -c $(CFLAGS) -o $@ $<

$(FLOWSTAR)/modelParser.tab.c: $(FLOWSTAR)/modelParser.y
	bison -o $(FLOWSTAR)/modelParser.tab.c -d -v $(FLOWSTAR)/modelParser.y

$(FLOWSTAR)/lex.yy.c: $(FLOWSTAR)/modelLexer.l $(FLOWSTAR)/modelParser.tab.c
	flex -o $(FLOWSTAR)/lex.yy.c $(FLOWSTAR)/modelLexer.l
	
clean: 
	rm -f $(OBJS) $(GENOBJS) $(FLOWSTAR)/modelParser.tab.c $(FLOWSTAR)/modelParser.tab.h $(FLOWSTAR)/modelParser.output $(FLOWSTAR)/lex.yy.c bin/flowstar.dll
