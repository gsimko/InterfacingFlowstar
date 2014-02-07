CXX = g++
include makefile.local
LIBS = -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(GMP_HOME) -I $(GSL_HOME) -g -O3 -Wl,--stack,16777216
LINK_FLAGS = -g -L$(GMP_LIB_HOME) -Wl,--stack,16777216
OBJS = $(patsubst %.o,../../obj/%.o,Interval.o Matrix.o Monomial.o Polynomial.o TaylorModel.o Continuous.o Geometry.o Constraints.o Hybrid.o modelParser.o)
GENOBJS = $(patsubst %.o,../../obj/%.o,lex.yy.o modelParser.tab.o)

all: ../../bin/flowstar

../../bin/flowstar: $(OBJS) $(GENOBJS)
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

$(OBJS): ../../obj/%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<
	
$(GENOBJS): ../../obj/%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<


modelParser.tab.c: modelParser.y
	bison -d -v modelParser.y

lex.yy.c: modelLexer.l modelParser.tab.c
	flex modelLexer.l

clean: 
	rm -f $(OBJS) $(GENOBJS) modelParser.tab.c modelParser.tab.h modelParser.output lex.yy.c ../../bin/flowstar
