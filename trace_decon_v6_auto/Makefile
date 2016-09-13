BIN=trace_decon
DATADIR=schemas/css3.0.ext
DATA=decon
SUBDIR=/contrib

PF=trace_decon.pf
CXX=c++

ldlibs=-L$(XMOTIFLIB) -lXm -lXt \
	 -lgclgrid -lseispp -lperf  -lm -llapack -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS) -L$(BOOSTLIB)  

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(XMOTIFINCLUDE)
CXXFLAGS += -I$(BOOSTINCLUDE)
CXXFLAGS += -fopenmp -O3

OBJS= trace_decon.o ComplexArray.o DeconOperator.o dpss.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
