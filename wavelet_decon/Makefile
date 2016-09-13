BIN=wavelet_decon
PF=wavelet_decon.pf
SUBDIR=/contrib

CXX=c++

ldlibs=-L$(XMOTIFLIB) -lXm -lXt \
	 -lgclgrid -lseispp -lperf  -lm -llapack -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS) -L$(BOOSTLIB)  -lgclgrid

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(XMOTIFINCLUDE)
CXXFLAGS += -I$(BOOSTINCLUDE)
CXXFLAGS += -fopenmp -O3

OBJS= wavelet_decon.o ComplexArray.o DeconOperator.o dpss.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
