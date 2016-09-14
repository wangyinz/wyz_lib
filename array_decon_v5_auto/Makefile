BIN=array_decon
PF=array_decon.pf
SUBDIR=/contrib
CXX=g++

ldlibs=-L$(XMOTIFLIB) -lseisppplot -lseisw -lsciplot -lXm -lXt \
	-lgclgrid -lseispp -lperf -L/Users/wangyinz/usr/local/lib -lgsl -lgslcblas -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS) -lseisw -lseisppplot -L$(BOOSTLIB) -lboost_thread -lboost_system

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(XMOTIFINCLUDE)
CXXFLAGS += -I$(BOOSTINCLUDE)
CXXFLAGS += -I/Users/wangyinz/usr/local/include -fopenmp
LDLIBS += -L/sw/lib/gcc4.7/lib 
LDFLAGS += -Wl,-rpath,/Users/wangyinz/usr/local/lib/ 

OBJS= array_decon.o ComplexArray.o DeconOperator.o dpss.o
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
