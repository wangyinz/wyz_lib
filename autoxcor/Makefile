BIN=autoxcor
PF=autoxcor.pf
SUBDIR=/contrib
CXX=c++

ldlibs=-L$(XMOTIFLIB) -lXm -lXt \
	-lgclgrid -lseispp -lperf -lseispp $(TRLIBS) $(DBLIBS)  \
       $(F77LIBS) $(X11LIBS) -L$(BOOSTLIB) -lboost_thread -lboost_system -lgclgrid

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)

CXXFLAGS += -I$(XMOTIFINCLUDE)
CXXFLAGS += -I$(BOOSTINCLUDE)

OBJS= autoxcor.o 
$(BIN) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
