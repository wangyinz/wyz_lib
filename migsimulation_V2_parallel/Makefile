BIN=migsimulation
PF=migsimulation.pf
#MAN1=migsimulation.1
SUBDIR=/contrib
ldlibs=-lseispp -ltrvltm -lbrttutil $(TRLIBS) $(DBLIBS) $(F77LIBS) \
               -lperf -lgclgrid -lseispp 

ANTELOPEMAKELOCAL = $(ANTELOPE)/contrib/include/antelopemake.local
include $(ANTELOPEMAKE) 
include $(ANTELOPEMAKELOCAL)

OBJS=migsimulation.o SimplePSPrimarySynthetic.o PointSourcePSSynthetic.o SphericalRayPathArray.o vectorcls.o
#OBJS=migsimulation.o SimplePSPrimarySynthetic.o

LDFLAGS += -L$(BOOSTLIB)

#LIBS=-I ~/boost_1_46_1
CXX=mpic++
CXXFLAGS += $(CCOMMON) $(DBG)  $(cflags) $(CEXTRA) -I$(ANTELOPE)/local/include -I$(ANTELOPE)/include -I$(XINCLUDE) $(DCFLAGS)  $(LIBS)

CXXFLAGS +=  -I$(BOOSTINCLUDE) -O3

$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

#clean:
#	-rm ${BIN}
#	-rm ${OBJS}
