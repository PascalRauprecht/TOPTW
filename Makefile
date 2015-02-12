# ---------------------------------------------------------------------
# CPLEX DIRECTORIES 
# ---------------------------------------------------------------------

#cluster (viadrina)
SYSTEM     = x86-64_sles10_4.1
BINDIST = x86-64_sles10_4.1

LIBFORMAT  = static_pic

#laptop, uni pc und cluster
CPLEXDIR       = /cluster/ILOG/CPLEX_Studio125/cplex
CONCERTDIR     = /cluster/ILOG/CPLEX_Studio125/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

#CCOPT = -O -fPIC -fexceptions -DNDEBUG -DIL_STD
CCOPT = -O3 -std=gnu++0x -fPIC -fexceptions -DNDEBUG -DIL_STD #-O3

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread 

# for the boost library
BOOSTDIR = /home/vschmid/boost_1_47_0

all:
	make all_cpp

execute: all
	make execute_cpp

CONCERTINCDIR = $(CONCERTDIR)/include/
CPLEXINCDIR   = $(CPLEXDIR)/include/

DEBUG		  = -g
GPROF      = -p

CCFLAGS = $(CCOPT) -std=c++0x -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -pedantic -Wall # $(DEBUG) 

CPP_EX = TOPTWSTDP

# ---------------------------------------------------------------------
# Compile and link
# ---------------------------------------------------------------------

all_cpp: $(CPP_EX)

execute_cpp: $(CPP_EX)
	TOPTWSTDP

SRCFILES = TOPTWSTDP.cpp \
data.cpp \
dynpro.cpp \
hssa.cpp \
mip.cpp 

OBJFILES = $(SRCFILES:%.cpp=%.o)

$(CPP_EX):  $(OBJFILES)
	$(CCC) $(CCFLAGS) $(OBJFILES)  -o $(CPP_EX) $(CCLNFLAGS) 
	
#$(BOOSTLIB)

%.o:%.cpp
	$(CCC) -c $(CCFLAGS) $< -o $(<:%.cpp=%.o)

clean:
	rm -f *.o
