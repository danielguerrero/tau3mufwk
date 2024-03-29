HDR   = ../interface/
OBJ   = ../obj/
LIB   = ../lib/
BIN   = ../bin/
PRG   = ../test/

ifdef CPP_BOOST_PATH
	BOOSTLIBFLAG = -L$(CPP_BOOST_PATH)/lib
else
	BOOSTLIBFLAG = 
endif

HdrSuf  =  .h
SrcSuf  =  .cc
ObjSuf  =  .o
PrgSuf  =  .cpp
BinSuf  =  .exe
LibSuf  =  .so

HDRS  =  $(wildcard $(HDR)*$(HdrSuf))
SRCS  =  $(wildcard *$(SrcSuf))
_OBJS =  $(patsubst %$(SrcSuf), %$(ObjSuf), $(SRCS))
OBJS  =  $(addprefix $(OBJ),$(_OBJS))
PRGS  =  $(wildcard $(PRG)*$(PrgSuf))

_BINS    =  $(wildcard $(PRG)*$(PrgSuf))
__BINS   =  $(_BINS:$(PrgSuf)=$(BinSuf))
___BINS  =  $(notdir $(__BINS))
BINS	 =  $(addprefix $(BIN),${___BINS})

LINKDEF   =  $(wildcard ${HDR}*LinkDef.h)
DICTHDRS  =  $(patsubst $(LINKDEF),,$(HDRS)) $(LINKDEF)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

ARCHL = -m64
CXX  =  g++
CXXFLAGS  = $(ROOTCFLAGS) -Wall -O -fPIC -I$(HDR)
CPP  =  g++
CPPFLAGS  = $(ROOTCFLAGS) -Wall  $(ARCHL) -I$(HDR)

F    = gfortran
FFLAGS = -c

LD       =  g++
LDFLAGS  =  -rdynamic -shared -O3 $(ARCHL)
SONAME	 =  libtau3muAnalysis.so
SONAME2  =  tau3muAnalysis
SOFLAGS  =  -Wl,-soname,$(SONAME)

#CMSLIBS =  -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lPhysicsToolsKinFitter -lJetMETCorrectionsModules
GLIBS   =  $(ROOTGLIBS) -lm -ldl -rdynamic -lboost_program_options -lGenVector -lFoam -lMinuit -lTMVA -lMLP -lXMLIO -lGpad  -lTreePlayer -lRooFit -lRooFitCore -lRooStats  $(CMSLIBS)

#################################################
#if mac 64
ARCH  =  $(shell root-config --arch)
ifeq ($(ARCH),macosx64)
LibSuf  =  .so

ARCHL = -m64

CPPFLAGS  = -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS) -std=c++14

CXXFLAGS  = -fPIC  -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS) -std=c++14

LDFLAGS  =  -dynamiclib -shared -single_module -undefined dynamic_lookup $(ARCHL)
#SONAME	 =  libPlotter.dylib
# SOFLAGS  =
SONAME	 =  libbbbbAnalysis.so
SONAME2  =  bbbbAnalysis
SOFLAGS  =  
endif
#################################################

.PHONY: all exe test clean exef

all: $(LIB)$(SONAME)

exe: $(BINS) 


test:
	@echo "..............................."
	@echo "HDRS          = $(HDRS)"
	@echo "....."
	@echo "DICTHDRS      = $(DICTHDRS)"
	@echo "....."
	@echo "SRCS          = $(SRCS)"
	@echo "....."
	@echo "OBJS          = $(OBJS)"
	@echo "....."
	@echo "PRGS          = $(PRGS)"
	@echo "....."
	@echo "BINS          = $(BIN)%$(BinSuf): $(PRG)%$(PrgSuf)"
	@echo "....."
	@echo "LIBS/SONAME   = $(LIB)$(SONAME)"
	@echo "....."
	@echo "CPP           = $(CPP)"
	@echo "....."
	@echo "CPPFLAGS      = $(CPPFLAGS)"
	@echo "....."
	@echo "LIB           = $(LIB)"
	@echo "....."
	@echo "BOOSTLIBFLAG  = $(BOOSTLIBFLAG)"
	@echo "....."
	@echo "GLIBS         = $(GLIBS)"
	@echo "..............................."


$(OBJ)%$(ObjSuf): %$(SrcSuf) $(HDRS)
	$(CXX) $< -c $(CXXFLAGS) -o $@  

## 28/12/2017 : LC - I removed the dictionary generation, generates a warning and is probably not needed as I am not reading / saving user classes in the input / output
# $(OBJ)mydict.cc: $(DICTHDRS)
# 	@echo "Generating dictionary for  ..."
# 	rootcling -f $(OBJ)mydict.cc -c -p ${CXXFLAGS} $(DICTHDRS)

# $(LIB)$(SONAME): $(OBJS) $(OBJ)mydict.o
# 	$(LD) $(LDFLAGS) $(OBJS) $(OBJ)mydict.o $(SOFLAGS) $(ROOTGLIBS) -o $(LIB)$(SONAME)
#         @echo "Linking $(SONAME):"

$(LIB)$(SONAME): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(SOFLAGS) $(ROOTGLIBS) -o $(LIB)$(SONAME)
        @echo "Linking $(SONAME):"


$(BIN)%$(BinSuf): $(PRG)%$(PrgSuf) $(HDRS) $(LIB)$(SONAME)
	$(CPP) $<  -l$(SONAME2) $(CPPFLAGS) $(GLIBS) -L$(LIB) $(BOOSTLIBFLAG) -o $@


clean:
	rm -f $(OBJ)*$(ObjSuf)  $(OBJ)mydict* $(BIN)*$(BinSuf) $(LIB)*$(LibSuf) 
