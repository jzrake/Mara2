-include Makefile.in

AR ?= ar rcu
RANLIB ?= ranlib
CFLAGS ?= -std=c++11 -Wall -O0 -g
CXX ?= $(HOME)/Software/mpich-3.2/bin/mpicxx
H5I ?= -I$(HOME)/Software/hdf5-1.10.1/include
H5L ?= -L$(HOME)/Software/hdf5-1.10.1/lib -lhdf5
SRC := $(wildcard src/*.cpp)
HDR := $(wildcard src/*.hpp)
OBJ := $(SRC:.cpp=.o)
EXE := mara
COW := Cow/libcow.a

default : $(EXE)

$(COW) : .FORCE
	$(MAKE) -C Cow

%.o : %.cpp $(HDR)
	$(CXX) $(CFLAGS) -o $@ -c $< -ICow/src

$(EXE) : $(OBJ) $(COW)
	$(CXX) $(CFLAGS) -o $@ $^ $(H5L)

clean :
	$(MAKE) -C Cow clean
	$(RM) $(EXE) $(OBJ)

.FORCE :
