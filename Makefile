-include Makefile.in

AR ?= ar rcu
RANLIB ?= ranlib
CFLAGS ?= -std=c++14 -Wall -O0 -g
CXX ?= $(HOME)/Software/mpich-3.2/bin/mpicxx
H5I ?= -I$(HOME)/Software/hdf5-1.10.1/include
H5L ?= -L$(HOME)/Software/hdf5-1.10.1/lib -lhdf5
SRC := $(wildcard src/*.cpp)
HDR := $(wildcard src/*.hpp)
OBJ := $(SRC:.cpp=.o)
EXE := mara
COW := Cow/src/libcow.a
LUA := Lua/src/liblua.a

default : $(EXE)

$(COW) : .FORCE
	$(MAKE) -C Cow

$(LUA) : .FORCE
	$(MAKE) -C Lua generic

%.o : %.cpp $(HDR)
	$(CXX) $(CFLAGS) -o $@ -c $< -ICow/src -ILua/src

$(EXE) : $(OBJ) $(COW) $(LUA)
	$(CXX) $(CFLAGS) -o $@ $^ $(H5L)

clean :
	$(MAKE) -C Cow clean
	$(MAKE) -C Lua clean
	$(RM) $(EXE) $(OBJ)

.FORCE :
