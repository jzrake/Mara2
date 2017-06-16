# =====================================================================
# Mara build system
# =====================================================================
#
#
# External library dependencies: HDF5, MPI
# Embedded library dependencies: Cow, Lua
#
#
# Notes
# -----
#
# - A useful resource for techniques to process Makefile dependencies:
# www.microhowto.info/howto/automatically_generate_makefile_dependencies.html
#
# - Using -O0 rather than -O3 during development may reduce compilation time
# significantly.


# Build configuration
# =====================================================================
#
# If a Makefile.in exists in this directory, then use it.
#
-include Makefile.in
#
# Any macros that are omitted receive these default values:
LUA_ARCH ?= generic
AR       ?= ar rcu
RANLIB   ?= ranlib
CXXFLAGS ?= -std=c++14 -Wall -O0 -g
CXX      ?= mpicxx
H5I      ?= -I/usr/include
H5L      ?= -L/usr/lib -lhdf5


# Build macros
# =====================================================================
COW      := Cow/src/libcow.a
LUA      := Lua/src/liblua.a
SRC      := $(wildcard src/*.cpp)
OBJ      := $(SRC:%.cpp=%.o)
DEP      := $(SRC:%.cpp=%.d)
CXXFLAGS += -MMD -MP
CXXFLAGS += -ICow/src -ILua/src
LDFLAGS  += $(H5L)


# Build rules
# =====================================================================
#
mara: $(OBJ) $(COW) $(LUA)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJ) $(DEP) mara

clean-all: clean
	$(MAKE) -C Cow clean
	$(MAKE) -C Lua clean

$(COW): FORCE
	$(MAKE) -C Cow

$(LUA): FORCE
	$(MAKE) -C Lua $(LUA_ARCH)

FORCE:

-include $(DEP)
