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
AR       ?= ar rcu
RANLIB   ?= ranlib
CXXFLAGS ?= -std=c++14 -Wall -O0 -g
CXX      ?= mpicxx
H5I      ?= -I/usr/include
H5L      ?= -L/usr/lib -lhdf5


# Build macros
# =====================================================================
COW      := Cow/src/libcow.a
SRC      := $(wildcard src/*.cpp) $(wildcard src/Problems/*.cpp)
OBJ      := $(SRC:%.cpp=%.o)
DEP      := $(SRC:%.cpp=%.d)
CXXFLAGS += -MMD -MP
CXXFLAGS += -ICow/src
LDFLAGS  += $(H5L)


# Build rules
# =====================================================================
#
mara: $(OBJ) $(COW)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(DEP) mara

clean-all: clean
	$(MAKE) -C Cow clean

$(COW): FORCE
	$(MAKE) -C Cow

FORCE:

-include $(DEP)
