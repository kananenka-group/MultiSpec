exe     = water_gen
CXX     = g++

SOFTDIR=/work/akanane/sw/MultiSpec/src/software/xdrfile/build
BUILDIR=/work/akanane/sw/MultiSpec/build

INC     = -I$(SOFTDIR)/include/xdrfile
FLAGS   = -O2 -march=native -Wall -std=c++11 -MMD -MP -fopenmp -DUSEOMP
DEBUGFLAGS = -g -std=c++11 -Wall -fopenmp -MMD -MP #-DUSEOMP #-DDEBUG -pg
LIBS    = -lm -ldl -lgomp -lpthread -l boost_program_options -lxdrfile
LIBDIRS = -L$(SOFTDIR)/lib

#MMD and MP generate .d files
OBJ_DIR = OBJ
SRCS  = main.cpp input.cpp ../traj/traj.cpp wmaps.cpp wmapb.cpp water.cpp ../util/util.cpp

OBJS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
DEPS := $(OBJS:%.o=%.d)
OBJS_DEBUG := $(SRCS:%.cpp=$(OBJ_DIR)/%_d.o)
DEPS_DEBUG := $(OBJS:%.o=%_d.d)

exed := "$(exe)_debug"

cpu: $(exe)
debug: $(exed)
all: cpu debug

#link
$(exe): $(OBJS)
	$(CXX) $(FLAGS) -o $(exe) $^ $(LIBDIRS) $(LIBS)
$(exed): $(OBJS_DEBUG)
	$(CXX) $(DEBUGFLAGS) -o $(exed) $^ $(LIBDIRS) $(LIBS)

#compile
$(OBJ_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(FLAGS) -c -o $@ $<
$(OBJ_DIR)/%_d.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(INC) $(DEBUGFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf "$(OBJ_DIR)"

-include $(DEPS)
