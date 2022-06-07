exe     = excp
CXX     = g++

SOFTDIR=/work/akanane/sw/MultiSpec/src/software

INC     = -I$(SOFTDIR)/xdrfile/build/include/xdrfile -I$(SOFTDIR)/fftw3/build/include
FLAGS   = -O2 -march=native -Wall -std=c++11 -MMD -MP -fopenmp -DUSEOMP -DMKL_ILP64
DEBUGFLAGS = -g -std=c++11 -Wall -fopenmp -MMD -MP #-DUSEOMP #-DDEBUG -pg
LIBS    = -ldl -lgomp -lpthread -l boost_program_options -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lfftw3 -lm
LIBDIRS = -L$(SOFTDIR)/fftw3/build/lib

#MMD and MP generate .d files
OBJ_DIR = OBJ
SRCS  = main.cpp input.cpp exc.cpp ../fft/dofft.cpp

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
