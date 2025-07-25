# -------------------------------------------------
#  Global compiler options
INCLUDES := -I ../../external/
CXX      := g++
CXXFLAGS := -std=c++20 -O3 -march=native -flto -funroll-loops -Wall $(INCLUDES)

LIBS     := -lpython3.13
PYBIND11_FLAGS := $(shell python3 -m pybind11 --includes)

#  Source file that contains the shared implementation
SRC_LIB := \
    ../../src/morse/_core/morse_sequence/morse_sequence.cpp \
    ../../src/morse/_core/morse_frame/union_find.cpp \
    ../../src/morse/_core/morse_frame/morse_frame_base.cpp \
	../../src/morse/_core/morse_frame/ref_map.cpp \
	../../src/morse/_core/morse_frame/coref_map.cpp 

OBJ_LIB  := $(SRC_LIB:.cpp=.o)   # → morse_sequence.o

# -------------------------------------------------
#  Default target: build every executable
.PHONY: all
all: test_m_seq timer_m_seq test_m_frm test_bind timer_m_frm

# -------------------------------------------------
#  Generic rule: compile any .cpp file to a .o file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PYBIND11_FLAGS) -c $< -o $@

# -------------------------------------------------
#  Executable targets
test_m_seq: test_m_sequence.o $(OBJ_LIB)
	$(CXX) $^ $(LIBS) -o $@

timer_m_seq: timer_m_sequence.o $(OBJ_LIB)
	$(CXX) $^ $(LIBS) -o $@

test_m_frm: test_m_frame.o $(OBJ_LIB)
	$(CXX) $^ $(LIBS) -o $@

test_bind: test_bindings.o $(OBJ_LIB)
	$(CXX) $^ $(LIBS) -o $@


#timer_m_frm : timer_m_frame.o $(OBJ_LIB)
#	$(CXX) $^ $(LIBS) -o $@

# -------------------------------------------------
#  Clean up all build artifacts
.PHONY: clean
clean:
	rm -f */*.o *.o test_m_seq timer_m_seq test_m_frm test_bind test_m_frm
