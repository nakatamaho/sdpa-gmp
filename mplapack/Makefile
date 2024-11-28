# Compiler and flags
CXX = g++
CXXFLAGS = -fopenmp -I.. -I. -I/home/docker/gmpxx_mkII/i/GMP-6.3.0/include
LDFLAGS = -L/home/docker/gmpxx_mkII/i/GMP-6.3.0/lib -lgmp

# Source files and target
SRCS = Rgemm_NN_omp.cpp Rgemm_NN_omp_benchmark.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = Rgemm_NN_omp_benchmark

# Build target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Dependencies (optional, requires GCC/Clang)
depend:
	$(CXX) $(CXXFLAGS) -MM $(SRCS) > .depend

-include .depend
