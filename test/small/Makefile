

CXX = g++-12
CXXFLAGS = -Wall -Wno-class-memaccess -Ofast -mcpu=native -std=c++23 -I../../src -I../../eigen-3.4.0

LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
LIBS = -lm

EXE = small_chapman

SRC = $(EXE).cpp
OBJ = $(SRC:.cpp=.o)

all: $(EXE)

$(EXE): $(OBJ)
	$(LD) -o $@ $(LDFLAGS) $^ $(LIBS)

%.o: %.cpp $(HDR)
	$(CXX) -o $@ $(CXXFLAGS) -c $<

clean:
	rm -f $(EXE) $(OBJ)