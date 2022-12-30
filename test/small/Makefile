

CXX = g++
CXXFLAGS = -Wall -Ofast -std=c++2b -I../../src -I../../eigen-3.4.0

LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
LIBS = -lm

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.hpp)
OBJ = $(SRC:.cpp=.o)

EXE = small_chapman

all: $(EXE)

$(EXE): $(OBJ)
	$(LD) -o $@ $(LDFLAGS) $^ $(LIBS)

%.o: %.cpp $(HDR)
	$(CXX) -o $@ $(CXXFLAGS) -c $<

clean:
	rm -f $(EXE) $(OBJ)