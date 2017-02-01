CXX = g++
CXXFLAGS = -O3 #-Wall -pedantic
LDFLAGS = -lblas -lm

all:
	$(CXX) $(CXXFLAGS) searchprod.cpp -o searchprod $(LDFLAGS)

clean:
	rm -f searchprod
