CXX ?= g++

CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGE_FILE_SOURCE -c -Wall $(shell pkg-config --cflags opencv)
LDFLAGS += -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGE_FILE_SOURCE $(shell pkg-config --libs --static opencv)

Farneback: Farneback.o dense_flow.o
	$(CXX) -o Farneback Farneback.o dense_flow.o $(LDFLAGS)

Farneback.o: Farneback.cpp 
	$(CXX) -c Farneback.cpp $(CXXFLAGS)

dense_flow.o: dense_flow.cpp dense_flow.h
	$(CXX) -c dense_flow.cpp $(CXXFLAGS)

clean: ; rm -f *.o Farneback
