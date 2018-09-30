CXX=g++
CXXFLAGS=-O3 -fopenmp


2d: 2d.o cell.o recorder2.o tissue.o log.o
	$(CXX) -o $@ $(CXXFLAGS) 2d.o cell.o recorder2.o tissue.o log.o

.cc.o:
	$(CXX) -c -MMD -MP $< $(CXXFLAGS)

clean:
	rm *.o
