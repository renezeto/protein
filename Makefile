#CXXFLAGS = -g -O3
CXXFLAGS = -G

#OBJECTS = grid.h showDensity.

all : protein

clean:
	rm -f *.o protein

#protein: test.o #foo.o

protein: protein.cpp protein.h
	g++ protein.cpp protein.h -o protein -I /usr/include/eigen2 -fpermissive
test.o: test.h

#foo.o: protein.h

.PHONY: clean all

#CXXFLAGS=-g -O3

#OBJECTS=Grid.o foo.o

#all: grid

#clean:
#	rm -f *.o grid

#grid: $(OBJECTS)
#	g++ -o grid $(OBJECTS)


#Grid.o: foo.h
#foo.o: foo.h

#.PHONY: clean all

