
CXX = gcc
CXXFLAGS = -MD -Wall -Os -march=native -ggdb -std=gnu++0x -I /usr/include/eigen2/
LDLIBS = -lstdc++

all: orthomesh example.svg

orthomesh: orthomesh.o omesh2d.o

example.svg: example.om orthomesh
	./orthomesh -vvv -o example.svg -f svg -s 200 -i example.om

test : gjk_test omesh2d_test
	./gjk_test
	./omesh2d_test

gjk_test: gjk_test.o
omesh2d_test: omesh2d.o omesh2d_test.o

clean:
	rm -f *.o *.d orthomesh gjk_test omesh2d_test
	rm -f omesh2d_test.svg omesh2d_test.html example.svg

-include *.d

