
CXX = gcc
CXXFLAGS = -MD -Wall -Os -march=native -ggdb -std=gnu++0x -I /usr/include/eigen2/
LDLIBS = -lstdc++

all: orthomesh
orthomesh: orthomesh.o omesh2d.o

example: example1.svg example2.svg

example1.svg: example1.om orthomesh
	./orthomesh -vvv -o example1.svg -f svg -s 200 -i example1.om

example2.svg: example2.om orthomesh
	./orthomesh -vvv -o example2.svg -f svg -s 200 -i example2.om

test : gjk_test omesh2d_test
	./gjk_test
	./omesh2d_test

gjk_test: gjk_test.o
omesh2d_test: omesh2d.o omesh2d_test.o

clean:
	rm -f *.o *.d orthomesh gjk_test omesh2d_test
	rm -f omesh2d_test.svg omesh2d_test.html example?.svg

-include *.d

