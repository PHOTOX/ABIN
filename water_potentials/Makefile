CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

LIB_OBJ = qtip4pf.o ps.o electrostatics.o gammq.o smear.o \
             ttm2f.o ttm3f.o ttm4f.o

TEST_OBJ = xyz-water-utils.o io-xyz.o test-ttm.o

all: libwater.a

libwater.a: $(LIB_OBJ)
	ar cru libwater.a $(LIB_OBJ) && ranlib libwater.a

test-ttm: libwater.a $(TEST_OBJ)
	$(CXX) -o test-ttm $(TEST_OBJ) libwater.a

test: test-ttm
	./test-ttm hexamers/prism-tschumper-ccsdt.xyz
	./test-ttm hexamers/cage-tschumper-ccsdt.xyz

clean:
	rm -f *.a *.o test-ttm

.PHONY: test
