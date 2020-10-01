CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

LIBTTM_OBJ = qtip4pf.o ps.o electrostatics.o gammq.o smear.o \
             ttm2f.o ttm3f.o ttm4f.o

all: water_interface.o libttm.a

water_interface.o: water_interface.cpp
	$(CXX) $(CXXFLAGS) water_interface.cpp -c

libttm.a: $(LIBTTM_OBJ)
	ar cru libttm.a $(LIBTTM_OBJ) && ranlib libttm.a

clean:
	rm -f libttm.a $(LIBTTM_OBJ) water_interface.o

