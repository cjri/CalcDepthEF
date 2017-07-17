CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/gsl/ -mmacosx-version-min=10.11
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas -lcblas -lm
CDEF		= process_sam.o utilities_sam.o alignment.o matchpairs.o call_snps.o likelihood.o optimisation.o io.o ddups.o

cdef: $(CDEF)
	$(CC) $(CC_FLAGS) $(CDEF) -o calcdepthef $(LD_FLAGS)
process_sam.o: process_sam.cpp
	$(CC) $(CC_FLAGS) -c process_sam.cpp
utilities_sam.o: utilities_sam.cpp
	$(CC) $(CC_FLAGS) -c utilities_sam.cpp
alignment.o: alignment.cpp
	$(CC) $(CC_FLAGS) -c alignment.cpp
matchpairs.o: matchpairs.cpp
	$(CC) $(CC_FLAGS) -c matchpairs.cpp
call_snps.o: call_snps.cpp
	$(CC) $(CC_FLAGS) -c call_snps.cpp
likelihood.o: likelihood.cpp
	$(CC) $(CC_FLAGS) -c likelihood.cpp
optimisation.o: optimisation.cpp
	$(CC) $(CC_FLAGS) -c optimisation.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
ddups.o: ddups.cpp
	 $(CC) $(CC_FLAGS) -c ddups.cpp
