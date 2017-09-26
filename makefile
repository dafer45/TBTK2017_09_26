#Uncomment the two following lines if compiling without CUDA support.
#CC = g++
#CFLAGS = -Wall -std=c++11 -fopenmp

#Comment out the two following lines if compiling without CUDA support.
CC = nvcc
CFLAGS = -std=c++11 --compiler-options "-fopenmp" -lcusparse

all:
	$(CC) $(CFLAGS) src/main.cpp -I$(TBTK_dir)/hdf5/hdf5-build/include -L$(TBTK_dir)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTBTK -lblas -llapack -lhdf5 -lhdf5_cpp

clean:
	rm -r build/*

