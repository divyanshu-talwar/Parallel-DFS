CPPFLAGS = -std=c++11

all: sequential_dfs parallel_dfs

sequential_dfs: sequential_dfs.cpp
	g++ -g -o sequential_dfs sequential_dfs.cpp ${CPPFLAGS}

parallel_dfs: parallel_dfs.cu
	nvcc -G -g -o parallel_dfs parallel_dfs.cu

clean:
	-rm -f sequential_dfs parallel_dfs
