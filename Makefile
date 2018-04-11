CPPFLAGS = -std=c++11

all: sequential_dfs parallel_dfs

sequential_dfs: sequential_dfs.cpp
	g++ -g -o sequential_dfs sequential_dfs.cpp ${CPPFLAGS}

parallel_dfs: parallel_dfs_cuda.cpp
	g++ -g -o parallel_dfs parallel_dfs_cuda.cpp ${CPPFLAGS} # need to make change for CUDA

clean:
	-rm -f sequential_dfs parallel_dfs
