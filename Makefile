CPPFLAGS = -std=c++11

all: sequential_dfs

sequential_dfs: sequential_dfs.cpp
	g++ -g -o sequential_dfs sequential_dfs.cpp ${CPPFLAGS}

clean:
	-rm -f sequential_dfs
