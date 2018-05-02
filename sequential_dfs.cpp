#include <bits/stdc++.h>

using namespace std;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;

int time_dfs = 0;

struct compressed_sparse_column {
	vi data;
	vi row;
	vi column;
	vi index_column;
};

struct graph {
	compressed_sparse_column* dataset;
	int vertices;
	int edges;	
};

graph* read_data(string file) {
	ifstream fin(file);
	compressed_sparse_column* dataset = new compressed_sparse_column;
	int rows, columns, nonzeros;
	while(fin.peek() == '%') {
		fin.ignore(2048, '\n');
	}
	fin >> rows >> columns >> nonzeros;
	for(int line = 0; line < nonzeros; line++) {
		int i, j, value;
		fin >> i >> j >> value;		
		dataset->row.push_back(i - 1);
		dataset->column.push_back(j - 1);		
		dataset->data.push_back(value);
	}

	int start = 0;
	for(int i = 0; i < rows; i++) {
		for(vi::iterator it = dataset->column.begin() + start; it != dataset->column.end(); it += 1) {
			if( *it == i ) {
				start = it - dataset->column.begin();
				break;
			}
		}
		dataset->index_column.push_back(start);
	}
	fin.close();
	assert(dataset->data.size() == nonzeros);
	graph* dataset_graph = new graph;
	dataset_graph->dataset = dataset;
	dataset_graph->vertices = rows;
	dataset_graph->edges = nonzeros;
	return dataset_graph;
}

void DFS_VISIT(graph* dataset_graph, int vertex, int* color, int* parents, int* discovery_time, int* finish_time) {
	time_dfs += 1;
	discovery_time[vertex] = time_dfs;
	color[vertex] = 1;
	for(vi::iterator it = dataset_graph->dataset->column.begin() + dataset_graph->dataset->index_column[vertex]; *it == vertex; it++ ) {
		int neighbor_index = it - dataset_graph->dataset->column.begin();
		int neighbor_vertex = dataset_graph->dataset->row[neighbor_index];
		if( color[neighbor_vertex] == 0 ) {
			parents[neighbor_vertex] = vertex;
			DFS_VISIT(dataset_graph, neighbor_vertex, color, parents, discovery_time, finish_time);
		}
	}
	color[vertex] = 2;
	time_dfs += 1;
	finish_time[vertex] = time_dfs;
}

void DFS(graph* dataset_graph, int* color, int* parents, int* discovery_time, int* finish_time) {
	time_dfs = 0;
	for(int i = 0; i < dataset_graph->vertices; i++) {
		if( color[i] == 0 ) {
			DFS_VISIT(dataset_graph, i, color, parents, discovery_time, finish_time);
		}
	}
}


int main() {

	// string filename = "./../Dataset/fl2010.mtx";
	string filename = "./../Dataset/testcase.mtx";
	// const char* filename = "./../Dataset/Stranke94.mtx";
	// string filename = "./../Dataset/cage3.mtx";

	// cout << "Please provide the dataset filename : ";
	// cin >> filename;
	
	clock_t begin = clock();
	graph* dataset_graph = read_data(filename);
	printf("Generated the csc matrix!\n");


	// 0 -> White, 1 -> Gray, 2 -> Black
	int* color = new int[dataset_graph->vertices];
	// -1 -> NIL
	int* parents = new int[dataset_graph->vertices];
	int* discovery_time = new int[dataset_graph->vertices];
	int* finish_time = new int[dataset_graph->vertices];

	fill(color, color + dataset_graph->vertices, 0);
	fill(parents, parents + dataset_graph->vertices, -1);
	fill(discovery_time, discovery_time + dataset_graph->vertices, -1);
	fill(finish_time, finish_time + dataset_graph->vertices, -1);
	
	DFS(dataset_graph, color, parents, discovery_time, finish_time);

	clock_t end = clock();

	double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	printf("Sequential code took %.3f sec for execution.\n", elapsed_time);

	for(int i = 0; i < dataset_graph->vertices; i++) {
		printf("Discovery of %d - %d\n", i, discovery_time[i]);
		printf("Finish of %d - %d\n", i, finish_time[i]);
	}

	return 0;

}