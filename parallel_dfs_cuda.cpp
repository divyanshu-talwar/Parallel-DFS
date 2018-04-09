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
	vi roots;
	vi leaves;
	vi singletons;
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
	bool find_roots[rows] = {false};
	bool find_leaves[rows] = {false};
	for(int line = 0; line < nonzeros; line++) {
		int i, j, value;
		fin >> i >> j >> value;
		dataset->row.push_back(i - 1);
		dataset->column.push_back(j - 1);
		find_roots[i - 1] = true;
		find_leaves[j - 1] = true;
		dataset->data.push_back(value);
	}

	int start = -1;
	for(int i = 0; i < rows; i++) {
		for(vi::iterator it = dataset->column.begin() + start; it != dataset->column.end(); it++) {
			if( *it == i ) {
				start = it - dataset->column.begin();
				break;
			}
			else if( *it > i ) {
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
	for(int i = 0; i < rows; i++){
		// if a node has no edges associated with it, it is considered a root node.
		if(find_roots[i] == false && find_leaves[i] == false){
			dataset_graph->singletons.push_back(i);
		}
		else if(find_roots[i] == false){
			dataset_graph->roots.push_back(i);
		}
		else if(find_leaves[i] == false){
			dataset_graph->leaves.push_back(i);
		}
	}
	return dataset_graph;
}



int main() {

	string filename = "./../Dataset/fl2010.mtx";
	// cout << "Please provide the dataset filename : ";
	// cin >> filename;

	graph* dataset_graph = read_data(filename);
	printf("Generated the csc matrix!\n");

	cout<<dataset_graph->roots.size()<<endl;
	cout<<dataset_graph->leaves.size()<<endl;
	cout<<dataset_graph->singletons.size()<<endl;
	clock_t begin = clock();

	// // 0 -> White, 1 -> Gray, 2 -> Black
	// int* color = new int[dataset_graph->vertices];
	// // -1 -> NIL
	// int* parents = new int[dataset_graph->vertices];
	// int* discovery_time = new int[dataset_graph->vertices];
	// int* finish_time = new int[dataset_graph->vertices];

	// fill(color, color + dataset_graph->vertices, 0);
	// fill(parents, parents + dataset_graph->vertices, -1);
	// fill(discovery_time, discovery_time + dataset_graph->vertices, -1);
	// fill(finish_time, finish_time + dataset_graph->vertices, -1);
	
	// DFS(dataset_graph, color, parents, discovery_time, finish_time);

	// clock_t end = clock();

	// double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	// printf("Sequential code took %.3f sec for execution.\n", elapsed_time);

	// return 0;

}