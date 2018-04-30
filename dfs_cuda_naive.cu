#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK_SIZE 1024
using namespace std;

struct compressed_sparse_column {
	int* data;
	int* row;
	int* column;
	int* index_column;
	int* index_row_start;
	int* index_row_end;
};

struct graph {
	compressed_sparse_column* dataset;
	bool* roots;
	bool* leaves;
	bool* singletons;
	int vertices;
	int edges;
};

__host__ graph* read_data(const char* file) {
	ifstream fin(file);
	int rows, columns, nonzeros;
	
	compressed_sparse_column* dataset = new compressed_sparse_column;
	
	while(fin.peek() == '%') {
		fin.ignore(2048, '\n');
	}
	
	fin >> rows >> columns >> nonzeros;
	
	bool find_roots[rows];
	memset(find_roots, false, rows * sizeof(bool));
	bool find_leaves[rows];
	memset(find_leaves, false, rows * sizeof(bool));
	
	dataset->row = new int[nonzeros];
	dataset->column = new int[nonzeros];
	dataset->data = new int[nonzeros];
	
	for(int line = 0; line < nonzeros; line++) {
		int i, j, value;
		fin >> i >> j >> value;
	
		dataset->row[line] = i - 1;
		dataset->column[line] = j - 1;
		dataset->data[line] = value;
	
		find_roots[i - 1] = true; // incoming edges
		find_leaves[j - 1] = true; // outgoing edges
	}

	dataset->index_column = new int[rows];
	
	for(int i = 0; i < rows; i++) {
		int start = -1;
		int row_start = -1;
		int row_end = -1;
		bool found = false;
		for(int j = 0; j < nonzeros; j++) {
			if( dataset->column[j] == i ) {
				start = j;
				break;
			}
			else if( dataset->column[j] > i ) {
				break;
			}

			if( dataset->row[j] == i && !found ) {
				found = true;
				row_start = j;
			}
			else if( dataset->row[j] && found ) {
				row_end = j;
			}
		}
		dataset->index_column[i] = start;
		dataset->index_row_start[i] = row_start;
		dataset->index_row_end[i] = row_end;
	}

	fin.close();
	
	graph* dataset_graph = new graph;
	dataset_graph->dataset = dataset;
	dataset_graph->vertices = rows;
	dataset_graph->edges = nonzeros;
	
	dataset_graph->singletons = new bool[dataset_graph->vertices];
	memset(dataset_graph->singletons, false, dataset_graph->vertices * sizeof(bool));
	dataset_graph->roots = new bool[dataset_graph->vertices];
	memset(dataset_graph->roots, false, dataset_graph->vertices * sizeof(bool));
	dataset_graph->leaves = new bool[dataset_graph->vertices];
	memset(dataset_graph->leaves, false, dataset_graph->vertices * sizeof(bool));
	
	for(int i = 0; i < dataset_graph->vertices; i++) {
		if(find_roots[i] == false && find_leaves[i] == false) {
			dataset_graph->singletons[i] = true;
		}
		else if(find_roots[i] == false) {
			dataset_graph->roots[i] = true;
		}
		else if(find_leaves[i] == false) {
			dataset_graph->leaves[i] = true;
		}
	}
	return dataset_graph;
}

__device__ char* my_strcpy(char* dest, char* src) {
	int i = 0;
	do {
		dest[i] = src[i];
	} while(src[i++] != '\0' );
	return dest;
}

__device__ char* my_strcat(char* dest, char* src) {
	int i = 0;
	while(dest[i] != '\0') {
		i++;
	}
	my_strcpy(dest + i, src);
	return dest;
}

__device__ int my_strcmp(char* a, char* b) {
	int i = 0;
	while( a[i] != '\0' ) {
		i++;
	}
	int a_size = i;
	i = 0;
	while( b[i] != '\0' ) {
		i++;
	}
	int b_size = i;
	if( a_size == b_size ) {
		while( a[i] != '\0' ) {
			if( a[i] > b[i] ) {
				return 1;
			}
			else if( a[i] < b[i] ) {
				return 2;
			}
			i++;
		}
		return 0;
	}
	else {
		return (a_size > b_size) ? 1 : 2;
	}
}

__device__ char* my_itoa(int number, char* str) {
	int i = 0;
	int counter = 0;
	if( number == 0 ) {
		str[i++] = '0';
		str[i] = '\0';
		return str;
	}
	while(number != 0) {
		int remainder = number % 10;
		str[i++] = remainder + '0';
		number = number / 10;
	}
	str[i] = '\0';
	i--;
	char* rev = new char[i];
	
	while(i >= 0) {
		rev[counter++] = str[i];		
		i--;
	}
	
	rev[counter] = '\0';
	
	return rev;
}

__device__ int exclusive_prefix_sum(int* zeta_tilde, int* zeta) {
	int n = sizeof(zeta) / sizeof(zeta[0]);
	zeta_tilde[0] = 0;
	for(int index = 1; index < n; index++){
		zeta_tilde[index] = zeta_tilde[index - 1] + zeta[index - 1];
	}
	return 1;
}

// remove the q list out of here and do this on the CPU.
// for this part we need to implement the exclusive prefix sum too; So see what can be done about the same.

__global__ void calculate_exclusive_prefix_sum(bool* c, int* zeta, int* zeta_tilde, graph* dataset_graph){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( c[i] ){
		for(int j = dataset_graph->dataset->index_column[i]; dataset_graph->dataset->column[j] == i; j++) { // not ordered now!
			int neighbor_vertex = dataset_graph->dataset->row[j];
			zeta[neighbor_vertex] += 1;
		}
		exclusive_prefix_sum(zeta_tilde, zeta);
	}
}

// no need for zeta_tilde in this definition now.
__global__ void subgraph_size(bool* queue, bool* c, int* zeta, graph* dataset_graph) { 
	
	bool* outgoing_edges = new bool[dataset_graph->edges];
	memset(outgoing_edges, false, dataset_graph->edges * sizeof(bool));
	
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(queue[i] == true){
		for(int j = dataset_graph->dataset->index_row_start[i]; j <= dataset_graph->dataset->index_row_end[i]; j++) {
			if( dataset_graph->dataset->row[j] == i ) {
				int neighbor_vertex = dataset_graph->dataset->row[j];
				
				outgoing_edges[j] = true;
				bool flag = true;
				for(int k = 0; k < dataset_graph->edges; k++) {
					if( dataset_graph->dataset->column[k] == neighbor_vertex && !outgoing_edges[k] ) {
						flag = false;
						break;
					}
				}
				if( flag ) {
					c[neighbor_vertex] = true;
				}
			}
		}
	}
}

__global__ void pre_post_order(int* depth, int* zeta, int* zeta_tilde, graph* dataset_graph) {
	int* pre = new int[dataset_graph->vertices];
	int* post = new int[dataset_graph->vertices];

	memset(pre, 0, dataset_graph->vertices * sizeof(int));
	memset(post, 0, dataset_graph->vertices * sizeof(int));

	bool* incoming_edges = new bool[dataset_graph->edges];
	memset(incoming_edges, false, dataset_graph->edges * sizeof(bool));

	bool* q = new bool[dataset_graph->vertices];
	memcpy(q, dataset_graph->roots, sizeof(bool) * dataset_graph->vertices);

	while(true) {
		bool* p = new bool[dataset_graph->vertices];
		memset(p, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		for(int i = 0; i < dataset_graph->vertices; i++) {
			if( q[i] ) {
				int pre_node = 	pre[i];
				int post_node = post[i];

				for(int j = dataset_graph->dataset->index_column[i]; dataset_graph->dataset->column[j] == i; j++) {
					int neighbor_vertex = dataset_graph->dataset->row[j];
					// zeta[i] = undefined!
					pre[neighbor_vertex] = pre_node + zeta_tilde[neighbor_vertex];
					post[neighbor_vertex] = post_node + zeta_tilde[neighbor_vertex];

					incoming_edges[j] = true;
					bool flag = true;
					for(int k = 0; k < dataset_graph->edges; k++) {
						if( dataset_graph->dataset->row[k] == neighbor_vertex && !incoming_edges[k] ) {
							flag = false;
							break;
						}
					}
					if( flag ) {
						global_check = true;
						p[neighbor_vertex] = true;
					}
				}
				pre[i] = pre_node + depth[i];
				post[i] = post_node + (zeta[i] - 1);
			}
		}
		q = p;
		if( !global_check ) {
			break;
		}
	}

}

__global__ void dag_to_dt(bool* queue, bool* p, int* depth, int* parent, char **global_path, graph* dataset_graph) {
	bool* incoming_edges = new bool[dataset_graph->edges];
	memset(incoming_edges, false, dataset_graph->edges * sizeof(bool));

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(queue[i] == true){
		for(int j = dataset_graph->dataset->index_column[i]; dataset_graph->dataset->column[j] == i; j++) {
			int neighbor_vertex = dataset_graph->dataset->row[j];
			
			char* old_copy = new char[1000];
			old_copy = my_strcpy(old_copy, global_path[i]);
			
			char* old_path = global_path[i];
			
			char* buffer = new char[1000];
			buffer = my_itoa(neighbor_vertex, buffer);
			char* new_path = my_strcat(old_copy, buffer);

			if( my_strcmp(old_path, new_path) <= 1 ) {
				global_path[i] = my_strcpy(global_path[i], new_path);
				parent[neighbor_vertex] = i;
				depth[i] += 1;
			}

			incoming_edges[j] = true;
			bool flag = true;
			for(int k = 0; k < dataset_graph->edges; k++) {
				if( dataset_graph->dataset->row[k] == neighbor_vertex && !incoming_edges[k] ) {
					flag = false;
					break;
				}
			}
			if( flag ) {
				p[neighbor_vertex] = true;
			}
		}
	}
}

int main() {

	const char* filename = "./../Dataset/fl2010.mtx";
	// cout << "Please provide the dataset filename : ";
	// cin >> filename;

	graph* dataset_graph = read_data(filename);
	printf("Generated the csc matrix!\n");

	char (*path)[1000];
	int* depth = new int[dataset_graph->vertices];
	int* parent = new int[dataset_graph->vertices];
	int* zeta = new int[dataset_graph->vertices];
	int* zeta_tilde = new int[dataset_graph->vertices];
	for(int i = 0; i < dataset_graph->vertices; i++) {
		path[i][0] = '/';
		path[i][1] = '\0';
	}
	fill(depth, depth + dataset_graph->vertices, 0);
	fill(parent, parent + dataset_graph->vertices, -1);
	fill(zeta, zeta + dataset_graph->vertices, 0);
	fill(zeta_tilde, zeta_tilde + dataset_graph->vertices, 0);

	char (*global_path)[1000];
	int* parent_gpu;
	int* zeta_gpu;
	int* zeta_tilde_gpu;
	int* depth_gpu;
	cudaMalloc((void**)&global_path, (dataset_graph->vertices * 1000) * sizeof(char));
	cudaMalloc((void**)&depth_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&parent_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&zeta_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&zeta_tilde_gpu, dataset_graph->vertices * sizeof(int));


	cudaMemcpy(global_path, path, (dataset_graph->vertices * 1000) * sizeof(char), cudaMemcpyHostToDevice);
	cudaMemcpy(parent_gpu, parent, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(zeta_gpu, zeta, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(zeta_tilde_gpu, zeta, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(depth_gpu, depth, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);
	
	// FINAL ALGORITHM

	// Part 1 - dag to dt (finding parents)
	bool* q = new bool[dataset_graph->vertices];
	bool* p = new bool[dataset_graph->vertices];
	memcpy(q, dataset_graph->roots, sizeof(bool) * dataset_graph->vertices);
	
	bool* Q;
	bool* P;
	cudaMalloc((void**)&Q, dataset_graph->vertices * sizeof(bool));
	cudaMalloc((void**)&P, dataset_graph->vertices * sizeof(bool));

	while(true) {
		cudaMemcpy(Q, q, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
		cudaMemset(P, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		dim3 grid, block;
		block.x = BLOCK_SIZE;	
		grid.x = dataset_graph->vertices/ block.x;
		dag_to_dt<<<grid, block>>>(Q, P, depth_gpu, parent_gpu, global_path, dataset_graph);
		cudaMemcpy(p, P, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);

		for(int i = 0; i <= dataset_graph->vertices; i++){
			if(p[i]){
				global_check = true;
				break;
			}
		}
		if( !global_check ) {
			break;
		}
		q = p;
	}
	// copy back to host what's required. (depth, parent and global_path)

	// Part 2 - subgraph size
	// bool* q = new bool[dataset_graph->vertices];
	bool* c = new bool[dataset_graph->vertices];
	memcpy(q, dataset_graph->leaves, sizeof(bool) * dataset_graph->vertices);
	
	// bool* Q;
	bool* C;
	// cudaMalloc((void**)&Q, dataset_graph->vertices * sizeof(bool));
	cudaMalloc((void**)&C, dataset_graph->vertices * sizeof(bool));

	while(true) {
		cudaMemcpy(Q, q, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
		cudaMemset(C, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		dim3 grid, block;
		block.x = BLOCK_SIZE;	
		grid.x = dataset_graph->vertices/ block.x;
		subgraph_size<<<grid, block>>>(Q, C, zeta_gpu, dataset_graph);

		// since kernel calls are asynchronous and I need to calculate zeta_gpu before launching the following kernel.
		cudaDeviceSynchronize();

		calculate_exclusive_prefix_sum<<<grid, block>>>(C, zeta_gpu, zeta_tilde_gpu, dataset_graph);

		cudaMemcpy(c, C, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);
		// copy back zeta_tilde, zeta (check - if required), zeta_tilde 

		for(int i = 0; i <= dataset_graph->vertices; i++){
			if(c[i]){
				global_check = true;
				break;
			}
		}
		if( !global_check ) {
			break;
		}
		q = c;
	}
	// Part 3 - pre_post_order

	return 0;
}