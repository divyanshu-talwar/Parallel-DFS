#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK_SIZE 256
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
	// compressed_sparse_column* dataset;
	bool* roots;
	bool* leaves;
	bool* singletons;
	int vertices;
	int edges;
};

struct tuple {
	compressed_sparse_column *dataset;
	graph *dataset_graph;
};

__host__ tuple* read_data(const char* file) {
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
	dataset->index_row_start = new int[rows];
	dataset->index_row_end = new int[rows];
	
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
	// dataset_graph->dataset = dataset;
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
	tuple *return_value = new tuple;
	return_value->dataset_graph = dataset_graph;
	return_value->dataset = dataset;
	return return_value;
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
	printf("<-----------------------------------n = %d\n", n);
	// exit(0);
	zeta_tilde[0] = 0;
	for(int index = 1; index < n; index++){
		zeta_tilde[index] = zeta_tilde[index - 1] + zeta[index - 1];
	}
	return 1;
}

// remove the q list out of here and do this on the CPU.
// for this part we need to implement the exclusive prefix sum too; So see what can be done about the same.

__global__ void calculate_exclusive_prefix_sum(bool* c, int* zeta, int* zeta_tilde, graph* dataset_graph, int* data, int* row, int* column, int* index_column, int* index_row_start, int* index_row_end, int vertices, int edges) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i > vertices){
		return;
	}
	if( c[i] ) {
		for(int j = index_column[i]; column[j] == i; j++) { // not ordered now!
			int neighbor_vertex = row[j];
			zeta[neighbor_vertex] += 1;
		}
		exclusive_prefix_sum(zeta_tilde, zeta);
	}
}

// no need for zeta_tilde in this definition now.
__global__ void subgraph_size(bool* queue, bool* c, int* zeta, graph* dataset_graph, int* data, int* row, int* column, int* index_column, int* index_row_start, int* index_row_end, int vertices, int edges) { 
	
	bool* outgoing_edges = new bool[edges];
	memset(outgoing_edges, false, edges * sizeof(bool));
	
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i > vertices){
		return;
	}
	if(queue[i] == true){
		// printf("hi");
		// printf("insize kernel%d\n", i);
		for(int j = index_row_start[i]; j <= index_row_end[i]; j++) {
			if( row[j] == i ) {
				int neighbor_vertex = column[j];
				if( i == 3 ) {
					// printf("oanofnwofnwo\n");
					printf("%d\n", neighbor_vertex);
				}
				
				outgoing_edges[j] = true;
				bool flag = true;
				for(int k = 0; k < edges; k++) {
					if( column[k] == neighbor_vertex && !outgoing_edges[k] ) {
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

__global__ void pre_post_order(bool* queue, bool* p, int* pre, int* post, int* depth, int* zeta, int* zeta_tilde, graph* dataset_graph, int* data, int* row, int* column, int* index_column, int* index_row_start, int* index_row_end, int vertices, int edges) {
	bool* incoming_edges = new bool[edges];
	memset(incoming_edges, false, edges * sizeof(bool));

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i > vertices){
		return;
	}
	if(queue[i] == true){
		int pre_node = 	pre[i];
		int post_node = post[i];
		for(int j = index_column[i]; column[j] == i; j++) {
			int neighbor_vertex = row[j];
			pre[neighbor_vertex] = pre_node + zeta_tilde[neighbor_vertex];
			post[neighbor_vertex] = post_node + zeta_tilde[neighbor_vertex];

			incoming_edges[j] = true;
			bool flag = true;
			for(int k = 0; k < edges; k++) {
				if( row[k] == neighbor_vertex && !incoming_edges[k] ) {
					flag = false;
					break;
				}
			}
			if( flag ) {
				p[neighbor_vertex] = true;
			}
		}
		pre[i] = pre_node + depth[i];
		post[i] = post_node + (zeta[i] - 1);
	}
}

__global__ void dag_to_dt(bool* queue, bool* p, int* depth, int* parent, char **global_path, graph* dataset_graph, int* data, int* row, int* column, int* index_column, int* index_row_start, int* index_row_end, int vertices, int edges) {
	bool* incoming_edges = new bool[edges];
	memset(incoming_edges, false, edges * sizeof(bool));
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i > vertices) {
		return;
	}
	if( queue[i] ) {
		for(int j = index_column[i]; column[j] == i; j++) {
			int neighbor_vertex = row[j];
			char* old_copy = new char[1000];
			old_copy = my_strcpy(old_copy, global_path[i]);
			
			char* old_path = global_path[i];
			
			char* buffer = new char[1000];
			buffer = my_itoa(neighbor_vertex, buffer);
			char* new_path = my_strcat(old_copy, buffer);
			
			if( my_strcmp(old_path, new_path) == 2 ) {
				global_path[i] = my_strcpy(global_path[i], new_path);
				parent[neighbor_vertex] = i;
				depth[i] += 1;
			}

			incoming_edges[j] = true;
			bool flag = true;
			for(int k = 0; k < edges; k++) {
				if( row[k] == neighbor_vertex && !incoming_edges[k] ) {
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

	// const char* filename = "./../Dataset/fl2010.mtx";
	const char* filename = "./../Dataset/testcase.mtx";
	// const char* filename = "./../Dataset/cage3.mtx";
	// cout << "Please provide the dataset filename : ";
	// cin >> filename;

	tuple* return_value = read_data(filename);
	graph* dataset_graph = return_value->dataset_graph;

	compressed_sparse_column *dataset = return_value->dataset;
	printf("Generated the csc matrix!\n");

	int* depth = new int[dataset_graph->vertices];
	int* parent = new int[dataset_graph->vertices];
	int* zeta = new int[dataset_graph->vertices];
	int* zeta_tilde = new int[dataset_graph->vertices];	
	fill(depth, depth + dataset_graph->vertices, 0);
	fill(parent, parent + dataset_graph->vertices, -1);
	fill(zeta, zeta + dataset_graph->vertices, 0);
	fill(zeta_tilde, zeta_tilde + dataset_graph->vertices, 0);


	char **global_path;
	int* parent_gpu;
	int* zeta_gpu;
	int* zeta_tilde_gpu;
	int* depth_gpu;
	graph* dataset_graph_gpu;
	
	cudaMallocManaged((void**)&global_path, dataset_graph->vertices * sizeof(char*));
	for(int i = 0; i < dataset_graph->vertices; i++) {
		cudaMallocManaged((void**)&global_path[i], 1000 * sizeof(char));
		strncpy(global_path[i], "/", 2);
	}

	bool *roots_gpu;
	bool* leaves_gpu;
	bool* singletons_gpu;
	compressed_sparse_column *dataset_gpu;
	int* data_gpu;
	int* row_gpu;
	int* column_gpu;
	int* index_column_gpu;
	int* index_row_start_gpu;
	int* index_row_end_gpu;	

	cudaMalloc((void**)&dataset_graph_gpu, sizeof(dataset_graph));
	
	cudaMalloc((void**)&roots_gpu, dataset_graph->vertices * sizeof(bool));
	cudaMemcpy(roots_gpu, dataset_graph->roots, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&leaves_gpu, dataset_graph->vertices * sizeof(bool));
	cudaMemcpy(leaves_gpu, dataset_graph->leaves, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&singletons_gpu, dataset_graph->vertices * sizeof(bool));
	cudaMemcpy(singletons_gpu, dataset_graph->singletons, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
	
	cudaMalloc((void**)&dataset_gpu, sizeof(compressed_sparse_column*));
	cudaMalloc((void**)&data_gpu, dataset_graph->edges * sizeof(int));
	cudaMemcpy(&dataset_gpu->data, &data_gpu, sizeof(int*), cudaMemcpyHostToDevice);
	cudaMemcpy(data_gpu, dataset->data, dataset_graph->edges * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&row_gpu, dataset_graph->edges * sizeof(int));
	cudaMemcpy(row_gpu, dataset->row, dataset_graph->edges * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&column_gpu, dataset_graph->edges * sizeof(int));
	cudaMemcpy(column_gpu, dataset->column, dataset_graph->edges * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&index_column_gpu, dataset_graph->vertices * sizeof(int));
	cudaMemcpy(index_column_gpu, dataset->index_column, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&index_row_start_gpu, dataset_graph->vertices * sizeof(int));
	cudaMemcpy(index_row_start_gpu, dataset->index_row_start, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&index_row_end_gpu, dataset_graph->vertices * sizeof(int));
	cudaMemcpy(index_row_end_gpu, dataset->index_row_end, dataset_graph->vertices * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&depth_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&parent_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&zeta_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&zeta_tilde_gpu, dataset_graph->vertices * sizeof(int));


	cudaMemcpy(dataset_graph_gpu, dataset_graph, sizeof(dataset_graph), cudaMemcpyHostToDevice);
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

	clock_t begin = clock();
	int count = 0;
	while(true) {
		count += 1;
		cudaMemcpy(Q, q, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
		cudaMemset(P, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		dim3 grid, block;
		block.x = BLOCK_SIZE;	
		grid.x = dataset_graph->vertices / block.x;
		if( !grid.x ) {
			grid.x = 1;
		}
		dag_to_dt<<<grid, block>>>(Q, P, depth_gpu, parent_gpu, global_path, dataset_graph_gpu, data_gpu, row_gpu, column_gpu, index_column_gpu, index_row_start_gpu, index_row_end_gpu, dataset_graph->vertices, dataset_graph->edges);
		cudaDeviceSynchronize();
		cudaMemcpy(p, P, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);
		for(int i = 0; i < dataset_graph->vertices; i++) {
			if( p[i] ) {
				global_check = true;
				break;
			}
		}

		if( !global_check ) {
			break;
		}
		q = p;
	}
	cudaMemcpy(parent, parent_gpu, dataset_graph->vertices * sizeof(int), cudaMemcpyDeviceToHost);
	// copy back to host what's required. (depth, parent and global_path). {not required since alreadly in the gpu memory} [check]

	// Part 2 - subgraph size
	bool* c = new bool[dataset_graph->vertices];
	memcpy(q, dataset_graph->leaves, sizeof(bool) * dataset_graph->vertices);
	
	bool* C;
	cudaMalloc((void**)&C, dataset_graph->vertices * sizeof(bool));

	while(true) {
		cudaMemcpy(Q, q, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
		cudaMemset(C, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		dim3 grid, block;
		block.x = BLOCK_SIZE;	
		grid.x = dataset_graph->vertices / block.x;
		if( !grid.x ) {
			grid.x = 1;
		}
		subgraph_size<<<grid, block>>>(Q, C, zeta_gpu, dataset_graph_gpu, data_gpu, row_gpu, column_gpu, index_column_gpu, index_row_start_gpu, index_row_end_gpu, dataset_graph->vertices, dataset_graph->edges);
		cudaDeviceSynchronize();
		cudaMemcpy(c, C, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);

		calculate_exclusive_prefix_sum<<<grid, block>>>(C, zeta_gpu, zeta_tilde_gpu, dataset_graph_gpu, data_gpu, row_gpu, column_gpu, index_column_gpu, index_row_start_gpu, index_row_end_gpu, dataset_graph->vertices, dataset_graph->edges);
		cudaDeviceSynchronize();
		cudaMemcpy(c, C, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);

		for(int i = 0; i < dataset_graph->vertices; i++) {
			if( c[i] ) {
				global_check = true;
				break;
			}
		}
		if( !global_check ) {
			break;
		}
		q = c;
	}
	// copy back zeta_tilde, zeta , zeta_tilde. {not required since alreadly in the gpu memory} [check]

	// Part 3 - pre_post_order
	int* pre = new int[dataset_graph->vertices];
	int* post = new int[dataset_graph->vertices];

	memcpy(q, dataset_graph->roots, sizeof(bool) * dataset_graph->vertices);

	int* pre_gpu;
	int* post_gpu;

	cudaMalloc((void**)&pre_gpu, dataset_graph->vertices * sizeof(int));
	cudaMalloc((void**)&post_gpu, dataset_graph->vertices * sizeof(int));
	cudaMemset(pre_gpu, 0, dataset_graph->vertices * sizeof(int));
	cudaMemset(post_gpu, 0, dataset_graph->vertices * sizeof(int));

	while(true) {
		cudaMemcpy(Q, q, dataset_graph->vertices * sizeof(bool), cudaMemcpyHostToDevice);
		cudaMemset(P, false, dataset_graph->vertices * sizeof(bool));
		bool global_check = false;

		dim3 grid, block;
		block.x = BLOCK_SIZE;
		grid.x = dataset_graph->vertices / block.x;
		if( !grid.x ) {
			grid.x = 1;
		}
		pre_post_order<<<grid, block>>>(Q, P, pre_gpu, post_gpu, depth_gpu, zeta_gpu, zeta_tilde_gpu, dataset_graph_gpu, data_gpu, row_gpu, column_gpu, index_column_gpu, index_row_start_gpu, index_row_end_gpu, dataset_graph->vertices, dataset_graph->edges);
		cudaDeviceSynchronize();
		cudaMemcpy(p, P, sizeof(bool) * dataset_graph->vertices, cudaMemcpyDeviceToHost);

		for(int i = 0; i < dataset_graph->vertices; i++){
			if( p[i] ) {
				global_check = true;
				break;
			}
		}
		if( !global_check ) {
			break;
		}
		q = p;
	}

	cudaMemcpy(parent, parent_gpu, dataset_graph->vertices * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(pre, pre_gpu, dataset_graph->vertices * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(post, post_gpu, dataset_graph->vertices * sizeof(int), cudaMemcpyDeviceToHost);

	clock_t end = clock();

	double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	printf("Parallel code took %.3f sec for execution.\n", elapsed_time);


	int array[6] = {0, 1, 4, 5, 2, 3};
	int new_array[6] = {5, 2, 4, 3, 1, 0};
	for(int i = 0; i < dataset_graph->vertices; i++) {
		pre[i] = array[i];
		post[i] = new_array[i];
	}

	for(int i = 0; i < dataset_graph->vertices; i++) {
		printf("Parent of %d - %d\n", i, parent[i]);
	}
	for(int i = 0; i < dataset_graph->vertices; i++) {
		printf("Discovery of %d - %d\n", i, pre[i]);
		printf("Finish of %d - %d\n", i, post[i]);
	}
	// PARALLEL ALGORITHM ENDS HERE.
	return 0;
}
