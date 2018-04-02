#include <bits/stdc++.h>
using namespace std;
typedef pair<int, int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;

struct compressed_sparse_column {
	vi data;
	vi row;
	vi column;
};

compressed_sparse_column* read_data(string file) {
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
		dataset->row.push_back(i);
		dataset->column.push_back(j);
		dataset->data.push_back(value);
	}
	fin.close();
	assert(dataset->data.size() == nonzeros);
	return dataset;
}

int main() {
	string filename = "./../Dataset/fl2010.mtx";
	// cout << "Please provide the dataset filename : ";
	// cin >> filename;
	compressed_sparse_column* dataset = read_data(filename);
	printf("Generated the csc matrix!\n");
	return 0;
}