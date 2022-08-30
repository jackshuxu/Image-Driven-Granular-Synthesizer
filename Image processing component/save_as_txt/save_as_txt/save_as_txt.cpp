#include <iostream>
#include <vector>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<istream>

using namespace std;
int main() {
	double array[7][2];
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 1; j++) {
			cin >> array[i][j];
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 1; j++) {
			cout << " " << array[i][j];
		}
	}
	fstream file;
	file.open("result.txt", ios_base::out);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 1; j++) {
			file << " " << array[i][j];
		}

		}
	cout << endl;
	
	return 0;
}