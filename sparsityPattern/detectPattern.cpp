#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <sstream>

#define NCOL 10
/*
void fillArray(float (&x)[n], int n){

	for (int i=0; i<10;i++){
		x[i] = -1;
	}


}
*/

//
void sampleMTXFile(std::string filename){
	
	std::string line;
	std::ifstream mtx_file (filename);
	
	//an integer to keep track of the current line being parsed
	int currentLine = 0; 
	int nrows, ncols, nnz;
	
	double sampleRate = 0.5;
	int numValues;


	static int nnzColIdx[16144];
	static int nnzRowIdx[16144];
	for (int i=0; i<16144; i++){
		nnzColIdx[i] = -1;
		nnzColIdx[i] = -1;
	}
	static int nnzVal[16144];
	int currentIndex = 0;
	int step;

	if (mtx_file.is_open()){
		while (std::getline(mtx_file, line)) {
			
			std::istringstream iss(line);
			
			if (currentLine == 0){
				std::cout << "processing first line" << std::endl;
				currentLine++;
				continue;

			}else if (currentLine == 1){
				std::cout << "processing second line" << std::endl;
				if (!(iss >> nrows >> ncols >> nnz)) {break;}
				
				std::cout << "number of rows is: " << nrows << std::endl;
				std::cout << "number of cols is: " << ncols << std::endl;
				std::cout << "nnz is: " << nnz << std::endl;	
				
				numValues = ceil(sampleRate * nnz);
				step = floor(nnz / numValues);
				currentLine++;
			
			}else { 
				if (currentLine % step == 0){
					
					int rowIndex, colIndex;
					float value;
					if (!(iss >> rowIndex >> colIndex >> value)) {break;}

					if (currentIndex < 16144){

						nnzColIdx[currentIndex] = colIndex;
						nnzRowIdx[currentIndex] = rowIndex; 
						nnzVal[currentIndex] = value;
				
						currentIndex++;
					}
					currentLine++;
				}else {
					currentLine++;
					continue;
				}


			
			}



		}



	}

	std::ofstream outputfile ("plot.py");
	if (outputfile.is_open())
	{
		outputfile << "import matplotlib.pyplot as plt\n";
		
		outputfile << "rows = [";
		for (int i=0; i<16144;i++){
			
			if (i != 16143) {

				outputfile << nnzRowIdx[i] << ", ";
			}
			else {
				outputfile << nnzRowIdx[i] << "]\n";
			}
		}
		outputfile << "cols = [";
		for (int i=0; i<16144;i++){
			
			if (i != 16143) {

				outputfile << nnzColIdx[i] << ", ";
			}
			else {
				outputfile << nnzColIdx[i] << "]\n";
			}
		}
		
		outputfile << "plt.figure(1)\n";
		outputfile << "plt.plot(cols, rows)\n";
		outputfile << "plt.savefig('mtx.png')\n";
	}

}

// helper method to fill a vector with a specific float value
void fillArray(float * array, size_t n, float value){

	for (int i = 0 ; i < n ; i++){

		array[i] = value;

	}

}

// helper method to fill a vector with float random values between r_start and r_end
void generateRandomItems(float * array, size_t n, float r_start, float r_end){
	
	std::random_device rd;

        std::mt19937 gen(rd());
        
	std::uniform_real_distribution<> dist(r_start, r_end);


        for (int i = 0; i < n; i++){
                array[i] = dist(gen);

        }

}

// main method to sample a 2D matrix 
// takes the following parameters:
// 	mat[][NCOL] : 2D matrix with NCOL columns
// 	n_row : number of rows of the matrix
// 	n_col : number of columns in the matrix
// 	row_sampleRate : sampling rate for the rows, eventually decides the number of rows to search
// 	col_sampleRate : sampling rate for the columns, eventually decides the number of nnz vals per row found
void sampleMatrix(float mat[][NCOL], size_t n_row, size_t n_col, float row_sampleRate, float col_sampleRate){
	
	std::cout << "Number of rows: " << n_row << std::endl;
	std::cout << "Number of columns: " << n_col << std::endl;
	std::cout << "Row sampling rate: " << row_sampleRate * 100 << "%" << std::endl;
	std::cout << "Column sampling rate: " << col_sampleRate * 100 << "%" << std::endl;

	int sampledRows = ceil(n_row * row_sampleRate);
	int step = floor(n_row / sampledRows);
	
	// determine the maximum number of nnz values per row to consider
	int nnz = ceil(n_col * col_sampleRate);


	// go over the rows to be sampled one by one
	for (int i = 0; i < n_row ; i+=step){
		
		std::cout << "Sampling row " << i << "...." << std::endl;		
		
		// a vector to capture the number of non-zero indices per row
		int nnzIdx[nnz];

		// initialize all elements in the indices vector to an invalid value (-1)
		for (int f=0; f<nnz; f++){
			nnzIdx[f] = -1;

		}

		// a vector to store the actual non-zero values
		float nnzVal[nnz];

		// start of the binary search-like technique
		int start = 0 ; 
		int end = n_col-1; 
		int lastIdx = 0;
		int checkAt = 0; 

		// FIRST PART: Searching the left side of the vector (row)
		//
		// keep looping while the distance between start and end is not too small
		// or until we have already reached the maximum number of non-zeros
		while ((end - start) > ceil(n_col/10) || lastIdx == nnz-1){

			// first start with the middle value
			checkAt = floor((end + start) / 2 );

			// if it happens to be a nnz, store the index and the value in the corresponding vectors
			if (mat[i][checkAt] != 0){
				if (nnzIdx[lastIdx] == -1){
					nnzIdx[lastIdx] = checkAt;
					nnzVal[lastIdx] = mat[i][checkAt];
					lastIdx++;

				}

			}
			// adjust the end index to be the middle index
			end = checkAt;

		}

		
		// reset the following values
		start = 0; 
		end = n_col-1; 
		checkAt = 0;

		// SECOND PART: Searching the right side of the vector (row)
		//
		while ((end - start) > ceil(n_col/10) || lastIdx == nnz-1){
			checkAt = ceil((end + start) / 2 );
			if (mat[i][checkAt] != 0){
				if (nnzIdx[lastIdx] == -1){
					nnzIdx[lastIdx] = checkAt;
					nnzVal[lastIdx] = mat[i][checkAt];
					lastIdx++;

				}

			}
			// adjust the start to be the middle value
			start = checkAt;

		}
		
		std::cout << "Non-zero sample map" << std::endl;
		std::cout << "(Column index, Non-zero value)" << std::endl;
		for (int k = 0 ; k < nnz; k++){
			if (nnzIdx[k] != -1) {


				std::cout << "(" << nnzIdx[k] << ", " << nnzVal[k] << ")" << std::endl;
		
			}
		}
	}


}

int main(){
	
	sampleMTXFile("fidap004.mtx");
	/*
	float mat[10][10];

	// initialize matrix with random float values
	for (int i = 0 ; i < 10; i++){
		generateRandomItems(mat[i], 10, 1.2, 400.3);

	}


	// print the generated matrix in order to verify the results
	
	std::cout << "Full matrix: " << std::endl;
	for (int i = 0 ; i < 10; i++){
		for (int j = 0 ; j < 10; j++){
			std::cout << mat[i][j] <<  " ";

		}
		std::cout << std::endl;

	}

	//Sample the matrix using a 0.2 sample rate for rows and 0.5 for columns
	sampleMatrix(mat, 10, 10, 0.2, 0.5);
	*/
}
