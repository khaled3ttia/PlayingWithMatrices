// By: Khaled Abdelaal
// khaled.abdelaal@ou.edu
// December 2020

// Utility to load MTX files and find the full sparsity pattern as a plot
#include <exception>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <queue>
#include "search.h"
#include "util.h"

//functions prototypes
//TODO: finish all function prototypes
void loadMTX(int *&, int *&, float *&, int *, int *, int *, std::string);
int flatIndex(int, int, int ,int);

//This function takes a MTX file as an input and loads it in the Coordinate format
// arguments:
//	*&nnzRowIdx : (in-out) an integer pointer reference to the output COO row indices array
//	*&nnzColIdx : (in-out) an integer pointer reference to the output COO col indices array
//	*&nnzVal    : (in-out) a float pointer reference to the output COO float values array
//	*nnz 	    : (in-out) an integer pointer to the number of non-zeros in the matrix
//	*nrows 	    : (in-out) an integer pointer to the number of rows
//	*ncols 	    : (in-out) an integer pointer to the number of columns
//	filename    : (in)     a string representing filename (file path)
void loadMTX(int *&nnzRowIdx, int *&nnzColIdx, float *&nnzVal, int* nnz, int *nrows, int *ncols,  std::string filename){

	std::cout << "Loading file " << filename << "....." <<  std::endl;
	
	
	int numValues = 0;
	int currentLine = 0;
	int currentIndex = 0;
	
	
	std::string line;
	std::ifstream mtx_file (filename);
	

	if (mtx_file.is_open()){
		while (std::getline(mtx_file, line)) {
			std::istringstream iss(line);

			if (currentLine == 0) {

				// First line in MTX foramt is about the format
				// We are only dealing with COO format so just ignore it
				// but increment the count of lines
				currentLine++;
				continue;
			}else if (currentLine == 1) {

				// Second line contains number of rows, cols, and nnz
				if (!(iss >> *nrows >> *ncols >> *nnz)) {break;}
				
				std::cout << "number of rows is: " << *nrows << std::endl;
				std::cout << "number of cols is: " << *ncols << std::endl;
				std::cout << "nnz is: " << *nnz << std::endl;	
				
				// dynamically allocate memory for COO arrays based on nnz				
				nnzColIdx = new int[*nnz];
				nnzRowIdx = new int[*nnz];
				nnzVal = new float[*nnz];

				// initialize all indices to -1
				
				//TOOD: = {-1} ?
				for (int i=0; i<*nnz; i++){
					nnzColIdx[i] = -1;
					nnzColIdx[i] = -1;
				}
	
				currentLine++;

			}else {
				
				// capture info from each line 
				// row index, column index, and nnz value

				int rowIndex, colIndex;
				float value;

				if (!(iss >> rowIndex >> colIndex >> value)) {break;}
				
				if (currentIndex < *nnz){
					
					nnzColIdx[currentIndex] = colIndex - 1;
					nnzRowIdx[currentIndex] = rowIndex - 1; 
					nnzVal[currentIndex] = value;
				
					currentIndex++;
				}
				currentLine++;

			}

		}


	}
	mtx_file.close();	
	
}

// This function plots the complete sparsity pattern of a matrix in the COO format
// Takes the following arguments:
//	nnz : (in) an integer representing the nnz 
//	*&nnzRowIdx	: (in) a pointer reference to the COO row indices array
//	*&nnzColIdx	: (in) a pointer reference to the COO column indices array
//	filename	: (in) a string representing the matrix file name
void fullSparsity(int nnz, int nrows, int ncols, int *&nnzRowIdx, int *&nnzColIdx, std::string filename){
		
	// create a Python script to plot 
	std::ofstream outputfile ("plot.py");
	if (outputfile.is_open())
	{
		outputfile << "import matplotlib.pyplot as plt\n";
		
		outputfile << "rows = [";
		for (int i=0; i<nnz;i++){
			
			if (i != nnz-1) {

				outputfile << nnzRowIdx[i] << ", ";
			}
			else {
				outputfile << nnzRowIdx[i] << "]\n";
			}
		}
		outputfile << "cols = [";
		for (int i=0; i<nnz;i++){
			
			if (i != nnz-1) {

				outputfile << nnzColIdx[i] << ", ";
			}
			else {
				outputfile << nnzColIdx[i] << "]\n";
			}
		}	
		outputfile << "plt.figure(1)\n";
		outputfile << "ax = plt.gca()\n";
		
		// put the x-axis on the top
		outputfile << "ax.xaxis.tick_top()\n";
		outputfile << "ax.xaxis.set_label_position('top')\n";
		
		// invert the y-axis
		outputfile << "ax.yaxis.tick_left()\n";
		outputfile << "ax.invert_yaxis()\n";
		outputfile << "plt.title('Sparsity Pattern for " + filename + "', y=-0.1)\n";
		outputfile << "plt.xlabel('column index')\n";
		//TODO rows and column indices on figure ?
		//outputfile << "plt.xticks([x for x in range(0," << ncols << ")])\n";
		//outputfile << "plt.yticks([x for x in range(0," << nrows << ")])\n";
		outputfile << "plt.ylabel('row index')\n";
		outputfile << "plt.plot(cols, rows, 's')\n";
		outputfile << "plt.savefig('sp_" + filename +   ".png')\n";
	}
	outputfile.close();

	// execute the script to generate the figure
	// NOTE: needs python3 installation
	system("python plot.py");
	
	std::cout << "Sparsity pattern figure saved at sp_" << filename << ".png !" << std::endl;
	

}

void cooToDense(int *&nnzRowIdx, int *&nnzColIdx, float *&nnzVal, int nnz, int nrows, int ncols, float *&denseMatrix){
	
		
	//denseMatrix = {0};
	for (int i = 0 ; i < nnz; i++){
		int row = nnzRowIdx[i];
		int col = nnzColIdx[i];
		int idx = flatIndex(row, col, nrows, ncols);
		denseMatrix[idx] = nnzVal[i];

	}


}

//TODO: debug
void sampleSparsity(float sampleRate, int *&nnzRowIdx, int *&nnzColIdx, int nnz, int nrows, int ncols, std::string filename){

	std::vector<int> nnzColVec(nnzColIdx, nnzColIdx + nnz);
	std::vector<int> nnzRowVec(nnzRowIdx, nnzRowIdx + nnz);

	int *sampleRowIdx, *sampleColIdx;
	sampleRowIdx = new int[nnz];
	sampleColIdx = new int[nnz];
	int currentIndex = 0; 
	
	int centerIdx = (int)ceil(nnz/2);
	std::cout << "centerIdx for nnzColIdx is: " << centerIdx << std::endl;

	int center_colIdx = nnzColIdx[centerIdx];
	int center_rowIdx = nnzRowIdx[centerIdx];
	
	std::queue<int> centerCols;
	centerCols.push(center_colIdx);



	sampleRowIdx[currentIndex] = center_rowIdx; 
	sampleColIdx[currentIndex] = center_colIdx;
	std::cout << "At centerIdx, we have:" << std::endl;
	std::cout << "Should have : (" << center_rowIdx << "," << center_colIdx << ")" << std::endl;

	//delete from original arrays
	//nnzRowIdx[centerIdx] = -1;
	//nnzColIdx[centerIdx] = -1;	
	nnzColVec.erase(nnzColVec.begin() + centerIdx);
	nnzRowVec.erase(nnzRowVec.begin() + centerIdx);	

	currentIndex++;

	//std::cout << "c_colIdx: " << c_colIdx << std::endl;
	//std::cout << "c_rowIdx: " << c_rowIdx << std::endl;

	int hstep = (int)ceil(nnz/ncols);
	/*
	auto col1_idx = std::lower_bound(nnzColIdx + centerIdx + 1 , nnzColIdx + nnz , c_colIdx + hstep);
	auto col2_idx = std::lower_bound(nnzColIdx , nnzColIdx + centerIdx, c_colIdx - hstep);
	auto col3_idx = std::lower_bound(nnzColIdx + centerIdx + 1 , nnzColIdx + nnz, c_colIdx + (int)ceil(hstep/2));
	auto col4_idx = std::lower_bound(nnzColIdx , nnzColIdx + centerIdx, c_colIdx - (int)ceil(hstep/2));
	*/
	//while (currentIndex <= (int)ceil(0.75 * nnz)){
	int expectedNNZ = (int)ceil(sampleRate * nnz);
	std::cout << "Expected number of NNZ elements to be sampled is : " << expectedNNZ << std::endl;
	while (currentIndex <= expectedNNZ){
	//while (currentIndex <= 5){
		//TODO make sure to fill the queue with items as long as we haven't hit the expectedNNZ threshold yet
		if (centerCols.empty()){

			break;
		}	
		int outerNNZ = 0;
		int innerNNZ = 0;

		int c_colIdx = centerCols.front();
		centerCols.pop();

		int col1 = c_colIdx + hstep;
		int col2 = c_colIdx - hstep;
		int col3 = c_colIdx + (int)ceil(hstep/2);
		int col4 = c_colIdx - (int)ceil(hstep/2);

		//Temporary
		//int row1, row2, row3, row4;

		
		//int col1_idx = binarySearch(nnzColIdx, centerIdx + 1 , nnz-1, col1);
		if (col1 >=0 && col1 < ncols){		
			int col1_idx = binarySearch(nnzColVec, centerIdx + 1 , nnz-1, col1);

			if (col1_idx != -1){
				int row1 = nnzRowIdx[col1_idx];
				sampleRowIdx[currentIndex] = row1;
				sampleColIdx[currentIndex] = col1;
				centerCols.push(col1);
				std::cout << "Should have: (" << row1 << "," << col1 << ")" << std::endl;

				currentIndex++;

				//delete from original arrays
				//nnzColIdx[col1_idx] = -1;
				//nnzRowIdx[col1_idx] = -1;
				nnzColVec.erase(nnzColVec.begin() + col1_idx);
				nnzRowVec.erase(nnzRowVec.begin() + col1_idx);
				//This variable will be local to the current iteration only
				outerNNZ++;
			}
		}

		
		//int col2_idx = binarySearch(nnzColIdx, 0, centerIdx - 1, col2);
		if (col2 >= 0 && col2 < ncols){

			int col2_idx = binarySearch(nnzColVec, 0, centerIdx - 1, col2);
			if (col2_idx != -1){
				int row2 = nnzRowIdx[col2_idx];
				sampleRowIdx[currentIndex] = row2; 
				sampleColIdx[currentIndex] = col2;
				centerCols.push(col2);

				std::cout << "Should have: (" << row2 << "," << col2 << ")" << std::endl;

				currentIndex++;

				//delete from original arrays
				//nnzColIdx[col2_idx] = -1;
				//nnzRowIdx[col2_idx] = -1;
				nnzColVec.erase(nnzColVec.begin() + col2_idx);
				nnzRowVec.erase(nnzRowVec.begin() + col2_idx);

				outerNNZ++;	
			
			}
		}		
		//std::cout << "COL[1840] right now is : " << nnzColIdx[1840] << std::endl;
		//int col3_idx = binarySearch(nnzColIdx, centerIdx + 1  , nnz-1, col3);
		if (col3 >= 0 && col3 < ncols){

			int col3_idx = binarySearch(nnzColVec, centerIdx + 1 , nnz-1, col3);
			if (col3_idx != -1){
				int row3 = nnzRowIdx[col3_idx];
				sampleRowIdx[currentIndex] = row3;
				sampleColIdx[currentIndex] = col3;
				centerCols.push(col3);
				
				std::cout << "Should have: (" << row3 << "," << col3 << ")" << std::endl;

				currentIndex++;
				
				//delete from original arrays
				//nnzColIdx[col3_idx] = -1;
				//nnzRowIdx[col3_idx] = -1;
				nnzColVec.erase(nnzColVec.begin() + col3_idx); 
				nnzRowVec.erase(nnzRowVec.begin() + col3_idx);

				innerNNZ++;

			}
		}
		//int col4_idx = binarySearch(nnzColIdx, 0, centerIdx - 1, col4);
		if (col4 >=0 && col4 < ncols){

			int col4_idx = binarySearch(nnzColVec, 0, centerIdx - 1, col4);
			if (col4_idx != -1){
				int row4 = nnzRowIdx[col4_idx];
				sampleRowIdx[currentIndex] = row4;
				sampleColIdx[currentIndex] = col4;
				centerCols.push(col4);			

				std::cout << "Should have: (" << row4 << "," << col4 << ")" << std::endl;

				currentIndex++;

				//delete from original arrays
				//nnzColIdx[col4_idx] = -1;
				//nnzRowIdx[col4_idx] = -1;
				nnzColVec.erase(nnzColVec.begin() + col4_idx);
				nnzRowVec.erase(nnzRowVec.begin() + col4_idx);

				innerNNZ++;

			}
		}
	}
		std::cout << "Number of sampled NNZ elements: " << currentIndex << std::endl;
		std::cout << "Total number of NNZ : " << nnz << std::endl;
		fullSparsity(currentIndex, nrows, ncols, sampleRowIdx, sampleColIdx, filename + "_sampled" );


}


void fullDense(float *&denseMatrix, int *&rows, int *&cols, int nrows, int ncols, std::string filename){
	
	int currentIndex = 0; 
	for (int i = 0; i < nrows; i++){
		for (int j = 0; j<ncols ; j++){
			if (denseMatrix[i*ncols + j] != 0){
				rows[currentIndex] = i;
				cols[currentIndex] = j;
				currentIndex++;

			}
		}
	}
	

	fullSparsity(currentIndex, nrows, ncols, rows, cols, filename);
	
}

int flatIndex(int i, int j, int nrows, int ncols){
	if (i >= 0 && i < nrows){
		if ( j >=0 && j < ncols){
			return i * ncols + j;
		}
	}

	return -1;
}

void unpackIndex(int index, int ncols, int& row, int& col){
	row = (int)(index/ncols);
	col = (int)(index%ncols);

}

//TODO: debug
void sampleDense(float sampleRate, float *&denseMatrix, int nrows, int ncols, std::string filename){
	
	std::queue<int> checkQ; 
	
	
	int numElements = nrows * ncols; 
	
	int *flatIndices = new int[(int)ceil(numElements*sampleRate)];	

	int centerIdx = (int)floor(numElements/2);
	
	int currentIndex = 0;
	if (denseMatrix[centerIdx] != 0) {

		flatIndices[currentIndex] = centerIdx;
		currentIndex++;
	}
	
	checkQ.push(centerIdx);

	
	int hstep = 2;
	int vstep = 2;	
	while (currentIndex < (int)ceil(numElements*sampleRate)){
		
		if (checkQ.empty()){
			break;
		}		

		int c0 = checkQ.front();
		checkQ.pop();
		int rowIdx, colIdx; 
		unpackIndex(c0, ncols, rowIdx, colIdx);

		int neighbors[8];

		neighbors[0] = flatIndex(rowIdx, colIdx + hstep, nrows, ncols);
		neighbors[1] = flatIndex(rowIdx, colIdx - hstep, nrows, ncols);
		neighbors[2] = flatIndex(rowIdx - vstep, colIdx, nrows, ncols);
		neighbors[3] = flatIndex(rowIdx + vstep, colIdx, nrows, ncols);
		neighbors[4] = flatIndex(rowIdx + (int)ceil(vstep/2) , colIdx + (int)ceil(hstep/2), nrows, ncols);
		neighbors[5] = flatIndex(rowIdx + (int)ceil(vstep/2) , colIdx - (int)ceil(hstep/2), nrows, ncols);
		neighbors[6] = flatIndex(rowIdx - (int)ceil(vstep/2) , colIdx + (int)ceil(hstep/2), nrows, ncols);
		neighbors[7] = flatIndex(rowIdx - (int)ceil(vstep/2) , colIdx - (int)ceil(hstep/2), nrows, ncols);
		
		int x = 0;
		
		std::cout << "Neighbors array: " << std::endl;
		printVector(neighbors, 8);		
	

		while (x < 8){
			if (neighbors[x] != -1){
				if (denseMatrix[neighbors[x]] != 0){
					flatIndices[currentIndex] = neighbors[x];
					denseMatrix[neighbors[x]] = 0;
					currentIndex++;
				}
				std::cout << "Queue size now is : " << checkQ.size() << std::endl;
				checkQ.push(neighbors[x]);
			
			}
			++x;
		}
	}
	printVector(flatIndices, currentIndex);	


}

// A function to print out how to use the tool
void usage(){
	std::cout << "Options:" << std::endl;
	std::cout << "-i <inputfile> : filename of a MTX file" << std::endl;
	std::cout << "-r <sampleRate> : a floating point number (e.g. 0.82) sepecifying sampling rate" << std::endl;
	std::cout << "-d : to specify a dense matrix" << std::endl;
	std::cout << "-g : to generate a dense matrix, instead of specifying an existing file" << std::endl;
	std::cout << "-n <number of rows> : (used only with -g) to specify the number of rows of the generated matrix" << std::endl;	
	std::cout << "-m <number of columns> : (used only with -g) to specify the number of columns of the generated matrix" << std::endl;	
	std::cout << "-z <non-zero rate> : (used only with -g) to specify the rate of non-zero elements of the generated matrix" << std::endl;	
	std::cout << "-h : prints usage information" << std::endl;

	std::cout << "Examples: " << std::endl;
	std::cout << "./dpd.exe -i fidap004.mtx -d -r 0.9" << std::endl;
	std::cout << "loads a matrix from filename fidap004.mtx in a dense format, and samples it with a sampling rate of 90%" << std::endl;
	std::cout << "==================================================================" << std::endl;
	std::cout << "./dpd.exe -g -n 10 -m 12 -z 0.2" << std::endl;
	std::cout << "generates a random matrix of 10 rows and 12 columns, having 20\% of the elements as non-zeros, and find the full sparsity pattern of the matrix" << std::endl; 
	exit(0);

}

// The main function

// TODO: clean code
int main(int argc, char** argv){

	int opt;
	float sampleRate;


	std::string filename; 

	bool sample = false;
	bool loadFile = false;
	bool generate = false;
	bool dcolsSet = false; 
	bool drowsSet = false; 
	bool nnzRateSet = false; 
	bool dense = false;

	int dcols, drows;
	float nnzRate;	

	// TODO re-adjust command line options 
	// parse command line arguments
	while ((opt = getopt(argc, argv, "i:r:n:m:z:hgd")) != -1){
		switch (opt) {
			case 'i': 
					filename = optarg;
				  	loadFile = true;
					 break;
			case 'r':
				 try{
					sampleRate = atof(optarg);
					sample = true;
				 } catch (...){
					std::cout << "Invalid value for sample rate" << std::endl;
					usage();
				
				 }
				break;
			case 'g': generate = true; break;
			case 'd': dense = true; break;
			case 'n': 
				try {
					drows = atoi(optarg);
					drowsSet = true;
					
				}catch (...){

					std::cout << "Invalid value for number of rows to generate matrix" << std::endl;
					usage();
				}
				break;
			case 'm': 
				try {
					dcols = atoi(optarg);
					dcolsSet = true;
					
				}catch (...){

					std::cout << "Invalid value for number of cols to generate matrix" << std::endl;
					usage();
				}
				break;
			case 'z': 
				try {
					nnzRate = atof(optarg);
					nnzRateSet = true;
					
				}catch (...){

					std::cout << "Invalid value for number of non-zero rate to generate matrix" << std::endl;
					usage();
				}
				break;


			case 'h': usage(); break;

		}

	}	

	// For dense:
	// - Load MTX , convert into dense, and sample
	bool sampleLoadDense = loadFile && !generate && dense && sample;

	// - Load MTX, convert into dense , find full pattern
	bool fullLoadDense = loadFile && !generate && dense && !sample;

	// - Generate a dense matrix (given number of rows and columns) and sample
	bool sampleGenDense = generate && !loadFile && drowsSet && dcolsSet && nnzRateSet && sample;

	// - Generate a dense matrix (given number of rows and columns) and find full pattern
	bool fullGenDense = generate && !loadFile && drowsSet && dcolsSet && nnzRateSet && ~sample;

	// For sparse

	// - Load MTX, convert into COO, and find full pattern

	bool fullSparse = loadFile && !generate &&  (!dense) && (!sample);

	// - Load MTX, convert into COO and sample
	bool sampleSparse = loadFile && !generate &&  (!dense) && (sample);

	if (sampleLoadDense) {
		std::cout << "You chose to load a dense matrix from a MTX file and this sample it" << std::endl;
			
		//COO arrays 
		int* nnzRowIdx;
		int* nnzColIdx;
		float* nnzVal;
		// --end of COO arrays
		
		int n, m, nnz; 

		//Load as COO
		auto load_start_time = std::chrono::high_resolution_clock::now();
		loadMTX(nnzRowIdx, nnzColIdx, nnzVal, &nnz, &n, &m, filename);
		auto load_end_time = std::chrono::high_resolution_clock::now();
		
		//Convert COO to Dense
		std::cout << "Converting from COO to Dense format...." << std::endl;
		float *denseMatrix = new float[n*m]; 
		auto convert_start = std::chrono::high_resolution_clock::now();
		cooToDense(nnzRowIdx, nnzColIdx, nnzVal, nnz, n, m, denseMatrix);
		auto convert_end = std::chrono::high_resolution_clock::now();
		
		//Sample Dense 
		std::cout << "Sampling the dense matrix at a sampling rate of: " << sampleRate * 100 << "%" << std::endl;
		int *rowsVec = new int[(int)ceil(nnz*sampleRate)];
		int *colsVec = new int[(int)ceil(nnz*sampleRate)];
		auto sample_start = std::chrono::high_resolution_clock::now();
		sampleDense(sampleRate, denseMatrix, n, m, filename);
		auto sample_end = std::chrono::high_resolution_clock::now();


		auto load_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end_time - load_start_time).count();
		auto convert_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(convert_end - convert_start).count();
		auto sample_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(sample_end-sample_start).count();
		
		std::cout << "File loading duration: " << load_duration << " nanoseconds\n";
		std::cout << "COO to Dense converstion duration: " << convert_duration << " nanoseconds\n";
		std::cout << "Sampling duration: " << sample_duration << " nanoseconds\n";	
	}

	if (fullLoadDense) {
		std::cout << "You chose to load a dense matrix from a MTX file and find its full sparsity (no sampling)" << std::endl;
		
		//COO arrays 
		int* nnzRowIdx;
		int* nnzColIdx;
		float* nnzVal;
		// --end of COO arrays
		
		int n, m, nnz; 

		//Load 
		auto load_start_time = std::chrono::high_resolution_clock::now();
		loadMTX(nnzRowIdx, nnzColIdx, nnzVal, &nnz, &n, &m, filename);
		auto load_end_time = std::chrono::high_resolution_clock::now();

		//Convert COO to Dense 
		std::cout << "Converting from COO to Dense Format...." << std::endl;
		float *denseMatrix = new float[n*m]; 
		auto convert_start = std::chrono::high_resolution_clock::now();
		cooToDense(nnzRowIdx, nnzColIdx, nnzVal, nnz, n, m, denseMatrix);
		auto convert_end = std::chrono::high_resolution_clock::now();
		
		//Find full sparsity pattern
		int *rowsVec = new int[nnz];
		int *colsVec = new int[nnz];
		auto full_start = std::chrono::high_resolution_clock::now();
		fullDense(denseMatrix, rowsVec, colsVec, n, m, filename);
		auto full_end = std::chrono::high_resolution_clock::now();
	
		auto load_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end_time - load_start_time).count();
		auto convert_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(convert_end - convert_start).count();
		auto full_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(full_end-full_start).count();
		
		std::cout << "File loading duration: " << load_duration << " nanoseconds\n";
		std::cout << "COO to Dense converstion duration: " << convert_duration << " nanoseconds\n";
		std::cout << "Full pattern detection duration: " << full_duration << " nanoseconds\n";

	
	}

	if (sampleGenDense) {

		std::cout << "You chose to generate a dense matrix, and sample it" << std::endl;

		std::cout << "Generating a random matrix of " << drows << " rows and " << dcols << " columns" << std::endl;
		std::cout << "with a non-zero rate of " << nnzRate * 100 << "%" << std::endl;

		//Generate a matrix
		float *denseMatrix = new float[drows*dcols]; 
		auto gen_start = std::chrono::high_resolution_clock::now();
		generateMatrix(drows, dcols, nnzRate, denseMatrix);
		auto gen_end = std::chrono::high_resolution_clock::now();

		std::cout << "Generated Dense Matrix: " << std::endl;
		
		//TODO :: print as a 2D matrix ?? 
		printVector(denseMatrix, drows*dcols);	

		
		int *rowsVec = new int[(int)ceil(sampleRate*nnzRate*drows*dcols)];
		int *colsVec = new int[(int)ceil(sampleRate*nnzRate*drows*dcols)];
		
		//Sample matrix
		std::cout << "Sampling matrix with a sampling rate of: " << sampleRate * 100 << "%" << std::endl;
		auto sample_start = std::chrono::high_resolution_clock::now();
		sampleDense(sampleRate, denseMatrix, drows, dcols, "Dense_mat");
		auto sample_end = std::chrono::high_resolution_clock::now();

		auto gen_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(gen_end - gen_start).count();
		auto sample_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(sample_end-sample_start).count();
		
		std::cout << "Matrix generation duration: " << gen_duration << " nanoseconds\n";
		std::cout << "Sampling duration: " << sample_duration << " nanoseconds\n";

	
	}

	if (fullGenDense) {

		std::cout << "You chose to generate a dense matrix, and find its full sparsity (no sampling)" << std::endl;

	
		std::cout << "Generating a random matrix of " << drows << " rows and " << dcols << " columns" << std::endl;
		std::cout << "with a non-zero rate of " << nnzRate * 100 << "%" << std::endl;

		//Generate a matrix	
		float *denseMatrix = new float[drows*dcols]; 
		auto gen_start = std::chrono::high_resolution_clock::now();
		generateMatrix(drows, dcols, nnzRate, denseMatrix);
		auto gen_end = std::chrono::high_resolution_clock::now();

		std::cout << "Generated Dense Matrix: " << std::endl;
		//TODO :: print as a 2D matrix ??
		printVector(denseMatrix, drows*dcols);	

		
		int *rowsVec = new int[(int)ceil(sampleRate*nnzRate*drows*dcols)];
		int *colsVec = new int[(int)ceil(sampleRate*nnzRate*drows*dcols)];
		
		//Full sparsity pattern
		auto full_start = std::chrono::high_resolution_clock::now();
		fullDense(denseMatrix, rowsVec, colsVec, drows, dcols, "Dense_mat");
		auto full_end = std::chrono::high_resolution_clock::now();
	
		auto gen_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(gen_end - gen_start).count();
		auto full_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(full_end-full_start).count();
		
		std::cout << "Matrix generation duration: " << gen_duration << " nanoseconds\n";
		std::cout << "Full sparsity pattern detection duration: " << full_duration << " nanoseconds\n";

	}

	if (fullSparse) {
		std::cout << "You chose to load a COO matrix from a MTX file, and find its full sparsity (no sampling)" << std::endl;
		
		//COO arrays 
		int* nnzRowIdx;
		int* nnzColIdx;
		float* nnzVal;
		// --end of COO arrays
		
		int n, m, nnz; 

		//Load as COO
		auto load_start = std::chrono::high_resolution_clock::now();
		loadMTX(nnzRowIdx, nnzColIdx, nnzVal, &nnz, &n, &m, filename);
		auto load_end = std::chrono::high_resolution_clock::now();

		//Full sparsity pattern
		auto full_start = std::chrono::high_resolution_clock::now();
		fullSparsity(nnz, n, m, nnzRowIdx, nnzColIdx, filename);
		auto full_end = std::chrono::high_resolution_clock::now();

		auto load_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end - load_start).count();
		auto full_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(full_end - full_start).count();

		std::cout << "File loading duration: " << load_duration << " nanoseconds\n";
		std::cout << "Full sparsity duration: " << full_duration << " nanoseconds\n";

	}

	if (sampleSparse){
		std::cout << "You chose to load a COO matrix from a MTX file, and sample it" << std::endl;
		//COO arrays 
		int* nnzRowIdx;
		int* nnzColIdx;
		float* nnzVal;
		// --end of COO arrays
		
		int n, m, nnz; 

		//Load as COO
		auto load_start = std::chrono::high_resolution_clock::now();
		loadMTX(nnzRowIdx, nnzColIdx, nnzVal, &nnz, &n, &m, filename);
		auto load_end = std::chrono::high_resolution_clock::now();

		//Full sparsity pattern
		auto sample_start = std::chrono::high_resolution_clock::now();
		sampleSparsity(sampleRate, nnzRowIdx, nnzColIdx, nnz, n, m, filename);
		auto sample_end = std::chrono::high_resolution_clock::now();

		auto load_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end - load_start).count();
		auto sample_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(sample_end - sample_start).count();

		std::cout << "File loading duration: " << load_duration << " nanoseconds\n";
		std::cout << "Sampling duration: " << sample_duration << " nanoseconds\n";

	}

	if (!(sampleLoadDense || fullLoadDense || sampleGenDense || fullGenDense || fullSparse || sampleSparse)){

		usage();
	}
	

	//freeing allocated memory
	//delete [] nnzColIdx;
	//delete [] nnzRowIdx;
	//delete [] nnzVal;

	/*
	delete [] denseMatrix;
	delete [] rowsVec;
	delete [] colsVec;
	*/
}
