// By: Khaled Abdelaal
// khaled.abdelaal@ou.edu
// December 2020

// Utility to load MTX files and find the full sparsity pattern as a plot

#include <math.h>
#include <unistd.h>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>

//This function takes a MTX file as an input and loads it in the Coordinate format
// arguments:
//	*&nnzRowIdx : (in-out) an integer pointer reference to the output COO row indices array
//	*&nnzColIdx : (in-out) an integer pointer reference to the output COO col indices array
//	*&nnzVal    : (in-out) a float pointer reference to the output COO float values array
//	*nnz 	    : (in-out) an integer pointer to the number of non-zeros in the matrix
//	filename    : (in)     a string representing filename (file path)
void loadMTX(int *&nnzRowIdx, int *&nnzColIdx, float *&nnzVal, int* nnz, std::string filename){

	std::cout << "Finding sparsity pattern for " << filename << "....." <<  std::endl;
	
	int nrows, ncols;
	
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
				if (!(iss >> nrows >> ncols >> *nnz)) {break;}
				
				std::cout << "number of rows is: " << nrows << std::endl;
				std::cout << "number of cols is: " << ncols << std::endl;
				std::cout << "nnz is: " << *nnz << std::endl;	
				
				// dynamically allocate memory for COO arrays based on nnz				
				nnzColIdx = new int[*nnz];
				nnzRowIdx = new int[*nnz];
				nnzVal = new float[*nnz];

				// initialize all indices to -1
				
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
					
					nnzColIdx[currentIndex] = colIndex;
					nnzRowIdx[currentIndex] = rowIndex; 
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
void fullSparsity(int nnz, int *&nnzRowIdx, int *&nnzColIdx, std::string filename){
		
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
		outputfile << "plt.ylabel('row index')\n";
		outputfile << "plt.plot(cols, rows)\n";
		outputfile << "plt.savefig('sp_" + filename +   ".png')\n";
	}
	outputfile.close();

	// execute the script to generate the figure
	// NOTE: needs python3 installation
	system("python plot.py");
	
	std::cout << "Sparsity pattern figure saved at sp_" << filename << ".png !" << std::endl;
	

}

// A function to print out how to use the tool
void usage(){
	std::cout << "Options:" << std::endl;
	std::cout << "-i <inputfile> : filename of a MTX file" << std::endl;
	std::cout << "-h : prints usage information" << std::endl;
	exit(0);

}

// The main function

int main(int argc, char** argv){

	int opt;
	int nnz;

	//COO arrays 
	int* nnzRowIdx;
	int* nnzColIdx;
	float* nnzVal;
	// --end of COO arrays

	std::string filename; 

	// parse command line arguments
	while ((opt = getopt(argc, argv, "i:h")) != -1){
		switch (opt) {
			case 'i': filename = optarg; break;
			case 'h': usage(); break;

		}

	}	
	if (filename.empty()){
		usage();
	}
	
	// load the MTX file and time the loading process
	auto load_start_time = std::chrono::high_resolution_clock::now();
	loadMTX(nnzRowIdx, nnzColIdx, nnzVal, &nnz, filename);
	auto load_end_time = std::chrono::high_resolution_clock::now();
	
	// generate the full sparsity plot and time the process
	auto full_start_time = std::chrono::high_resolution_clock::now();
	fullSparsity(nnz, nnzRowIdx, nnzColIdx, filename);
	auto full_end_time = std::chrono::high_resolution_clock::now();
	
	auto load_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end_time-load_start_time).count();

	auto full_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(full_end_time-full_start_time).count();

	std::cout << "File loading duration: " << load_duration << " nanoseconds\n";
	std::cout << "Full pattern detection duration: " << full_duration << " nanoseconds\n";
}
