/* By Khaled Abdelaal  khaled.abdelaal@ou.edu
 *
 * This program converts dense matrices into 3 popular
 * spare matrix representations: CSR, COO, and LIL
 *
 * Matrix B of size n * m is randomly generated 
 * range of values in B is between lower and upper
 * at least 80% of the values of B will be zeros 
 * to motivate the need for sparse representations
 *
 * The dense representation is then converted into :
 * [1] Compressed Sparse Row: Three arrays are maintained
 * 		      (a) vals: contains the non-zero values in row order
 * 		      (b) b_crd: contains column indices for non-zero values
 * 		      (c) pos: contains the start index for non-zero values for this row in b_crd 
 *
 * [2] Coordinate representation: Three arrays are maintained
 * 		      (a) vals: contains the non-zero values in row order (doesn't have to be
 * 		      in any order, but to re-use vals from CSR, the assumption here is order)
 * 		      (b) coo_row: row indices for non-zero vals (from vals)
 * 		      (c) coo_col: column indices for non-zero vals (from vals)
 *
 * [3] List of Lists: Each row is represented as a list of tuples for non-zero values
 * 		      each tuple consists of column index and the non-zero value
 *
 *
 *
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#define n 3
#define m 3
#define lower 0
#define upper 100

//Tuples for LIL format
//each row is a list of tuple
//each tuple contains the column index
//and the value 
typedef struct tuple{
	int col;
	int val;

};

int main() {

	srand(time(0));

	int B[n][m];
	int nnz = 0;
	int idx = 0; 
	
	//Keeping track of nnz per row
	//useful for formats such as LIL
	int nnz_row[n] = {0};

	//Populating Matrix B with random values
	for (int x = 0; x < n; x++){

		for (int y = 0; y < m; y++){
					
			int randNum = (rand() % (upper - lower + 1)) + lower;
			if (randNum < 0.8 * upper){
				B[x][y] = 0;
			}
			else {
				B[x][y] = randNum;
			}
			if (B[x][y] != 0){
				
				nnz_row[x]++;			
				nnz++;

			}
		}
	}

	int pos[n+1] = {0};
	pos[n] = nnz;
	int b_crd[nnz];
	int vals[nnz];
	
	int coo_row[nnz];
	int coo_col[nnz];

	//LIL represtnation, we have n rows, each row is a list of size = nnz per row
	//each entry in a row list will be a tuple of (col, val)
	struct tuple* lil_rows[n];
	//dynamically allocating the tuples for each row depending on the nnz per row
	for (int i = 0; i < n; i++){

		lil_rows[i] = (struct tuple*)malloc(sizeof(struct tuple)* nnz_row[i]);

	}

	for (int x = 0 ; x < n; x++) {
		
		//pointer to the current row in LIL representation
		struct tuple* currentRow = lil_rows[x];

		for (int y=0; y<m; y++){
			
			if (B[x][y] != 0) { 
				
				//CSR representation
				vals[idx] = B[x][y];
				b_crd[idx] = y;
				
				//COO representation
				coo_row[idx] = x;
				//for COO col array, we can use the same b_crd from CSR
				//(since here we are maintaing order for NNZ values in
				//COO, while it's not required)
			

				//LIL representation
				//Each tuple of the row will have a (col, val) pair	
				currentRow->col = y;
				currentRow->val = B[x][y];
				currentRow++;

					
				idx++;
			}
	
			
			
		}
		
	}

	//CSR related 
	//Constructing the pos array for CSR
	int pos_idx = 1;
	bool row_start_found = false;
	int lastpos = 0;
	int crd_idx = 0;
	for (int i = 0; i < n; i++){
		row_start_found = false;
		for (int j = 0 ; j<m; j++){
			if ((row_start_found==false) && B[i][j] != 0){
				for (int o=lastpos; o < nnz; o++){
					if (b_crd[o] == j){
						pos[i] = o;
						lastpos = o+1;
						row_start_found = true;
						break;	
					}
				}		
			}
		}
	}


	//Printing out results to visualize all different representations
	//along with the original dense matrix
	printf("Matrix B:\n");
	
	for (int x = 0; x < n; x++){

		for (int y = 0; y<m ;y++){

			printf("%d ",B[x][y]);
		}
		printf("\n");
	}
	
	printf("CSR representation: \n");
	printf("B_pos:[");
	for (int i=0; i < n+1; i++){
		if (i == n){
			printf("%d]\n", pos[i]);
		}else {
			printf("%d, ",pos[i]);
		}
	}

	printf("B_crd:[");
	for (int i=0; i<nnz; i++){
		if (i == nnz-1){
			printf("%d]\n", b_crd[i]);
		}else {
			printf("%d, ", b_crd[i]);
		}
	}


	printf("vals:[");
	for (int i=0; i<nnz; i++){
		if (i == nnz-1){
			printf("%d]\n", vals[i]);
		}else {
			printf("%d, ", vals[i]);
		}
	}

	printf("==========================\n");
	printf("COO representation: \n");
	printf("row_idx:[");
	for (int i = 0 ; i < nnz; i++) {
		if (i == nnz-1){
			printf("%d]\n", coo_row[i]);

		}else { 
			printf("%d, " , coo_row[i]);
		}
	}

	printf("col_idx:[");
	for (int i = 0; i < nnz; i++){
		if (i == nnz-1){
			printf("%d]\n", b_crd[i]);
		}else {
			printf("%d, ", b_crd[i]);
		}
	}

	printf("vals:[");
	for (int i=0; i<nnz; i++){
		if (i == nnz-1){
			printf("%d]\n", vals[i]);
		}else {
			printf("%d, ", vals[i]);
		}
	}


	printf("==========================\n");
	printf("LIL representation: \n");
	int ff=0;
	for (int i = 0 ; i < n ; i++){
		struct tuple* r = lil_rows[i];
		printf("row %d:", i);
		for (int j = 0; j < nnz_row[ff]; j++){
			printf("(%d,%d)", r->col, r->val);
			r++;
		}
		printf("\n");
		ff++;
		lil_rows[i]++;
	}



}
