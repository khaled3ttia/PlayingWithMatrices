#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#define n 3
#define m 3
#define k 2
#define lower 0
#define upper 100
int main() {

	srand(time(0));

	int B[n][m];
	int C[m][k];
	int A[n][k] = {0};
	int nnz = 0;
	int idx = 0; 
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


	for (int x = 0 ; x < n; x++) {

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

					
				idx++;
			}
	
			
			
		}
		
	}

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




}
