/* By Khaled Abdelaal 	khaled.abdelaal@ou.edu
 * 
 * This program performs (dense) matrix multiplication 
 * C = A * B
 * Where: A, B, C are n * n dense matrices
 * 	  
 *
 * Multiplication is done in three different variants:
 * [1] Naive implementation: three nested loops, all data in row-major format
 * [2] Column-major for B: A in row major while B in col major 
 * [3] Tiling: a BLOCK_SIZE of n/10 is used for tiling
 *
 *
 * variants 2, and 3 objetive is to exploit cache locality (spatial and temporal)
 * to improve performance
 *
 * A and B are populated with random values,
 * then each variant is executed NUM_SAMPLES times
 * average execution time for each variant is captured and printed
 *
 *
 * Program doesn't take any arguments
 *
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define n 100
#define m 100
#define k 100
#define BLOCK_SIZE n/10
#define NUM_SAMPLES 5
#define lower 0
#define upper 10


int main() {
		
	srand(time(0));
	
	
	double time_naive, time_colmaj, time_tiling;

	static int B[n][m];
	static int C[m][k];

	//C transpose
	static int Ct[k][m];

	static int A[n][k] = {0};
	//int A[n][k];

	for (int x = 0; x < n; x++){

		for (int y = 0; y < m; y++){
					
			B[x][y] = (rand() % (upper - lower + 1)) + lower;
		}
	}

	for (int x = 0; x < m; x++){

		for (int y = 0; y <k ; y++){
					
			C[x][y] = (rand() % (upper - lower + 1)) + lower;

			//calculate D which is C Transpose
			Ct[y][x] = C[x][y];
		}
	}
	
	/*
	printf("Matrix B:\n");
	for (int x = 0; x < n; x++){

		for (int y = 0; y<m ;y++){

			printf("%d ",B[x][y]);
		}
		printf("\n");
	}


	printf("Matrix C:\n");
	for (int x = 0; x < m; x++){

		for (int y = 0; y<k ;y++){

			printf("%d ",C[x][y]);
		}
		printf("\n");
	}
	
	*/

	
	for (int ss=0; ss < NUM_SAMPLES; ss++){


		clock_t start_instance = clock();

		for (int i =0; i< n; i++)
		{
			for (int j = 0 ; j < k ; j++){
			
				for (int l = 0; l < m; l++){

					A[i][j] += B[i][l] * Ct[j][l];
				}
			}
		}	

		clock_t end_instance = clock();
		time_colmaj += ((double) (end_instance - start_instance)) / CLOCKS_PER_SEC;
	}

	double avg_time_colmaj = time_colmaj / NUM_SAMPLES;
	
		
	for (int ss = 0 ; ss < NUM_SAMPLES; ss++){
		clock_t start_naive = clock();

		for (int i =0; i< n; i++)
		{
			for (int j = 0 ; j < k ; j++){
				
				for (int l = 0; l < m; l++){

					A[i][j] += B[i][l] * C[l][j];
				}
			}
		}	

		clock_t end_naive = clock();
		time_naive += ((double) (end_naive - start_naive)) / CLOCKS_PER_SEC;
	}
	double avg_time_naive = time_naive / NUM_SAMPLES;
	
	/*
	printf("Result A (original):\n");
	for (int x = 0; x < n; x++){

		for (int y = 0; y<k ;y++){

			printf("%d ",A[x][y]);
		}
		printf("\n");
	}
	*/
	

	for (int i = 0 ; i < n ; i++){

		for (int j = 0; j < k; j++){
			A[i][j] = 0;
		}
	}

	int r;

	for (int ss = 0 ; ss < NUM_SAMPLES; ss++){
		clock_t start_tiling = clock();

		for (int i =0; i< n; i += BLOCK_SIZE)
		{
			for (int j = 0 ; j < k ; j += BLOCK_SIZE){
				
				for (int d = 0; d < m; d++){

					for (int ii = i ; ii < i+BLOCK_SIZE; ii++){
						r = 0;
						for (int jj = j; jj < j+BLOCK_SIZE; jj++){

							r = r + B[d][jj] * C[jj][ii];
						}
						A[d][ii] = A[d][ii] + r;
					}
				}
			}
		}	

		clock_t end_tiling = clock();
		time_tiling += ((double) (end_tiling - start_tiling)) / CLOCKS_PER_SEC;
	}
	double avg_time_tiling = time_tiling / NUM_SAMPLES;

	/*
	printf("Result A (tiled):\n");
	for (int x = 0; x < n; x++){

		for (int y = 0; y<k ;y++){

			printf("%d ",A[x][y]);
		}
		printf("\n");
	}
	*/
	printf("\nExecution time when C is stored in column major format: %f seconds", avg_time_colmaj);
	printf("\nExecution time when both B, and C are stored row-major: %f seconds\n", avg_time_naive);
	printf("\nExecution time using tiling (blocking): %f seconds\n", avg_time_tiling);
	
}
