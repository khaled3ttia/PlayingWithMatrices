#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define n 2
#define m 2 
#define k 2
#define lower 0
#define upper 10
int main() {

	srand(time(0));

	int B[n][m];
	int C[m][k];
	int A[n][k] = {0};

	for (int x = 0; x < n; x++){

		for (int y = 0; y < m; y++){
					
			B[x][y] = (rand() % (upper - lower + 1)) + lower;
		}
	}

	for (int x = 0; x < m; x++){

		for (int y = 0; y < k; y++){
					
			C[x][y] = (rand() % (upper - lower + 1)) + lower;
		}
	}

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


	for (int i =0; i< n; i++)
	{
		for (int j = 0 ; j < k ; j++){
			
			for (int l = 0; l < m; l++){

				A[i][j] += B[i][l] * C[l][j];
			}
		}
	}	



	printf("Result A:\n");
	for (int x = 0; x < n; x++){

		for (int y = 0; y<k ;y++){

			printf("%d ",A[x][y]);
		}
		printf("\n");
	}


}
