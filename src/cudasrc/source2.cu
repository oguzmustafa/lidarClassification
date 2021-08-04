#include <cuda_runtime.h>
#include <source2.cuh>
#include <device_launch_parameters.h>
#include <cuda.h>
#include <helper_functions.h>
#include <c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h>
#include "c:/Users/mtf_d/Desktop/src/IPoint.h"
#include "c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h"
#include "c:/Users/mtf_d/Desktop/src/StackedPoints.h"
#include <c:/Users/mtf_d/Desktop/src/UnclassifiedPoints.h>
//#include "c:/Users/mtf_d/Desktop/src/PointVector.h"


namespace mcc {
	__global__ void basari(double ***A, int Asize) {
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		//printf("thread id : %d", threadID);
		//printf("	%f\n", A[threadID][threadID][threadID]);
		for (int i = threadID; i < Asize; i++)
		{
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("X:%f %d", A[i][j][0],threadID);
				//std::cout << A[i][j][1] << " ";
				printf("Y:%f %d", A[i][j][1], threadID);
				//std::cout << A[i][j][3] << " ";
				printf("Z:%f %d", A[i][j][2], threadID);
				//std::cout << A[i][j][4] << " ";
				printf("S:%f %d", A[i][j][3], threadID);
				printf("\n");
			}
			printf("\n");
		}
	}

	__global__ void luDeSpline(double***A, double***mtx_v, double***mtx_l, double***cells) {

	}
	__global__ void cels(double ***cells) {
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		printf("tid : %d cell.x : %f cell.y : %f\n", threadID, cells[threadID][0][1], cells[threadID][0][2]);
	}

	int source2::bos(double ***A, double ***cellsize, int Asize) {
		//allocate control_points
		std::cout << "begin d_C" << std::endl;
		double*** h_c = (double***)malloc(Asize * sizeof(double**));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			h_c[i] = (double**)malloc(twoSize * sizeof(double*));
			for (int j = 0; j < twoSize; j++) {
				cudaMalloc((void**)&h_c[i][j], 4 * sizeof(double));
				cudaMemcpy(h_c[i][j], A[i][j], 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}
		double ***h_c1 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			cudaMalloc((void***)&(h_c1[i]), twoSize * sizeof(double*));
			cudaMemcpy(h_c1[i], h_c[i], twoSize * sizeof(double*), cudaMemcpyHostToDevice);
		}
		
		double*** d_c;
		cudaMalloc((void****)&d_c, Asize * sizeof(double**));
		cudaMemcpy(d_c, h_c1, Asize * sizeof(double**), cudaMemcpyHostToDevice);
		std::cout << "end d_C" << std::endl;
		

		//allocate mtx_v
		std::cout << "begin mtx_v" << std::endl;
		double*** mtx_v1 = (double***)malloc(Asize * sizeof(double**));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3]+3;
			mtx_v1[i] = (double**)malloc(twoSize * sizeof(double*));
			for (int j = 0; j < twoSize; j++) {
				cudaMalloc((void**)&mtx_v1[i][j], 1 * sizeof(double));
			}
		}
		double ***mtx_v = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			cudaMalloc((void***)&(mtx_v[i]), twoSize * sizeof(double*));
			cudaMemcpy(mtx_v[i], h_c[i], twoSize * sizeof(double*), cudaMemcpyHostToDevice);
		}
		double*** d_mtx_v;
		cudaMalloc((void****)&d_mtx_v, Asize * sizeof(double**));
		cudaMemcpy(d_mtx_v, mtx_v, Asize * sizeof(double**), cudaMemcpyHostToDevice);
		std::cout << "end mtx_v" << std::endl;


		//allocate mtx_l
		std::cout << "begin mtx_l" << std::endl;
		double*** mtx_l1 = (double***)malloc(Asize * sizeof(double**));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3] + 3;
			mtx_l1[i] = (double**)malloc(twoSize * sizeof(double*));
			for (int j = 0; j < twoSize; j++) {
				cudaMalloc((void**)&mtx_l1[i][j], twoSize * sizeof(double));
			}
		}
		double ***mtx_l = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3] + 3;
			cudaMalloc((void***)&(mtx_v[i]), twoSize * sizeof(double*));
			cudaMemcpy(mtx_l1[i], mtx_l[i], twoSize * sizeof(double*), cudaMemcpyHostToDevice);
		}
		double*** d_mtx_l;
		cudaMalloc((void****)&d_mtx_l, Asize * sizeof(double**));
		cudaMemcpy(d_mtx_l, mtx_l, Asize * sizeof(double**), cudaMemcpyHostToDevice);
		std::cout << "end mtx_l" << std::endl;

		//allocate cells
		std::cout << "begin cells" << std::endl;
		double*** cells1 = (double***)malloc(Asize * sizeof(double**));
		for (int i = 0; i < Asize; i++) {
			int twoSize = cellsize[i][0][0];
			cells1[i] = (double**)malloc(twoSize * sizeof(double*));
			for (int j = 0; j < twoSize; j++) {
				cudaMalloc((void**)&cells1[i][j], 3 * sizeof(double));
				cudaMemcpy(cells1[i][j], cellsize[i][j], 3 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}
		double ***cells = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++) {
			int twoSize = cellsize[i][0][0];
			cudaMalloc((void***)&(mtx_v[i]), twoSize * sizeof(double*));
			cudaMemcpy(cells[i], cells1[i], twoSize * sizeof(double*), cudaMemcpyHostToDevice);
		}
		double*** d_cells;
		cudaMalloc((void****)&d_cells, Asize * sizeof(double**));
		cudaMemcpy(d_cells, cells, Asize * sizeof(double**), cudaMemcpyHostToDevice);
		std::cout << "end cells" << std::endl;
























		/*double ***d_A;
		A[0][0][0] = 999;
		printf("dA 0 0 0 : %f", A[0][0][0]);
		cudaMalloc(&d_A, Asize * sizeof(double));
		int x;
		scanf("%d", &x);
		printf("allocate started");
		for (int i = 0; i < Asize; i++) {
			int tSize = A[i][0][3];
			printf("2");
			cudaMalloc(&d_A[i], tSize * sizeof(double));
			printf("3");
			for (int j = 0; j < tSize; j++) {
				printf("4");
				cudaMalloc(&d_A[i][j], 4 * sizeof(double));
				cudaMemcpy(d_A[i][j], A[i][j], 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}
		printf("allocate completed");
		cudaMemcpy(&d_A, A, Asize, cudaMemcpyHostToDevice);
		kk << <10, 10 >> > (d_A, 3);*/

		/*double **B;
		cudaMalloc(&B, 2 * sizeof(double));
		cudaMalloc(&B[0], 2 * sizeof(double));
		cudaMalloc(&B[1], 2 * sizeof(double));
		double *C = new double[2];
		C[0] = 1;
		C[1] = 2;
		double *D = new double[2];
		D[0] = 2;
		D[1] = 3;
		cudaMemcpy(B[0], C, 2 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(B[1], D, 2 * sizeof(double), cudaMemcpyHostToDevice);
		aa << <2, 2 >> > (B);*/
		/*double *C = new double[2];
		C[0] = 1;
		C[1] = 2;
		double *D = new double[2];
		D[0] = 3;
		D[1] = 4;
		double **B;
		cudaMalloc(&B[0], 2 * sizeof(double));
		cudaMalloc(&B[1], 2 * sizeof(double));
		cudaMemcpy(B[0], C, 2 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(B[1], D, 2 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMalloc(&B, 2 * sizeof(double));
		cudaMemcpy(B, B, 2 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(B, B, 2 * sizeof(double), cudaMemcpyHostToDevice);
		aa << <2, 2 >> > ();*/







		/*int val = 0;
		for (int i = 0; i < Asize; i++) {
			val += A[i][0][3];
		}
		printf("val: %d\n", val);
		int d = val * 5;
		printf("val*5: %d\n", d);
		double *X = (double*)malloc(sizeof(double)*d);

		double a = 0, b = 0;

		for (int i = 0, int c = 0; i < Asize; i++) {
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize;) {
				if (c % 5 == 0) {
					X[c] = A[i][j][0];
					//printf("X : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 1) {
					X[c] = A[i][j][1];
					//printf("Y : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 2) {
					X[c] = A[i][j][2];
					//printf("Z : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 3) {
					X[c] = A[i][j][3];
					//printf("S : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 4) {
					a++;
					X[c] = b;
					//printf("V--->> : %f	i: %d\n", X[c], c);
					if (a == A[i][j][3]) {
						//printf("V--->> : %f	i: %d\n", X[c], c);
						a = 0;
						b++;
					}
					j++;
				}
				//printf("	a: %f", a);
				//printf("	A: %f\n", A[i][j][3]);
				if (c == d)
					break;
				c++;
			}
		}
		printf("%f\n", X[d - 1]);
		double *d_A;

		cudaMalloc(&d_A, d * sizeof(double));

		cudaMemcpy(d_A, X, d * sizeof(double), cudaMemcpyHostToDevice);

		kk << <2,3 >> > (d_A, 3);

		/*double***d_AA;

		for (int j = 0; j < Asize; j++) {
			for (int i = A[0][0][0]; i < Asize; i++) {
				cudaMalloc(&A[i][j], Asize * sizeof(float));
			}
		}
		printf("allocate started");*/
		/*for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			//printf("ilk döngü size: %d",twoSize);
			for (int j = 0; j < twoSize; j++) {
				//int arrSize = A[i][0][0];
				//printf("ikinci dongu");
				cudaMalloc(&d_AA[i][j], 4 * sizeof(double));
			}
		}*/
		
		/*double*** d_AAA = new double**[Asize];
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			d_AAA[i] = new double*[twoSize];
			for (int j = 0; j < twoSize; j++) {
				d_AAA[i][j] = new double[4];
			}
		}
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			for (int j = 0; j < twoSize; j++) {
				cudaMalloc(&d_AAA[i][j], 4 * sizeof(double));
			}
		}
		printf("allocate completed");

		printf("%d", Asize);

		for (int a = 0; a < Asize; a++) {
			int twoArrSize = A[a][0][3];
			//printf("ilk döngü size: %d", twoArrSize);
			for (int b = 0; b < twoArrSize; b++) {
				//printf("ikinci dongu");
				cudaMemcpy(d_AAA[a][b], A[a][b], 4 * sizeof(double), cudaMemcpyHostToDevice);
				//printf("atama sonrasý d_A: %f", d_A[0][0][3]);
			}
		}
		printf("copy completed");

		kk << <2, 2 >> > (d_AAA, 3);*/

		return 1;
	}
}

