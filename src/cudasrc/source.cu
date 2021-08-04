#include <cuda_runtime.h>
#include <source.cuh>
#include <device_launch_parameters.h>
#include <cuda.h>
#include <helper_functions.h>
#include <c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h>
#include "c:/Users/mtf_d/Desktop/src/IPoint.h"
#include "c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h"
#include "c:/Users/mtf_d/Desktop/src/StackedPoints.h"
#include <c:/Users/mtf_d/Desktop/src/UnclassifiedPoints.h>
//#include "c:/Users/mtf_d/Desktop/src/PointVector.h"


namespace mcc{

	#include <cstdio>
		inline void GPUassert(cudaError_t code, char * file, int line, bool Abort = true)
		{
			if (code != 0) {
				fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
				if (Abort) exit(code);
			}
		}

	#define GPUerrchk(ans) { GPUassert((ans), __FILE__, __LINE__); }


	
	__global__ void testkernel(double *A, int arraySize) {
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		printf("threadid : %d", threadID);
		if (threadID < 10) {
			printf("threadid : %d", threadID);
		}

		/*if (threadID < 10) {
			for (int i = 0, int c = 0; i < 2; i++) {
				for (int j = 0; j < 25;) {
					if (c % 5 == 0) {
						printf("X : %f	i: %d	", A[c], c);
					}
					if (c % 5 == 1) {
						printf("Y : %f	i: %d	", A[c], c);
					}
					if (c % 5 == 2) {
						printf("Z : %f	i: %d	", A[c], c);
					}
					if (c % 5 == 3) {
						printf("S : %f	i: %d	", A[c], c);
					}
					if (c % 5 == 4) {
						printf("V--->> : %f	i: %d\n", A[c], c);
						j++;
					}
					c++;
				}
			}
		}*/



		/*if (threadID < 10) {
			/*for (int d = 0; d < 5; d++) {
				A[0] += A[0];
				printf("in device value of -->>> %f\n", A[0]);
			}*/
			/*printf("%f", A[0][0][0]);
			int size = (int)A[0][0][3];
			for (int j = 0; j < size; j++)
			{
				printf("%f ", A[threadID][j][0]);
				printf("%f ", A[threadID][j][1]);
				printf("%f ", A[threadID][j][2]);
				printf("%f \n", A[threadID][j][3]);
			}

		}*/
			
	}
	__global__ void vectorAdditionKernel(double* A, double* B, double* C, int arraySize) {
		// Get thread ID.
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;

		// Check if thread is within array bounds.
		if (threadID < arraySize) {
			// Add a and b.
			C[threadID] = A[threadID] + B[threadID];
			printf("global ici");
		}
	}
	__global__ void vkadd(double* A, double* B, int arraySize) {
		// Get thread ID.
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;

		// Check if thread is within array bounds.
		if (threadID < arraySize) {
			// Add a and b.
			printf("global ici vkADD");
		}
	}

	void source::kernel(double* A, double* B, double* C, int arraySize) {

		/*int numElements = 5000000;
		size_t size = numElements * sizeof(float);
		float *h_A = (float *)malloc(size);

		CUdeviceptr d_X;
		cudaMalloc((void**)&d_X, size * sizeof(double));*/

		// Initialize device pointers.

		double *xa = (double*)malloc(sizeof(double) * 3);
		xa[0] = 0;
		xa[1] = 1;
		xa[2] = 2;
		double *xb = (double*)malloc(sizeof(double) * 3);
		xb[0] = 0;
		xb[1] = 1;
		xb[2] = 2;
		double*xx;
		double*xy;
		cudaMalloc(&xx, 3 * sizeof(double));
		cudaMalloc(&xy, 3 * sizeof(double));

		cudaMemcpy(xx, xa, 3 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(xy, xb, 3 * sizeof(double), cudaMemcpyHostToDevice);

		double* d_A, *d_B, *d_C;

		// Allocate device memory.
		cudaMalloc(&d_A, arraySize * sizeof(double));
		cudaMalloc(&d_B, arraySize * sizeof(double));
		cudaMalloc(&d_C, arraySize * sizeof(double));

		// Transfer arrays a and b to device.
		cudaMemcpy(d_A, A, 5 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(d_B, B, arraySize * sizeof(double), cudaMemcpyHostToDevice);
		//printf("kernel deger icinde d_A degeri : %f", d_A);

		// Calculate blocksize and gridsize.
		dim3 blockSize(512, 1, 1);
		dim3 gridSize(512 / arraySize + 1, 1);

		// Launch CUDA kernel.

		vectorAdditionKernel << <1,1 >> > (xx, xy, d_C, 3);
		//vectorAdditionKernel << <gridSize, blockSize >> > (d_A, d_B, d_C, arraySize);

		// Copy result array c back to host memory.
		cudaMemcpy(C, d_C, arraySize * sizeof(double), cudaMemcpyDeviceToHost);

		/*int x = 0;
		scanf("%d",&x);
		printf("cuda kernel matris toplamý\n");*/

	}

	void source::clas(IUnclassifiedPoints & points) {
		std::cout << "cu icine points atma, sayýsý:--> " << points.count() << std::endl;
	}

	void source::test(double ***A, int Asize) {
		int arraySize = 3;
		double *xa = (double*)malloc(sizeof(double) * 3);
		xa[0] = 0;
		xa[1] = 1;
		xa[2] = 2;
		double *xb = (double*)malloc(sizeof(double) * 3);
		xb[0] = 0;
		xb[1] = 1;
		xb[2] = 2;
		double*xx;
		double*xy;
		cudaMalloc(&xx, 3 * sizeof(double));
		cudaMalloc(&xy, 3 * sizeof(double));

		cudaMemcpy(xx, xa, 3 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(xy, xb, 3 * sizeof(double), cudaMemcpyHostToDevice);

		double* d_A, *d_B, *d_C;

		// Allocate device memory.
		cudaMalloc(&d_A, arraySize * sizeof(double));
		cudaMalloc(&d_B, arraySize * sizeof(double));
		cudaMalloc(&d_C, arraySize * sizeof(double));

		// Transfer arrays a and b to device.
		//cudaMemcpy(d_A, A, 5 * sizeof(double), cudaMemcpyHostToDevice);
		//cudaMemcpy(d_B, B, arraySize * sizeof(double), cudaMemcpyHostToDevice);
		//printf("kernel deger icinde d_A degeri : %f", d_A);

		// Calculate blocksize and gridsize.
		dim3 blockSize(512, 1, 1);
		dim3 gridSize(512 / arraySize + 1, 1);

		// Launch CUDA kernel.

		testkernel << <1, 1 >> > (xx, 3);
		/*int val = 0;
		for (int i = 0; i < Asize; i++) {
			val += A[i][0][3];
		}
		printf("val: %d\n", val);
		int d = val * 5;
		printf("val*5: %d\n", d);
		double *X = (double*)malloc(sizeof(double)*d);

		double a = 0, b = 0;

		for (int i = 0,int c = 0; i < Asize; i++) {
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
		printf("%f\n", X[d-1]);

		double *xa = (double*)malloc(sizeof(double) * 3);
		xa[0] = 0;
		xa[1] = 1;
		xa[2] = 2;
		double *xb = (double*)malloc(sizeof(double) * 3);
		xb[0] = 0;
		xb[1] = 1;
		xb[2] = 2;
		double*xx;
		double*xy;
		cudaMalloc(&xx, 3 * sizeof(double));
		cudaMalloc(&xy, 3 * sizeof(double));
		double *d_C;
		cudaMalloc(&d_C, 3 * sizeof(double));

		cudaMemcpy(xx, xa, 3 * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(xy, xb, 3 * sizeof(double), cudaMemcpyHostToDevice);

		vkadd << <10, 10 >> > (xx, xy, d_C, 3);

		/*double *d_A;

		cudaMalloc(&d_A, 3 * sizeof(double));

		cudaMemcpy(d_A, xa, 3 * sizeof(double), cudaMemcpyHostToDevice);
		
		testkernel <<<10, 10 >>> (d_A, 3);
		
		for (int i = 0, int c = d-120; i < 2; i++) {
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize;) {
				if (c % 5 == 0) {
					printf("X : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 1) {
					printf("Y : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 2) {
					printf("Z : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 3) {
					printf("S : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 4) {
					printf("V--->> : %f	i: %d\n", X[c], c);
					j++;
				}
				c++;
			}
		}

		for (int i = d-25; i < d; i++)
		{
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("%f ", A[i][j][0]);
				//std::cout << A[i][j][1] << " ";
				printf("%f ", A[i][j][1]);
				//std::cout << A[i][j][3] << " ";
				printf("%f ", A[i][j][2]);
				//std::cout << A[i][j][4] << " ";
				printf("%f ", A[i][j][3]);
				std::cout << std::endl;
			}
		}*/
























		/*for (int i = 0; i < 2; i++)
		{
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("%f ", A[i][j][0]);
				//std::cout << A[i][j][1] << " ";
				printf("%f ", A[i][j][1]);
				//std::cout << A[i][j][3] << " ";
				printf("%f ", A[i][j][2]);
				//std::cout << A[i][j][4] << " ";
				printf("%f ", A[i][j][3]);
				std::cout << std::endl;
			}
		}*/
		/*for (int i = 0, int c = 0; i < 2; i++) {
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize;) {
				if (c % 5 == 0) {
					printf("X : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 1) {
					printf("Y : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 2) {
					printf("Z : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 3) {
					printf("S : %f	i: %d	", X[c], c);
				}
				if (c % 5 == 4) {
					printf("V--->> : %f	i: %d\n", X[c], c);
					j++;
				}
				c++;
			}
		}*/





		//printf("X[0]: %f --- A: %f", X[0]);

		/*for (int i = 0; i < 10; i++) {
			if (i % 5 == 0) {
				printf("bu x: %d", X[i]);
				printf("	bu da Ax: %d\n", A[i][0][0]);
			}
			if (i % 5 == 1) {
				printf("bu y:%d", X[i]);
				printf("	bu da Ay: %d\n", A[i][0][1]);
			}
			if (i % 5 == 2) {
				printf("bu z:%d", X[i]);
				printf("	bu da Az: %d\n", A[i][0][2]);
			}
			if (i % 5 == 3) {
				printf("bu size:%d", X[i]);
				printf("	bu da Asize: %d\n", A[i][0][3]);
			}
			if (i % 5 == 4) {
				printf("bu val: %d\n", X[i]);
				//printf("bu da Ax: ", A[i][0][4]);
			}
		}*/

		/*for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("%f ", A[i][j][0]);
				//std::cout << A[i][j][1] << " ";
				printf("%f ", A[i][j][1]);
				//std::cout << A[i][j][3] << " ";
				printf("%f ", A[i][j][2]);
				//std::cout << A[i][j][4] << " ";
				printf("%f ", A[i][j][3]);
				std::cout << std::endl;
			}
		}*/




		/*double ***d_A;
		cudaMalloc(&d_A, Asize * sizeof(double));

		for (int i = 0; i < Asize; i++) {
			int tSize = A[i][0][3];
			cudaMalloc(&d_A[i], tSize * sizeof(double));
			for (int j = 0; j < tSize; j++) {
				cudaMalloc(&d_A[i][j], 4 * sizeof(double));
				cudaMemcpy(&d_A[i][j], A[i][j], 4 * sizeof(double), cudaMemcpyHostToDevice);
			}
		}
		cudaMemcpy(&d_A, A, Asize, cudaMemcpyHostToDevice);
		testkernel << <10, 10 >> > (d_A, 3);*/


		/*printf("malloc begin");
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			for (int j = 0; j < twoSize; j++) {
				GPUerrchk(cudaMalloc(&A[i][j], 4 * sizeof(double)));
			}
		}
		printf("malloc comp");

		for (int i = 0; i < Asize; i++) {
			int s = A[i][0][3];
			GPUerrchk(cudaMalloc(&A[i], s * sizeof(double)));
			GPUerrchk(cudaMemcpy(A[i], A[i], s * sizeof(double), cudaMemcpyHostToDevice));
		}
		printf("memcpy comp ");

		

		GPUerrchk(cudaMalloc(&d_A, Asize * sizeof(double)));
		GPUerrchk(cudaMemcpy(d_A, A, Asize * sizeof(double), cudaMemcpyHostToDevice));
		testkernel << <10, 10 >> > (d_A, 3);*/
		//double ***d_A;
		/*for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			//printf("ilk döngü size: %d",twoSize);
			for (int j = 0; j < twoSize; j++) {
				//int arrSize = A[i][0][0];
				//printf("ikinci dongu");
				cudaMalloc(&d_A[i][j], 4 * sizeof(double));
				cudaMemcpy(d_A[i][j], A[i][j], 4 * sizeof(double),cudaMemcpyHostToDevice);
			}
		}
		for (int j = 0; j < Asize; j++) {
			int twoSize = A[j][0][3];
			for (int i = 0; i < twoSize; i++) {
				cudaMalloc(&d_A[i], twoSize * sizeof(double));
			}
		}

		cudaMalloc(&d_A, Asize * sizeof(double));*/

		/*
		for (int a = 0; a < Asize; a++) {
			int twoArrSize = A[a][0][3];
			//printf("ilk döngü size: %d", twoArrSize);
			for (int b = 0; b < twoArrSize; b++) {
				//printf("ikinci dongu");
				cudaMemcpy(d_A[a][b], A[a][b], 4 * sizeof(double), cudaMemcpyHostToDevice);
				//printf("atama sonrasý d_A: %f", d_A[0][0][3]);
			}
		}*/
		


		/*int val = 0;
		for (int i = 0; i < Asize; i++) {
			val += A[i][0][3];
			//val += val;
		}
		printf("val %d\n", val);
		int a[10];
		for (int i = 1; i < 11; i++) {
			a[i] = i;
		}
		int deg = 0;
		for (int i = 1; i < 11; i++) {
			deg += a[i];
			printf("%d\n", a[i]);
		}
		printf("deg: %d\n", deg);*/



		//printf("cu icinde test %f", A[0][0][0]);
		//std::cout << "cu icinde test" << A[0] << std::endl;
		/*
		double* d_A;
		double* P;
		P = A;
		cudaMalloc((void**)&d_A,sizeof(double));
		cudaMemcpy(d_A, P, 1 * sizeof(double), cudaMemcpyHostToDevice);

		test << <2,2 >> > (d_A);
		*/
		
		/*for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				//std::cout << A[i][j][0] << " ";
			printf("%f ", A[i][j][0]);
			//std::cout << A[i][j][1] << " ";
			printf("%f ", A[i][j][1]);
			//std::cout << A[i][j][3] << " ";
			printf("%f ", A[i][j][2]);
			//std::cout << A[i][j][4] << " ";
			printf("%f ", A[i][j][3]);
			}
			std::cout << std::endl;
		}*/
		/*for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("%f ", A[i][j][0]);
				//std::cout << A[i][j][1] << " ";
				printf("%f ", A[i][j][1]);
				//std::cout << A[i][j][3] << " ";
				printf("%f ", A[i][j][2]);
				//std::cout << A[i][j][4] << " ";
				printf("%f ", A[i][j][3]);
				std::cout << std::endl;
			}
		}*/

		/*printf("A dizisi boyutu : %d\n", Asize);


		// Initialize device pointers.
		//double*** d_A[Asize];
		double*** d_A = new double**[Asize];
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			d_A[i] = new double*[twoSize];
			for (int j = 0; j < twoSize; j++) {
				d_A[i][j] = new double[4];
			}
		}*/

		//double ***d_A = A;
		/*printf("d_A deger : %f\n", d_A[0][0][3]);
		printf("A value : %f", A[0][0][3]);
		d_A[0][0][3] = 5;
		printf("d_A deger : %f\n", d_A[0][0][3]);
		printf("A value : %f\n", A[0][0][3]);
		*/
		// Allocate device memory.
		//cudaMalloc((void****)&d_A, arraySize * sizeof(double));

		/*for (int j = 0; j < Asize; j++) {
			for (int i = A[0][0][0]; i < Asize; i++) {
				cudaMalloc(&A[i][j], Asize * sizeof(float));
			}
		}*/
		/*printf("allocate started");
		for (int i = 0; i < Asize; i++) {
			int twoSize = A[i][0][3];
			//printf("ilk döngü size: %d",twoSize);
			for (int j = 0; j < twoSize; j++) {
				//int arrSize = A[i][0][0];
				//printf("ikinci dongu");
				cudaMalloc(&d_A[i][j], 4 * sizeof(double));
			}
		}*/
		//printf("allocate completed");
		/*
		int x = 0;
		scanf("%d", &x);
		*/
		/*for (int i = 0; i < n; i++) {
			cudaMemcpy(temph[i], a[i], n * sizeof(float), cudaMemcpyHostToDevice);
		}*/

		//printf("%d", Asize);

		/*for (int a = 0; a < Asize; a++) {
			int twoArrSize = A[a][0][3];
			//printf("ilk döngü size: %d", twoArrSize);
			for (int b = 0; b < twoArrSize; b++) {
				//printf("ikinci dongu");
				cudaMemcpy(d_A[a][b], A[a][b], 4 * sizeof(double), cudaMemcpyHostToDevice);
				//printf("atama sonrasý d_A: %f", d_A[0][0][3]);
			}
		}
		printf("copy completed");*/


		// Transfer arrays a and b to device.
		//cudaMemcpy(d_A, A, arraySize * sizeof(double), cudaMemcpyHostToDevice);

		// Calculate blocksize and gridsize.
		/*dim3 blockSize(512, 1, 1);
		dim3 gridSize(512 / 3 + 1, 1);

		// Launch CUDA kernel.
		testkernel << <blockSize,gridSize >> > (d_A,3);*/
		printf("kernel alti");
}
}