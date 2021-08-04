#include <cuda_runtime.h>
#include <source3.cuh>
#include <device_launch_parameters.h>
#include <cuda.h>
#include <helper_functions.h>
#include <c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h>
#include "c:/Users/mtf_d/Desktop/src/IPoint.h"
#include "c:/Users/mtf_d/Desktop/src/IUnclassifiedPoints.h"
#include "c:/Users/mtf_d/Desktop/src/StackedPoints.h"
#include <c:/Users/mtf_d/Desktop/src/UnclassifiedPoints.h>
#include <windows.h>
#include <cstdio>
#include <ctime>
#include <curand.h>

//#include "c:/Users/mtf_d/Desktop/src/PointVector.h"

namespace mcc
{
	inline void GPUassert(cudaError_t code, char *file, int line, bool Abort = true)
	{
		if (code != 0)
		{
			fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
			if (Abort)
				exit(code);
		}
	}

#define GPUerrchk(ans)                        \
	{                                         \
		GPUassert((ans), __FILE__, __LINE__); \
	}

	__global__ void basarii(double ***A, int Asize)
	{
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		//printf("thread id : %d", threadID);
		//printf("	%f\n", A[threadID][threadID][threadID]);
		for (int i = threadID; i < Asize; i++)
		{
			int tSize = A[i][0][3];
			for (int j = 0; j < tSize; j++)
			{
				//std::cout << A[i][j][0] << " ";
				printf("X:%f %d", A[i][j][0], threadID);
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

	/*__global__ void cel(double ***A)
	{
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		if (threadID == 0)
		{
			printf("in global");
		}
		int t = A[threadID][0][3];
		/*for (int j = 0; j < t; j++) {
			//std::cout << A[i][j][0] << " ";
			printf("X:%f %d", A[threadID][j][0], threadID);
			//std::cout << A[i][j][1] << " ";
			printf("Y:%f %d", A[threadID][j][1], threadID);
			//std::cout << A[i][j][3] << " ";
			printf("Z:%f %d", A[threadID][j][2], threadID);
			//std::cout << A[i][j][4] << " ";
			printf("S:%f %d", A[threadID][j][3], threadID);
			printf("\n");
			}*/
	//printf("global ici double pointer dec\n");
	/*double **res = new double *[t];
		for (int i = 0; i < t; i++)
		{
			res[i] = new double[4];
		}
		//printf("global ici double pointer id\n");
		for (int i = 0; i < t; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				res[i][j] = A[threadID][i][j];
			}
		}
		for (int i = 0; i < t; i++)
		{
			printf("global ici res x : %f", res[i][0]);
			printf("global ici res y : %f", res[i][1]);
			printf("global ici res z : %f", res[i][2]);
			printf("global ici res s : %f", res[i][3]);
			printf("threadID : %d", threadID);
			printf("\n");
		}

		printf("global	");
		printf("%d\n", threadID);
		printf("cells x: %f cells y: %f", cells[threadID][0][1], cells[threadID][0][2]);*/
	/*int twos = cells[threadID][0][0];
		printf("global ici double pointer dec\n");
		double **res = new double*[twos];
		for (int i = 0; i < 8; i++) {
		res[i] = new double[3];
		}
		printf("global ici double pointer id\n");
		for (int i = 0; i < twos; i++) {
		for (int j = 0; j < 3; j++) {
		res[i][j] = cells[threadID][i][j];
		}
		}
		for (int i = 0; i < twos; i++) {
		for (int j = 0; j < 3; j++) {
		printf("global ici res : %f\n", res[i][j]);
		}
		}
	}*/
	/*if (threadID == 630) {
		for (int i = 630; i < 631; i++) {
			int twoS = cells[i][0][0];
			for (int j = 0; j < twoS; j++) {
				printf("\nres %d : %f		", i, cells[i][j][1]);
				printf("res %d : %f\n", i, cells[i][j][2]);
			}
			printf("\n");
		}
		printf("global ici 631.deger\n");
		for (int j = 0; j < points[631][0][3]; j++) {
			//std::cout << A[i][j][0] << " ";
			printf("X:%f		", points[631][j][0]);
			//std::cout << A[i][j][1] << " ";
			printf("Y:%f		", points[631][j][1]);
			//std::cout << A[i][j][3] << " ";
			printf("Z:%f		", points[631][j][2]);
			//std::cout << A[i][j][4] << " ";
			printf("S:%f		", points[631][j][3]);
			printf("\n");
		}
	}*/

	__global__ void spline(double ***points, double ***cells, double ***res, double ***mtx_l, double **mtx_v, int size)
	{
		int threadID = blockDim.x * blockIdx.x + threadIdx.x;
		//if (threadID > size) return;

		if (threadID < size)
		{
			int p = points[threadID][0][3];

			//double *mtx_v = new double[p + 3];

			/*if (threadID == 512) {
				printf("end dec\n");
			}*/
			//

			//double a = 0.0;

			for (int i = 0; i < p; ++i)
			{
				for (int j = i + 1; j < p; ++j)
				{
					//double pt_x = points[threadID][i][0] - points[threadID][j][0];
					//double pt_z = points[threadID][i][2] - points[threadID][j][2];
					
					//double x = pt_x * pt_x;
					//double z = pt_z * pt_z;
					double elen = sqrt(((points[threadID][i][0] - points[threadID][j][0])*(points[threadID][i][0] - points[threadID][j][0])) + ((points[threadID][i][2] - points[threadID][j][2])*(points[threadID][i][2] - points[threadID][j][2])));
					if (elen == 0)
					{
						mtx_l[threadID][i][j] = mtx_l[threadID][j][i] = 0.0;
					}
					else
					{
						mtx_l[threadID][i][j] = mtx_l[threadID][j][i] = elen * elen * log(elen);
					}
					//a += elen * 2;

					if (i == 0 && j == 1)
					{
						//printf("%f", elen);
					}
				}
			}

			//printf("first for\n");

			//a /= (double)(p * p);
			//printf("a : %f\n", a);
			for (int i = 0; i < p; ++i)
			{
				// diagonal: reqularization parameters (lambda * a^2)
				/*if (i == 0) {
					std::cout << "sirayla i mtx_l(0,1)" << i << mtx_l(0, 1) << std::endl;
				}*/
				mtx_l[threadID][i][i] = 0.0;

				// P (p x 3, upper right)
				mtx_l[threadID][i][p + 0] = 1.0;
				mtx_l[threadID][i][p + 1] = points[threadID][i][0];
				mtx_l[threadID][i][p + 2] = points[threadID][i][2];

				// P transposed (3 x p, bottom left)
				mtx_l[threadID][p + 0][i] = 1.0;
				mtx_l[threadID][p + 1][i] = points[threadID][i][0];
				mtx_l[threadID][p + 2][i] = points[threadID][i][2];
			}
			//printf("second for\n");

			// O (3 x 3, lower right)
			for (int i = p; i < p + 3; ++i)
				for (int j = p; j < p + 3; ++j)
					mtx_l[threadID][i][j] = 0.0;

			//printf("third for\n");
			// Fill the right hand vector V
			for (int i = 0; i < p; ++i)
				mtx_v[threadID][i] = points[threadID][i][1];
			mtx_v[threadID][p + 0] = mtx_v[threadID][p + 1] = mtx_v[threadID][p + 2] = 0.0;
			//printf("fourth for\n");

			/*if (threadID == 512) {
				printf("mtx_v\n");
				for (int i = 0; i < p+3; i++) {
					printf("%f\n", mtx_v[i]);
				}
				printf("mtx_l\n");
				for (int i = 0; i < p+3; i++) {
					for (int j = 0; j < p+3; j++) {
						printf("%f	", mtx_l[i][j]);
					}
					printf("\n");
				}
			}*/

			int m = p + 3, n = p + 3;

			//int pivsign = 0;
			int *piv = new int[m];

			for (int i = 0; i < m; ++i)
				piv[i] = i;
			//pivsign = 1;

			for (int j = 0; j < n; ++j)
			{
				double *col = new double[m];
				for (int i = 0; i < m; i++)
				{
					col[i] = mtx_l[threadID][i][j];
				}
				double *row = new double[n];
				for (int i = 0; i < m; ++i)
				{
					for (int l = 0; l < n; l++)
					{
						row[l] = mtx_l[threadID][i][l];
					}
					int kmax = fminf(i, j);
					double s = 0.0;
					for (int k = 0; k < kmax; k++)
					{
						s += row[k] * col[k];
					}
					row[j] = col[i] -= s;
					for (int l = 0; l < m; l++)
					{
						mtx_l[threadID][l][j] = col[l];
					}
					for (int l = 0; l < m; l++)
					{
						mtx_l[threadID][i][l] = row[l];
					}
				}
				free(row);

				int p = j;
				for (int i = j + 1; i < m; i++)
				{
					if (fabs(col[i]) > fabs(col[p]))
					{
						p = i;
					}
				}
				free(col);


				if (p != j)
				{
					for (int k = 0; k < n; k++)
					{
						double t = mtx_l[threadID][p][k];
						mtx_l[threadID][p][k] = mtx_l[threadID][j][k];
						mtx_l[threadID][j][k] = t;
					}
					int k = piv[p];
					piv[p] = piv[j];
					piv[j] = k;
					//pivsign = -pivsign;
				}

				if (j < m && mtx_l[threadID][j][j] != 0.0)
				{
					for (int i = j + 1; i < m; i++)
					{
						mtx_l[threadID][i][j] /= mtx_l[threadID][j][j];
					}
				}
			}

			double y = 0;
			for (int i = 0; i < m; ++i)
			{
				if (piv[i] != i)
				{
					y = mtx_v[threadID][i];
					mtx_v[threadID][i] = mtx_v[threadID][piv[i]];
					mtx_v[threadID][piv[i]] = y;
				}
				for (int j = i; j < m; ++j)
					if (piv[j] == i)
					{
						piv[j] = piv[i];
						break;
					}
			}
			free(piv);
			
			for (int k = 0; k < n; k++)
			{
				for (int i = k + 1; i < n; i++)
				{
					mtx_v[threadID][i] -= mtx_v[threadID][k] * mtx_l[threadID][i][k];
				}
			}

			/*printf("mtx_v\n");
			for (int i = 0; i < m; i++) {
				printf("%f\n", mtx_v[i]);
			}*/

			for (int k = n - 1; k >= 0; k--)
			{
				mtx_v[threadID][k] /= mtx_l[threadID][k][k];
				//printf("1. %f\n", mtx_v[k]);
				for (int i = 0; i < k; i++)
				{
					mtx_v[threadID][i] -= mtx_v[threadID][k] * mtx_l[threadID][i][k];
				}
				//printf("2. %f\n", mtx_v[k]);
			}
			/*printf("mtx_v\n");
			for (int i = 0; i < m; i++) {
				printf("%f\n", mtx_v[i]);
			}*/

			int trn = cells[threadID][0][0];
			//printf("trn : %d\n", trn);
			for (int j = 0; j < trn; j++)
			{
				/*printf("mtx_v[p+0] : %f", mtx_v[p + 0]);
				printf("mtx_v[p + 1] : %f", mtx_v[p + 1]);
				printf("cells[threadID][j][1] : %f", cells[threadID][j][1]);
				printf("mtx_v[p + 2] : %f", mtx_v[p + 2]);
				printf("cells[threadID][j][2] : %f", cells[threadID][j][2]);*/

				double h = mtx_v[threadID][p + 0] + mtx_v[threadID][p + 1] * cells[threadID][j][1] + mtx_v[threadID][p + 2] * cells[threadID][j][2];
				//printf("global ici h : %f\n", h);
				double cx = cells[threadID][j][1];
				double cy = cells[threadID][j][2];

				for (int i = 0; i < p; i++)
				{
					double x = points[threadID][i][0] - cx;
					double z = points[threadID][i][2] - cy;

					double xx = x * x;
					double zz = z * z;
					double elen = sqrtf(xx + zz);
					//printf("elen : %f  ", elen);
					//printf("sqrtf ile : %f  ", sqrtf((x*x) + (z*z)));

					if (elen == 0.0)
					{
						h += mtx_v[threadID][i] * 0.0;
					}
					else
					{
						h += mtx_v[threadID][i] * (elen * elen * logf(elen));
					}
					//printf("res h: %f\n", h);
				}
				//printf("\nson res h: %f", h);
				res[threadID][j][0] = h;
				//printf("\nres[threadID][j][0] = %f", res[threadID][j][0]);
				/*if (threadID == 0)
				{
					//printf("\n%d. thread h: %f res[0] : %f \n", threadID, h, res[threadID][j][0]);
					//res[threadID + 1][j][0] = res[threadID][j][0];
					//printf("\n%d. thread h: %f res[631] : %f \n", threadID + 1, h, res[threadID + 1][j][0]);
				}*/
			}
			//__syncthreads();
			free(mtx_v);
		}
	}

	void source3::cells(double ***A, double ***cellsize, double ***h_res, int Asize, int celS){
		//allocate control_points
		std::cout << "allocating mainArray	" << std::endl;
		double ***h_c = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = A[i][0][3];
			h_c[i] = (double **)malloc(twoSize * sizeof(double *));
			for (int j = 0; j < twoSize; j++)
			{
				GPUerrchk(cudaMalloc((void **)&h_c[i][j], 4 * sizeof(double)));
				GPUerrchk(cudaMemcpy(h_c[i][j], A[i][j], 4 * sizeof(double), cudaMemcpyHostToDevice));
			}
		}
		double ***h_c1 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = A[i][0][3];
			GPUerrchk(cudaMalloc((void ***)&(h_c1[i]), twoSize * sizeof(double *)));
			GPUerrchk(cudaMemcpy(h_c1[i], h_c[i], twoSize * sizeof(double *), cudaMemcpyHostToDevice));
		}

		double ***d_c;
		GPUerrchk(cudaMalloc((void ****)&d_c, Asize * sizeof(double **)));
		GPUerrchk(cudaMemcpy(d_c, h_c1, Asize * sizeof(double **), cudaMemcpyHostToDevice));

		//allocate cells
		std::cout << "allocating cellsArray" << std::endl;
		double ***cells1 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = cellsize[i][0][0];
			cells1[i] = (double **)malloc(twoSize * sizeof(double *));
			for (int j = 0; j < twoSize; j++)
			{
				GPUerrchk(cudaMalloc((void **)&cells1[i][j], 3 * sizeof(double)));
				GPUerrchk(cudaMemcpy(cells1[i][j], cellsize[i][j], 3 * sizeof(double), cudaMemcpyHostToDevice));
			}
		}
		double ***cells = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = cellsize[i][0][0];
			GPUerrchk(cudaMalloc((void ***)&(cells[i]), twoSize * sizeof(double *)));
			GPUerrchk(cudaMemcpy(cells[i], cells1[i], twoSize * sizeof(double *), cudaMemcpyHostToDevice));
		}
		double ***d_cells;
		GPUerrchk(cudaMalloc((void ****)&d_cells, Asize * sizeof(double **)));
		GPUerrchk(cudaMemcpy(d_cells, cells, Asize * sizeof(double **), cudaMemcpyHostToDevice));

		//allocate res
		std::cout << "allocating resArray" << std::endl;
		double ***res1 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = cellsize[i][0][0];
			res1[i] = (double **)malloc(twoSize * sizeof(double *));
			for (int j = 0; j < twoSize; j++)
			{
				GPUerrchk(cudaMalloc((void **)&res1[i][j], 1 * sizeof(double)));
			}
		}
		double ***res2 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = cellsize[i][0][0];
			GPUerrchk(cudaMalloc((void ***)&res2[i], twoSize * sizeof(double *)));
			GPUerrchk(cudaMemcpy(res2[i], res1[i], twoSize * sizeof(double *), cudaMemcpyHostToDevice));
		}
		double ***d_res;
		GPUerrchk(cudaMalloc((void ****)&d_res, Asize * sizeof(double ***)));
		GPUerrchk(cudaMemcpy(d_res, res2, Asize * sizeof(double **), cudaMemcpyHostToDevice));

		std::cout << "allocating mtx_l" << std::endl;
		double ***mtx_l1 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = A[i][0][3] + 3;
			mtx_l1[i] = (double **)malloc(twoSize * sizeof(double *));
			for (int j = 0; j < twoSize; j++)
			{
				cudaMalloc((void **)&mtx_l1[i][j], twoSize * sizeof(double));
			}
		}
		double ***mtx_l2 = (double ***)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = A[i][0][3] + 3;
			cudaMalloc((void ***)&(mtx_l2[i]), twoSize * sizeof(double *));
			cudaMemcpy(mtx_l2[i], mtx_l1[i], twoSize * sizeof(double *), cudaMemcpyHostToDevice);
		}
		double ***d_mtx_l;
		cudaMalloc((void ****)&d_mtx_l, Asize * sizeof(double **));
		cudaMemcpy(d_mtx_l, mtx_l2, Asize * sizeof(double **), cudaMemcpyHostToDevice);

		std::cout << "allocating mtx_v" << std::endl;
		double **mtx_v1 = (double **)malloc(Asize * sizeof(double **));
		for (int i = 0; i < Asize; i++)
		{
			int twoSize = A[i][0][3];
			cudaMalloc(&mtx_v1[i], twoSize * sizeof(double));
		}
		double **d_mtx_v = (double **)malloc(Asize * sizeof(double **));
		cudaMalloc(&d_mtx_v, Asize * sizeof(double **));
		cudaMemcpy(d_mtx_v, mtx_v1, Asize * sizeof(double), cudaMemcpyHostToDevice);


		unsigned int numberOfThreads = Asize;
		unsigned int requiredNumberOfBlocks = (numberOfThreads / 1024) + 1;
		dim3 block = dim3(1024, 1, 1);
		dim3 grid = dim3(requiredNumberOfBlocks, 1, 1);
		std::clock_t start;
		double duration;
		start = std::clock();
		//cudaDeviceSetLimit(cudaLimitMallocHeapSize, 32*1024*1024);
		printf("launch kernel	");
		spline<<<grid, block>>>(d_c, d_cells, d_res, d_mtx_l, d_mtx_v, Asize); //---->>>>>><<<<<<<>>>>>>><<<<<<<>>>>>>><<<<<<>>>>
		printf("end kernel \n");
		cudaDeviceSynchronize();
		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "kernel process time :  " << duration << " seconds" << std::endl;
		std::cout << std::endl;


		printf("copy gpu to ram");
		for (int i = 0; i < Asize; i++)
		{
			int twoS = cellsize[i][0][0];
			for (int j = 0; j < twoS; j++)
			{
				cudaMemcpy(&h_res[i][j][0], res1[i][j], 1 * sizeof(double), cudaMemcpyDeviceToHost);
			}
		}

	}
}