#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H
/*
 *  Copyright (C) 2003, 2004 by Jarno Elonen
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// Solve a linear equation system a*x=b using inplace LU decomposition.
//
// Stores x in 'b' and overwrites 'a' (with a pivotted LUD).
//
// Matrix 'b' may have any (>0) number of columns but
// must contain as many rows as 'a'.
//
// Possible return values:
//  0=success
//  1=singular matrix
//  2=a.rows != b.rows
template <typename T>
int LU_Solve(
	boost::numeric::ublas::matrix<T> &a,
	boost::numeric::ublas::matrix<T> &b)
{
	// This routine is originally based on the public domain draft for JAMA,
	// Java matrix package available at http://math.nist.gov/javanumerics/jama/

	typedef boost::numeric::ublas::matrix<T> Matrix;
	typedef boost::numeric::ublas::matrix_row<Matrix> Matrix_Row;
	typedef boost::numeric::ublas::matrix_column<Matrix> Matrix_Col;

	if (a.size1() != b.size1())
		return 2;

	int m = a.size1(), n = a.size2();
	//std::cout << "mtx_l size 1: " << m << "mtx_l size 2: " << n << std::endl;
	int pivsign = 0;
	int *piv = (int *)alloca(sizeof(int) * m);

	// PART 1: DECOMPOSITION
	//
	// For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
	// unit lower triangular matrix L, an n-by-n upper triangular matrix U,
	// and a permutation vector piv of length m so that A(piv,:) = L*U.
	// If m < n, then L is m-by-m and U is m-by-n.
	{
		// Use a "left-looking", dot-product, Crout/Doolittle algorithm.
		for (int i = 0; i < m; ++i)
			piv[i] = i;
		pivsign = 1;

		// Outer loop.
		for (int j = 0; j < n; ++j)
		{
			// Make a copy of the j-th column to localize references.
			//Matrix_Col LUcolj(a, j);

			// Apply previous transformations.
			/*for (int i = 0; i < m; ++i)
		{
			Matrix_Row LUrowi(a,i);

			// This dot product is very expensive.
			// Optimize for SSE2?
			int kmax = (i<=j)?i:j;
			typename Matrix_Row::const_iterator ri_ite( LUrowi.begin());
			typename Matrix_Col::const_iterator cj_ite( LUcolj.begin());
			typename Matrix::value_type sum = 0.0;
			while( kmax-- > 0 )
			  sum += (*(ri_ite++)) * (*(cj_ite++));
			LUrowi[j] = LUcolj[i] -= sum;
		}*/

			double *col = new double[m];
			for (int i = 0; i < m; i++)
			{
				col[i] = a(i, j);
			}
			double *row = new double[n];
			for (int i = 0; i < m; ++i)
			{
				for (int l = 0; l < n; l++)
				{
					row[l] = a(i, l);
				}

				int kmax = std::min(i, j);
				double s = 0.0;
				for (int k = 0; k < kmax; k++)
				{
					s += row[k] * col[k];
				}
				row[j] = col[i] -= s;
				for (int i = 0; i < m; i++)
				{
					a(i, j) = col[i];
				}
				for (int l = 0; l < m; l++)
				{
					a(i, l) = row[l];
				}
			}

			// Find pivot and exchange if necessary.
			//
			// Slightly optimized version of:
			//  for (int i = j+1; i < m; ++i)
			//    if ( fabs(LUcolj[i]) > fabs(LUcolj[p]) )
			//      p = i;
			/*int p = j;
		typename Matrix::value_type coljp_abs = fabs(LUcolj[p]);
		for ( typename Matrix_Col::const_iterator
				beg = LUcolj.begin(),
				ite = beg + j+1,
				end = LUcolj.end();
			  ite < end;
			  ++ite )
		{
		  if (fabs(*ite) > coljp_abs)
		  {
			p = ite-beg;
			coljp_abs = fabs(LUcolj[p]);
		  }
		}*/
			int p = j;
			for (int i = j + 1; i < m; i++)
			{
				if (std::abs(col[i]) > std::abs(col[p]))
				{
					p = i;
				}
			}

			/*if (p != j)
		{
		  Matrix_Row raj(a, j);
		  Matrix_Row(a, p).swap(raj);

		  int tmp = piv[p];
		  piv[p] = piv[j];
		  piv[j] = tmp;
		  pivsign = -pivsign;
		}*/
			if (p != j)
			{
				for (int k = 0; k < n; k++)
				{
					double t = a(p, k);
					a(p, k) = a(j, k);
					a(j, k) = t;
				}
				int k = piv[p];
				piv[p] = piv[j];
				piv[j] = k;
				pivsign = -pivsign;
			}

			// Compute multipliers.
			/*if (j < m && a(j, j) != 0.0)
		  for (int i = j + 1; i < m; ++i)
			LUcolj[i] /= LUcolj[j];*/

			if (j < m && a(j, j) != 0.0)
				for (int i = j + 1; i < m; i++)
					a(i, j) /= a(j, j);

			/*for (int i = 0; i < m; i++)
		{
		  for (int j = 0; j < m; j++)
		  {
			printf("%f	", a(i, j));
		  }
		  printf("\n");
		}*/
		}
	}

	/*std::cout << "c matrix" << std::endl;
  for (int i = 0; i < m; i++) {
	  for (int j = 0; j < n; j++) {
		  std::cout << c(i, j) << "    ";
	  }
	  std::cout << std::endl;
  }*/

	/*std::cout << "a matrix" << std::endl;
  for (int i = 0; i < m; i++) {
	  for (int j = 0; j < m; j++) {
		  printf("%f	", a(i, j));
	  }
	  printf("\n");
  }*/
	/*printf("piv\n");
  for (int i = 0; i < m; i++) {
	  printf("%d\n", piv[i]);
  }*/

	// PART 2: SOLVE

	// Check singluarity
	for (int j = 0; j < n; ++j)
		if (a(j, j) == 0)
			return 1;

	double *z = new double[m];
	for (int i = 0; i < m; i++)
	{
		//printf("%f\n", b(i));
		z[i] = b(i);
	}
	// Reorder b according to pivotting
	for (int i = 0; i < m; ++i)
	{
		if (piv[i] != i)
		{
			Matrix_Row b_ri(b, i);
			Matrix_Row(b, piv[i]).swap(b_ri);
			for (int j = i; j < m; ++j)
				if (piv[j] == i)
				{
					piv[j] = piv[i];
					break;
				}
		}
	}
	/*printf("b\n");
  for (int i = 0; i < m; i++) {
	  printf("%f\n", b(i));
  }*/
	/*printf("b\n");
  for (int i = 0; i < m; i++) {
	  printf("%f\n", b(i));
  }
  printf("piv\n");
  for (int i = 0; i < m; i++) {
	  printf("%d\n", piv[i]);
  }*/
	/*for (int i = 0; i < m; i++) {
	  printf("%f\n", z[i]);
  }*/
	/*for (int i = 0; i < m; ++i) {
	  if (piv[i] != i) {
		  y = z[i];
		  z[i] = z[piv[i]];
		  z[piv[i]] = y;
	  }
	  for(int j=i;j<m;++j)
		  if (piv[j] == i) {
			  piv[j] = piv[i];
			  break;
		  }
  }
 */
	// Solve L*Y = B(piv,:)
	for (int i = 0; i < m; i++)
	{
		z[i] = b(i);
	}
	for (int k = 0; k < n; ++k)
	{
		const Matrix_Row &b_rk = Matrix_Row(b, k);
		for (int i = k + 1; i < n; ++i)
		{
			const typename Matrix_Row::value_type aik = a(i, k);
			Matrix_Row(b, i) -= b_rk * aik;
		}
	}

	/*printf("b \n");
	for (int i = 0; i < m; i++)
	{
		printf("%f\n", z[i]);
	}
	printf("a \n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			printf("%f	", a(i, j));
		}
		printf("\n");
	}*/


	for (int k = 0; k < n; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			z[i] -= z[k] * a(i, k);
		}
	}

	/*printf("b \n");
	for (int i = 0; i < m; i++)
	{
		printf("%f\n", z[i]);
	}*/

	// Solve U*X = Y;
	for (int k = n - 1; k >= 0; --k)
	{
		Matrix_Row(b, k) *= 1.0 / a(k, k);

		const Matrix_Row &b_rk = Matrix_Row(b, k);
		for (int i = 0; i < k; ++i)
		{
			const typename Matrix_Row::value_type aik = a(i, k);
			Matrix_Row(b, i) -= b_rk * aik;
		}
	}

	for (int k = n - 1; k >= 0; k--)
	{
		z[k] /= a(k, k);
		for (int i = 0; i < k; i++)
		{
			z[i] -= z[k] * a(i, k);
		}
	}

	/*printf("b \n");
	for (int i = 0; i < m; i++)
	{
		printf("%f\n", z[i]);
	}*/

	return 0;
}

#endif // LUDECOMPOSITION_H
