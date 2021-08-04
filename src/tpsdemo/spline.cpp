/*
 *  Copyright (C) 2003,2005 by Jarno Elonen
 *
 *  TPSDemo is Free Software / Open Source with a very permissive
 *  license:
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

#include "spline.h"
#include "ludecomposition.h"

#include <vector>
#include <cmath>

using namespace boost::numeric::ublas;

//-----------------------------------------------------------------------------

static double tps_base_func(double r)
{
  if ( r == 0.0 )
    return 0.0;
  else
    return r*r * log(r);
}

//-----------------------------------------------------------------------------

/*
 *  Calculate Thin Plate Spline (TPS) weights from
 *  control points.
 */
int abc = 0;
tpsdemo::Spline::Spline(const std::vector<Vec> & control_pts, double regularization)
  : p(control_pts.size()),
    control_points(control_pts),
    mtx_v(p+3, 1),
    mtx_orig_k(p, p)
{
	/*std::cout << "len: " << control_points[0].len() <<
		"	norm: " << control_points[0].norm() <<
		"	x: " << control_points[0].x <<
		"	y: " << control_points[0].y <<
		"	z: " << control_points[0].z << std::endl;*/
	/*std::cout << "spline x and y" << std::endl;
	for (int i = 0; i < p; i++) {
		std::cout << control_points[i].x << "	" << control_points[i].z << std::endl;
	}*/

  // You We need at least 3 points to define a plane
  if ( control_points.size() < 3 )
    throw std::runtime_error("need at least 3 points for thin plate spline");

  //unsigned p = control_points.size();

  // Allocate the matrix and vector
  matrix<double> mtx_l(p+3, p+3);
  //matrix<double> mtx_v(p+3, 1);
  //matrix<double> mtx_orig_k(p, p);

  // Fill K (p x p, upper left of L) and calculate
  // mean edge length from control points
  //
  // K is symmetrical so we really have to
  // calculate only about half of the coefficients.
  double a = 0.0;
  for ( unsigned i=0; i<p; ++i )
  {
    for ( unsigned j=i+1; j<p; ++j )
    {
      Vec pt_i = control_points[i];
      Vec pt_j = control_points[j];
      pt_i.y = pt_j.y = 0;
      double elen = (pt_i - pt_j).len();
	  //std::cout << " spline icinde 85.satir elen: " <<elen<< std::endl;
	  Vec test = (pt_i - pt_j);
	  test.x = test.x*test.x;
	  test.z = test.z*test.z;
	  double testt = sqrt(test.x+test.z);



      mtx_l(i,j) = mtx_l(j,i) =
        mtx_orig_k(i,j) = mtx_orig_k(j,i) =
          tps_base_func(elen);
      a += elen * 2; // same for upper & lower tri
	  /*if (i == 0 && j == 1) {
		  printf("%f", elen);
		  std::cout << "testt len: " << testt << std::endl;
	  }*/
    }
  }
  /*for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < p + 3; j++) {
		  printf("%f	", mtx_l(i,j));
	  }
	  printf("\n");
  }*/

  a /= (double)(p*p);

  // Fill the rest of L
  for ( unsigned i=0; i<p; ++i )
  {
    // diagonal: reqularization parameters (lambda * a^2)
	  /*if (i == 0) {
		  std::cout << "sirayla i mtx_l(0,1)" << i << mtx_l(0, 1) << std::endl;
	  }*/
    mtx_l(i,i) = mtx_orig_k(i,i) =
      0.0 * (a*a);

    // P (p x 3, upper right)
    mtx_l(i, p+0) = 1.0;
    mtx_l(i, p+1) = control_points[i].x;
    mtx_l(i, p+2) = control_points[i].z;

    // P transposed (3 x p, bottom left)
    mtx_l(p+0, i) = 1.0;
    mtx_l(p+1, i) = control_points[i].x;
    mtx_l(p+2, i) = control_points[i].z;
  }
  // O (3 x 3, lower right)
  for ( unsigned i=p; i<p+3; ++i )
    for ( unsigned j=p; j<p+3; ++j )
      mtx_l(i,j) = 0.0;


  // Fill the right hand vector V
  for ( unsigned i=0; i<p; ++i )
    mtx_v(i,0) = control_points[i].y;
  mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;
  /*printf("mtx_l\n");
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < p + 3; j++) {
		  printf("%f	", mtx_l(i,j));
	  }
	  printf("\n");
  }*/

  /*std::cout << "mtx_l" << std::endl;
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < p+3; j++) {
		  std::cout << mtx_l(i, j) <<"	";
	  }
	  std::cout << std::endl;
  }*/
  /*std::cout << "mtx_v" << std::endl;
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < 1; j++) {
		  std::cout << mtx_v(i, j) << "	";
	  }
	  std::cout << std::endl;
  }*/
  /*std::cout << "before lu solve mtx_v" << std::endl;
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < 1; j++) {
		  std::cout << mtx_v(i,j) << std::endl;
	  }
  }*/

  // Solve the linear system "inplace"

  double **A = new double*[p+3];
  for (int i = 0; i < p+3; i++) {
	  A[i] = new double[p + 3];
  }
  for (int i = 0; i < p+3; i++) {
	  for (int j = 0; j < p+3; j++) {
		  A[i][j] = mtx_l(i, j);
	  }
  }
  /*printf("A spline \n");
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < p + 3; j++) {
		  printf("%f	", A[i][j]);
	  }
	  printf("\n");
  }*/


  int m = p+3, n = p+3;

  int pivsign = 0;
  int* piv = new int[m];

  for (int i = 0; i < m; ++i)
	  piv[i] = i;
  pivsign = 1;

  for (int j = 0; j < n; ++j) {
	  double *col = new double[m];
	  for (int i = 0; i < m; i++) {
		  col[i] = A[i][j];
	  }
	  double *row = new double[n];
	  for (int i = 0; i < m; ++i) {
		  for (int l = 0; l < n; l++) {
			  row[l] = A[i][l];
		  }
		  int kmax = fminf(i, j);
		  double s = 0.0;
		  for (int k = 0; k < kmax; k++) {
			  s += row[k] * col[k];
		  }
		  A[i][j] = col[i] -= s;
	  }

	  int p = j;
	  for (int i = j + 1; i < m; i++) {
		  if (fabs(col[i]) > fabs(col[p])) {
			  p = i;
		  }
	  }

	  //printf("p : %d",p);

	  if (p != j) {
		  for (int k = 0; k < n; k++) {
			  double t = A[p][k];
			  A[p][k] = A[j][k];
			  A[j][k] = t;
		  }
		  int k = piv[p];
		  //printf("pp : %d",piv[p]);
		  piv[p] = piv[j];
		  //printf("pj : %d",piv[j]);
		  piv[j] = k;
		  //printf("k : %d",k);
		  pivsign = -pivsign;
	  }

	  /*for(int i=0; i<3; i++) {
		  printf("piv%d: %d ",i,piv[i]);
	  }*/

	  if (j < m && A[j][j] != 0.0) {
		  for (int i = j + 1; i < m; i++) {
			  A[i][j] /= A[j][j];
		  }
	  }
  }
  /*printf("A lude\n");
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < p + 3; j++) {
		  printf("%f	", A[i][j]);
	  }
	  printf("\n");
  }*/


  if (0 != LU_Solve(mtx_l, mtx_v))
  {
    throw SingularMatrixError();
  }
  /*std::cout << "after lu solve mtx_v" << std::endl;
  for (int i = 0; i < p + 3; i++) {
	  for (int j = 0; j < 1; j++) {
		  std::cout << mtx_v(i, j) << std::endl;
	  }
  }*/
}

//-----------------------------------------------------------------------------

double tpsdemo::Spline::interpolate_height(double x, double z) const
{

  double h = mtx_v(p+0, 0) + mtx_v(p+1, 0)*x + mtx_v(p+2, 0)*z;
  //std::cout << "height mtx_v(0,0): " << mtx_v(0, 0) << std::endl;
  //printf("host ici h : %f\n", h);
  Vec pt_i, pt_cur(x,0,z);
  /*std::cout << "pt_cur x y z " << pt_cur.x << std::endl;
  std::cout << pt_cur.y << std::endl;
  std::cout << pt_cur.z << std::endl;*/

  for ( unsigned i=0; i<p; ++i )
  {
    pt_i = control_points[i];
    pt_i.y = 0;
	double elen = (pt_i - pt_cur).len();
	/*if (abc > 0 || abc < 3) {
		std::cout << elen << std::endl;
	}*/
    h += mtx_v(i,0) * tps_base_func( ( pt_i - pt_cur ).len());
  }
  return h;
}

//-----------------------------------------------------------------------------

double tpsdemo::Spline::compute_bending_energy() const
{
  matrix<double> w( p, 1 );
  for ( unsigned i=0; i<p; ++i )
    w(i,0) = mtx_v(i,0);
  matrix<double> be = prod( prod<matrix<double> >( trans(w), mtx_orig_k ), w );
  return be(0,0);
}
