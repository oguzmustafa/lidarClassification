// Copyright 2009-2010 Green Code LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cmath>
#include <iostream>

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "DisjointRegions.h"
#include "IInterpolationRegion.h"
#include "IPointVector.h"
#include "IRegionGenerator.h"
#include "LineIndent.h"
#include "Point.h"
#include "ProgressBar.h"
#include "RasterSurface.h"
#include "RegularizedSpline.h"
#include "SplineExceptions.h"
#include "SurfaceInterpolation.h"
#include "XYCoordinates.h"
#include <source.cuh>
#include <source2.cuh>
#include <source3.cuh>




namespace mcc
{
  // A point selector that selects every point.
  bool useEveryPoint(const IPoint & point)
  {
    return true;
  }

  //---------------------------------------------------------------------------

  SurfaceInterpolation::SurfaceInterpolation()
    : prevCellResolution_(0)
  {
  }

  //---------------------------------------------------------------------------

  boost::shared_ptr<IRasterSurface> SurfaceInterpolation::operator()(const IPointVector & points,
                                                                     double               cellResolution,
                                                                     double               tension)
  {
    return this->operator()(points, &useEveryPoint, cellResolution, tension);
  }

  //---------------------------------------------------------------------------
  
  boost::shared_ptr<IRasterSurface> SurfaceInterpolation::operator()(const IPointVector & points,
                                                                     PointSelector        pointSelector,
                                                                     double               cellResolution,
                                                                     double               tension)
  {
    // If the cell resolution is different from the previous call, then
    // create a new raster based on the new cell size.
    if (cellResolution != prevCellResolution_) {
      // Determine the desired raster dimensions by adding a 1/2 cell-wide
      // margin around the boundaries of the point cloud read by the program.
      double margin = cellResolution / 2;
      double desiredWidth  = (inputExtent_.maxX - inputExtent_.minX) + 2 * margin; // left & right margins
      double desiredHeight = (inputExtent_.maxY - inputExtent_.minY) + 2 * margin; // top & bottom margins

      // Compute the numer of rows & columns needed to minimally cover the
      // desired raster dimensions.
      unsigned int cols = int(std::ceil(desiredWidth  / cellResolution));
      unsigned int rows = int(std::ceil(desiredHeight / cellResolution));

      // Determine the lower-left corner of the raster by centering its actual
      // area around the original input extent.
      double actualWidth = cols * cellResolution;
      double horizontalMargin = (actualWidth - desiredWidth) / 2;
      Coordinate x0 = Coordinate(inputExtent_.minX - horizontalMargin);

      double actualHeight = rows * cellResolution;
      double verticalMargin = (actualHeight - desiredHeight) / 2;
      Coordinate y0 = Coordinate(inputExtent_.minY - verticalMargin);

      XYCoordinates lowerLeft(x0, y0);

      rasterSurface_ = boost::make_shared<RasterSurface>(rows, cols, lowerLeft, Coordinate(cellResolution));
      prevCellResolution_ = cellResolution;
    }


    // Determine where splines will be interpolated for the points and the
    // raster.
	//int b = 0;
    /*boost::shared_ptr<IRegionGenerator> reg = boost::make_shared<DisjointRegions>();
	int k = 0;
	double *arr = new double();

	while (k<10) {
		const IInterpolationRegion * regg = reg->getNextRegion();
		std::vector<const IPoint *> points = regg->points();
		arr[k] = points.operator[](k)->x();
		k++;
	}
	while (k < 20) {
		printf("%f\n", arr[k-10]);
		k++;
	}*/
	boost::shared_ptr<IRegionGenerator> regions = boost::make_shared<DisjointRegions>();
    int nRegions = regions->subdivide(points, pointSelector, *rasterSurface_);
	//printf("n Regions packet size : %d\n", nRegions);
    LineIndent indent("  ");

    // For each region, compute a spline for its points and then interpolate
    // heights for its cells.
    //std::cout << indent << "Computing splines for regions and cell heights for raster surface:" << std::endl
    //          << indent << "  ";
    //ProgressBar progressBar(std::cout, nRegions);
	std::cout << std::endl;
    
	double*** A = new double**[nRegions];
	double*** cellsize = new double**[nRegions];
	double*** res = new double**[nRegions];
	
	source x;
	int nSplinesComputed = 0, a = 0, i = 0, counter = 0;
    while (const IInterpolationRegion * region = regions->getNextRegion()) {
		//if (a>0||a<10){
			//printf("SurfaceInterpolation while giris");//---------->>>>>>
			//std::cout << region->points().size() << std::endl;
		//
		//printf("region points size: %d\n", region->points().size());
			//a++;
		//}
		
		std::vector<const IPoint *> points = region->points();
		
		A[i] = new double*[points.size()];
		for (int j = 0; j < points.size(); j++) {
			A[i][j] = new double[4];
			A[i][j][0] = 0;
			A[i][j][1] = 0;
			A[i][j][2] = 0;
			A[i][j][3] = 0;
		}
		
		for (int j = 0; j < points.size(); j++)
		{
			A[i][j][0] = points.operator[](j)->x();
			A[i][j][1] = points.operator[](j)->z();
			A[i][j][2] = points.operator[](j)->y();
			A[i][j][3] = points.size();
		}
		//printf("cel size : %d", region->cells().size());
		int cels = region->cells().size();

		cellsize[i] = new double*[cels];
		for (int j = 0; j < cels; j++) {
			//printf("for ici cel size : %d", region->cells().size());
			cellsize[i][j] = new double[3];
			cellsize[i][j][0] = 0;
			cellsize[i][j][1] = 0;
			cellsize[i][j][2] = 0;
		}
		//cellsize[0][0][0] = 5.0;
		//printf("while disi cellsize[i][0][0] : %f\n ", cellsize[0][0][0]);


			

		
	
		/*for (int j = 0; j < points.size(); j++)
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
		}*/
	


		//std::cout << "points deger x: " << points.operator[](a)->x();
		/*printf("sonraki next x: %f", points.operator[](0)->x());
		printf(" y: %f", points.operator[](0)->y());
		printf(" z: %f\n", points.operator[](0)->z());
		*/
		//std::cout << "Surface points size: " << points.size() << std::endl;
		//("%f\n", points.size());
		double abcd = 0.1;
		bool splineComputed = false;
		int j = 0;
		while (!splineComputed) {
			try {
				/*if (counter == 0) {
					RegularizedSpline spline(region->points(), 0.0);
					BOOST_FOREACH(const Cell & cell, region->cells()) {
						abcd = spline.interpolateHeight(cell.x(), cell.y());
						//printf("interpolation h : %f", abcd);
						(*rasterSurface_)[cell] = abcd;
					}
					counter++;
				}*/
				//RegularizedSpline spline(region->points(), 0.0);
				//std::cout << region->points() << std::endl;
				splineComputed = true;
				//std::cout << "dongu ustu cell size : " << region->cells().size() << std::endl;
				//int cels = region->cells().size();
				cellsize[i][j][0] = cels;
				//printf("cellsize [i][j][0] : %f", cellsize[i][j][0]);
				BOOST_FOREACH(const Cell & cell, region->cells()) {
					cellsize[i][j][1] = cell.x();
					cellsize[i][j][2] = cell.y();
					/*abcd = spline.interpolateHeight(cell.x(), cell.y());
					if (counter == 1) {
						printf("interpolation h : %f", abcd);
					}
					(*rasterSurface_)[cell] = abcd;*/
					counter++;
					j++;
				}
			}
			catch (SingularMatrixException) {
				// Add another neighboring point and try the spline calculation again.
				//printf("exceptiona girdi...--->>>>\n");
				regions->addNeighborPointsToCurrentRegion(1);
				// A safety check to prevent an endless loop from consuming all the
				// point cloud.
				if (region->points().size() >= 300)
					throw;  // Bail
			}
		}
		//std::cout << "dongu sonu cell size : "<<counter << std::endl;
		//counter = 0;
		a++;
		i++;
		nSplinesComputed++;
		//progressBar.update(nSplinesComputed);
	}
	int celS = 0;
	for (int i = 0; i < nRegions; i++) {
		int twoSize = cellsize[i][0][0];
		celS += twoSize;
		res[i] = new double *[twoSize];
		for (int j = 0; j < twoSize; j++) {
			res[i][j] = new double[1];
		}
	}

	printf("cells size : %d\n", celS);
	/*printf("surfaceint end for\n");
	for (int i = 0; i < nRegions; i++) {
		int twoSize = cellsize[i][0][0];
		for (int j = 0; j < twoSize; j++) {
			printf("%d. cell %d.pack x: %f y: %f	", i, j, cellsize[i][j][1], cellsize[i][j][2]);
		}
		printf("\n");
		int ts = A[i][0][3];
		for (int j = 0; j < ts; j++) {
			printf("%f	%f	%f	%f\n", A[i][j][0], A[i][j][1], A[i][j][2], A[i][j][3]);
		}
	}*/


	source3 ab;
	ab.cells(A, cellsize, res, nRegions, celS);

	//printf("begin spline \n");
	int regCo = 0;
	boost::shared_ptr<IRegionGenerator> regis = boost::make_shared<DisjointRegions>();
	int n = regis->subdivide(points, pointSelector, *rasterSurface_);
	while (const IInterpolationRegion * regi = regis->getNextRegion()) {
		bool splineComputed = false;
		while (!splineComputed) {
			splineComputed = true;
			int j = 0;
			BOOST_FOREACH(const Cell & cell, regi->cells()) {
				//if(regCo < 10)
					//printf("res %d : %f", regCo, res[regCo][j][0]);
				(*rasterSurface_)[cell] = res[regCo][j][0];
				j++;
			}
		}
		regCo++;
	}
	//printf("regCo: %d", regCo);



	std::cout << std::endl;
	//x.test(A,nRegions);
	//double ***XAA;
	//int s = 999;
	//x.test(XAA, s);
	//printf("while cikis a: %d, i: %d", a, i);
	//std::cout << std::endl;

	//if (d == 0) {
		return rasterSurface_;
	//}
	
  }
}
