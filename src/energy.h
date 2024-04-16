#pragma once
#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include "lsd.h"
#include <Eigen/Eigen>
#include "fline.h"
#include "tool.h"
using namespace std;
using namespace cv;
using namespace Eigen;

namespace energy
{
    SparseMatrix<double> shape_energy(vector<vector<Point > > mesh);
    SparseMatrix<double> line_energy(vector<vector<vector<Line > > > mesh_line,vector<vector<Point > > mesh,double*bins,int num);
    pair<SparseMatrix<double>,VectorXd > bound_energy(vector<vector<Point > > mesh,double inf,int n,int m);
    double* Line_rotate_count(int *bins,vector<vector<vector<Line > > > mesh_line,vector<vector<Point> >mesh,vector<vector<Point > > V,int number);
}