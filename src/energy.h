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
    SparseMatrix<double> shape_energy(vector<vector<point > > mesh);
    SparseMatrix<double> line_energy(vector<vector<vector<Line > > > mesh_line,vector<vector<point > > mesh,vector<vector<point > > V,double*bins,int num,int n,int m);
    pair<SparseMatrix<double>,VectorXd > bound_energy(vector<vector<point > > mesh,double inf,int n,int m);
    void Line_rotate_count(double *bins,int *num,vector<vector<vector<Line > > >&mesh_line,vector<vector<point> >mesh,vector<vector<point > > V,int number);
}