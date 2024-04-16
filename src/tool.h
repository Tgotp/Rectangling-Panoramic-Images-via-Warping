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
using namespace std;
using namespace cv;
using namespace Eigen;

void output(MatrixXd x,string name = "No name mentioned");
void output(SparseMatrix<double> x,string name = "No name mentioned");
SparseMatrix<double> merge_matrix(SparseMatrix<double> a,SparseMatrix<double> b);
VectorXd merge_matrix(VectorXd a,VectorXd b);