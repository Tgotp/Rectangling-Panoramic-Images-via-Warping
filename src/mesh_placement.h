#pragma once
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
#include "tool.h"

using namespace std;
using namespace cv;
using namespace Eigen;

vector<vector<point > > Mesh_Placement(Mat img,Mat image,vector<vector<point > > U,vector<vector<point > >&pos,bool is_show=1);
VectorXd get_mesh_point(vector<vector<point> >pos);
vector<vector<point > > vec_to_mesh(VectorXd V);
SparseMatrix<double> merge_matrix(SparseMatrix<double> a,SparseMatrix<double> b);
VectorXd merge_matrix(VectorXd a,VectorXd b);
void img_mesh(Mat&img,vector<vector<point > > pos);
