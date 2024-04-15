#pragma once
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include <Eigen/Eigen>
using namespace std;
using namespace cv;
using namespace Eigen;

vector<vector<Point > > Mesh_Placement(Mat img,Mat image,vector<vector<Point > > U,vector<vector<Point > >&pos,bool is_show=1);
VectorXd get_mesh_point(vector<vector<Point> >pos);
vector<vector<Point > > vec_to_mesh(VectorXd V);
SparseMatrix<double> merge_matrix(SparseMatrix<double> a,SparseMatrix<double> b);
VectorXd merge_matrix(VectorXd a,VectorXd b);