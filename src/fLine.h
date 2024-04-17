#pragma once
#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include "lsd.h"
#include "tool.h"
#include <Eigen/Eigen>

using namespace std;
using namespace cv;
using namespace Eigen;
#define PI acos(-1)

namespace fline
{
    void check_line(Mat img,vector<vector<vector<Line > > >line);

    double* line_img(Mat img,int *n_line,Mat&line_img,bool is_show=1);
    double* solve_img(Mat img);
    double get_gray(Mat&img,int x,int y);
    vector<vector<vector<Line > > > init_line(vector<vector<point > > mesh,double *line,int n,int*line_n);
    void init_bins(int *num,vector<vector<vector<Line > > >&mesh_line,int number);
}
