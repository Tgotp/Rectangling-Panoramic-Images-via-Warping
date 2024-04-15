#pragma once
#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include "lsd.h"
#include <Eigen/Eigen>
using namespace std;
using namespace cv;
using namespace Eigen;
#define PI acos(-1)

struct Line
{
    Point x,y; double theta; int pos;
    Line(Point _x,Point _y,double p = 0,double _k = 0)
    {
        x = _x,y = _y; theta = p; pos = _k;
    }
};

namespace fline
{
    void check_line(Mat img,vector<vector<vector<Line > > >line);

    double* line_img(Mat img,int *n_line,bool is_show=1);
    double* solve_img(Mat img);
    double get_gray(Mat&img,int x,int y);
    vector<vector<vector<Line > > > init_line(vector<vector<Point > > mesh,double *line,int n,int*line_n);

}
