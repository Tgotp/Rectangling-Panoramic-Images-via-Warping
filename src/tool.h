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

struct point
{
    double x,y;
    point():x(0),y(0){};
    point(double _x,double _y) { x = _x,y = _y; }

    point operator + (point a)
    {
        return point(x + a.x,y + a.y);
    }
    point operator - (point a)
    {
        return point(x - a.x,y - a.y);
    }

    double len()
    {
        return sqrt(x * x + y * y);
    }

    bool operator == (point a)
    {
        return point(x - a.x,y - a.y).len() < 0.00001;
    }
 
    friend ostream &operator<<( ostream &output ,const point&z)
    { 
            output << "x : " << z.x << " y : " << z.y << endl;
            return output;            
    }
};

struct Line
{
    point x,y; int pos;double theta;bool solve;
    Line(point _x,point _y,int p = 0)
    {
        x = _x,y = _y; pos = p;
    }
};

void output(MatrixXd x,string name = "No name mentioned");
void output(SparseMatrix<double> x,string name = "No name mentioned");
SparseMatrix<double> merge_matrix(SparseMatrix<double> a,SparseMatrix<double> b);
VectorXd merge_matrix(VectorXd a,VectorXd b);
void show_new_line(MatrixXd K,vector<vector<Point > > V,int i,int j);