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

namespace energy
{
    double shape_energy(vector<vector<Point > > mesh,vector<vector<Point > > V);
    double line_energy(vector<vector<vector<Line > > > mesh,vector<vector<Point > > V);
}