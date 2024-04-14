#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/opencv.hpp>
#include "lsd.h"
#include <Eigen>
using namespace std;
using namespace cv;
using namespace Eigen;

namespace energy
{
    double* line_img(Mat img,bool is_show=1);
    double* solve_img(Mat img);
    double get_gray(Mat&img,int x,int y);

    MatrixXf shape_energy(Mat img,vector<vector<pair<int,int> > > mesh);
}