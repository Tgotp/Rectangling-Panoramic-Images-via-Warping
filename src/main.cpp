#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "seam_carving.h"
#include "mesh_placement.h"
#include "energy.h"
#include "fLine.h"
#include <Eigen/Eigen>
#include "GLout.h"
using namespace cv;
using namespace std;


int main()
{
    double st = clock();
    Mat img = imread("input/3_input.jpg");
    namedWindow("origin image",WINDOW_FREERATIO);
    imshow("origin image", img);
    SeamC::init(img,0); //op
    img = SeamC::get_rec_img();

    double st1 = clock();
    cout << fixed << setprecision(3)<< "seam carving cost time: " << (st1 - st) / CLOCKS_PER_SEC << endl; 

    vector<vector<Point > > U = SeamC::get_U();
    // waitKey(0);
    vector<vector<Point > > mesh = Mesh_Placement(SeamC::get_seam_carving(),img,U,0); // 网格坐标 op
    // double st2 = clock();
    // cout << fixed << setprecision(3)<< "mesh cost time: " << (st2 - st1) / CLOCKS_PER_SEC << endl; 

    // int *n_line = new int; // 线段数量
    // double* line = fline::line_img(img,n_line,0); // 线段 and op
    // vector<vector<vector<Line > > > mesh_line = fline::init_line(mesh,line,*n_line); // 每个网格中的线段
    // // fline::check_line(img,mesh_line);// op

    // int M = 50;
    // int*num;
    // num = new int[M];
    // double*bins = fline::Line_rotate_count(num,mesh_line,M);

    // for(int iter = 0;iter < 10; ++ iter)
    // {
    //     double EL = energy::line_energy(mesh_line,mesh); 
    // }
    // // energy::shape_energy(mesh,energy::V);
    // double st3 = clock();
    // cout << fixed << setprecision(3)<< "energy cost time: " << (st3 - st2) / CLOCKS_PER_SEC << endl; 

    GLout(img,mesh);
    waitKey(0);
    return 0;
}