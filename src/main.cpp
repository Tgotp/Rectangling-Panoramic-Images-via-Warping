#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "seam_carving.h"
#include "mesh_placement.h"
using namespace cv;
using namespace std;


int main()
{
    double st = clock();
    Mat img = imread("input//2_input.jpg");
    namedWindow("origin image",WINDOW_FREERATIO);
    imshow("origin image", img);
    SeamC::init(img,0);
    img = SeamC::get_seam_carving();
    imshow("carvd image", img);
    double st1 = clock();
    cout << fixed << setprecision(3)<< "seam carving cost time: " << (st1 - st) / CLOCKS_PER_SEC << endl; 

    vector<vector<pair<int,int> > > U = SeamC::get_U();
    // waitKey(0);
    Mat x = Mesh_Placement(img,SeamC::get_rec_img(),U,0);
    double st2 = clock();
    cout << fixed << setprecision(3)<< "mesh cost time: " << (st2 - st1) / CLOCKS_PER_SEC << endl; 

    

    waitKey(0);
    return 0;
}