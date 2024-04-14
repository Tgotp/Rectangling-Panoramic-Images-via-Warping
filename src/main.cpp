#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "seam_carving.h"
#include "mesh_placement.h"
#include "energy.h"
#include <Eigen>
using namespace cv;
using namespace std;


int main()
{
    double st = clock();
    Mat img = imread("input//3_input.jpg");
    namedWindow("origin image",WINDOW_FREERATIO);
    imshow("origin image", img);
    SeamC::init(img,0);
    img = SeamC::get_rec_img();

    double st1 = clock();
    cout << fixed << setprecision(3)<< "seam carving cost time: " << (st1 - st) / CLOCKS_PER_SEC << endl; 

    vector<vector<pair<int,int> > > U = SeamC::get_U();
    // waitKey(0);
    vector<vector<pair<int,int> > > mesh = Mesh_Placement(SeamC::get_seam_carving(),img,U,0);
    double st2 = clock();
    cout << fixed << setprecision(3)<< "mesh cost time: " << (st2 - st1) / CLOCKS_PER_SEC << endl; 

    energy::line_img(img);
    energy::shape_energy(img,mesh);

    waitKey(0);
    return 0;
}