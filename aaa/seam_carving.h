#pragma once
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp>
#include <opencv2\opencv.hpp>
#include <iostream>
#include "tool.h"

using namespace cv;
using namespace std;

namespace SeamC
{
    void init(Mat&img,bool is_show=1);
    
    void seam_carving();
    void get_gray();
    void get_mask();
    void fill_hole();
    void seam_insert();
    void get_rec();
    void get_energy(int sx,int sy,int tx,int ty,bool dir);
    void get_cost(int sx,int sy,bool dir);
    void seam_insert(int sx,int sy,int tx,int ty,int dir);
    Mat get_seam_carving();
    Mat get_rec_img();
    vector<vector<point > > get_U();
};