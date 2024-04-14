
    
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp>
#include <opencv2\opencv.hpp>
#include <iostream>
#include "seam_carving.h"
using namespace cv;
using namespace std;

int A[6] = {1,1,1,-1,-1,-1};
int B[6] = {-1,0,1,-1,0,1};
int W[6] = {1,2,1,-1,-2,-1};
Mat rimg,img,cimg,mask,maskx,masky;
int n,m; bool cshow;
vector<vector<double > > G,E,C;
vector<vector<Point > > U; // 位移场

Mat SeamC::get_seam_carving() { return img; }
Mat SeamC::get_rec_img() { return rimg; }
vector<vector<Point > > SeamC::get_U() { return U; }
void SeamC::init(Mat&image,bool is_show)
{
    img = image;
    n = img.rows; m = img.cols;
    // cout << n << 'x' << m << endl;
    get_gray();
    get_mask();
    get_rec();

    cshow = is_show;
    if(cshow) img.copyTo(cimg);
    seam_carving();

    if(cshow)
    {
        imshow("carvd image", img);
    }
}

void SeamC::get_rec()
{
    int sx,sy,tx,ty;bool flag;
    for(int i = 0,flag = 1;i < n && flag; ++ i)
        for(int j = 0;j < m;++ j)
            if(mask.at<uchar>(i,j))
            {
                sx = i;flag = 0;
                break;
            }
    for(int i = n - 1,flag = 1;i >= 0 && flag; -- i)
        for(int j = 0;j < m;++ j)
            if(mask.at<uchar>(i,j))
            {
                tx = i;flag = 0;
                break;
            }
    
    for(int j = m - 1,flag = 1;j >= 0 && flag; -- j)
        for(int i = 0;i < n;++ i)
            if(mask.at<uchar>(i,j))
            {
                ty = j;flag = 0;
                break;
            }
    for(int j = 0,flag = 1;j < m && flag; ++ j)
        for(int i = 0;i < n;++ i)
            if(mask.at<uchar>(i,j))
            {
                sy = j;flag = 0;
                break;
            }

    cout << n << 'x' << m << endl;
    cout << sx << ' ' << sy << ' ' << tx << ' ' << ty << endl;

    Mat imgROI(img,Rect(sy,sx,ty - sy + 1,tx - sx + 1));
    imgROI.copyTo(img);
    Mat maskROI(mask,Rect(sy,sx,ty - sy + 1,tx - sx + 1));
    maskROI.copyTo(mask);
    // imshow("origin image",img);
    // waitKey(0);
    img.copyTo(rimg);

    
    n = img.rows; m = img.cols;
    mask.copyTo(maskx);
    mask.copyTo(masky);
    get_gray();

    vector<Point > u;
    for(int i = 0;i < m; ++ i)
        u.push_back(Point(0,0));
    for(int i = 0;i < n; ++ i)
        U.push_back(u);
    
    // Mat mask_img;
    // img.copyTo(mask_img,mask);
    // namedWindow("erode_img",WINDOW_FREERATIO);
    // imshow("erode_img",mask);
}

void SeamC::get_mask()
{
    mask = Mat :: zeros(n,m,CV_8UC1);
    for(int i = 0;i < n;++ i)
        for(int j = 0;j < m;++ j)
            mask.at<uchar>(i,j) = abs(G[i][j] - 255) < 3 ? 255 : 0;
    fill_hole();

    Mat erode_out,dilate_out;//腐蚀,膨胀
    Mat element = getStructuringElement(MORPH_ELLIPSE,Size(5, 5));
	cv::dilate(mask, dilate_out, element);
	cv::dilate(dilate_out, dilate_out, element);
	cv::dilate(dilate_out, dilate_out, element);
    erode(dilate_out,erode_out,getStructuringElement(MORPH_ELLIPSE,Size(5, 5)));    
    // namedWindow("erode_out");
    // imshow("erode_out",erode_out);
    
    mask = ~erode_out;
    
}

void SeamC::fill_hole()
{
    int sx[n],sy[m],tx[n] = {},ty[m] = {};
    memset(sx,0x3f,sizeof sx);
    memset(sy,0x3f,sizeof sy);
    for(int i = 0;i < n; ++ i)
        for(int j = 0;j < m;++ j) if(mask.at<uchar>(i,j) == 0)
        {
            sx[i] = min(sx[i],j);
            sy[j] = min(sy[j],i);

            tx[i] = max(tx[i],j);
            ty[j] = max(ty[j],i);
        }
    
    for(int i = 0;i < n; ++ i)
        for(int j = 0;j < m;++ j) if(mask.at<uchar>(i,j) == 255)
            if(j < tx[i] && j > sx[i] && i < ty[j] && i > sy[j])
                mask.at<uchar>(i,j) = 0;
}
void SeamC::seam_carving()
{
    int dir = 0,loop = 0;
    do
    {
        // dir = 0表示不存在，1表示第一行，3表示最后一行，2第一列，4最后一列
        int x=0,y=0,lst = 0;
        for(int i = 0;i < m;++ i)
            if(!mask.at<uchar>(0,i))
            {
                if(i - lst > y - x)
                    y = i,x = lst,dir = 1;
            }else lst = i + 1;
        lst = 0;
        for(int i = 0;i < m;++ i)
            if(!mask.at<uchar>(n-1,i))
            {
                if(i - lst > y - x)
                    y = i,x = lst,dir = 3;
            }else lst = i + 1;
        lst = 0;
        for(int i = 0;i < n;++ i)
            if(!mask.at<uchar>(i,0))
            {
                if(i - lst > y - x)
                    y = i,x = lst,dir = 2;
            }else lst = i + 1;
        lst = 0;
        for(int i = 0;i < n;++ i)
            if(!mask.at<uchar>(i,m - 1))
            {
                if(i - lst > y - x)
                    y = i,x = lst,dir = 4;
            }else lst = i + 1;
        
        // cout << "direction: "<< dir << " len: " << y - x << " x: "<<x << " y: " << y << endl;

        if(y - x == 0) break;

        if(dir & 1)
        {
            get_energy(0,x,n,y,dir&1);
            // cout << E.size() << endl;
            get_cost(0,x,dir&1);
            // cout << C.size() << endl;
            seam_insert(0,x,n,y,dir);
        }
        else
        {
            get_energy(x,0,y,m,dir&1);
            get_cost(x,0,dir&1);
            seam_insert(x,0,y,m,dir);
        }
        if(cshow)
        {
            namedWindow("mask",WINDOW_FREERATIO);
            imshow("mask", mask);
            namedWindow("carved",WINDOW_FREERATIO);
            imshow("carved", img);
            namedWindow("carving",WINDOW_FREERATIO);
            imshow("carving", cimg);
            // namedWindow("maskx",WINDOW_FREERATIO);
            // imshow("maskx", maskx);
            // namedWindow("masky",WINDOW_FREERATIO);
            // imshow("masky", masky);
            char key = waitKey(5);
            if(key == 'q') exit(0);
        }
        if(y - x <= 8)break;
    }while(true);
}

void SeamC::get_gray()
{
    G.clear();
    vector<double > g;
    for(int i = 0 ;i < n;++ i)
    {
        g.clear();
        for(int j = 0;j < m;++ j)
            g.push_back(img.at<Vec3b>(i,j)[0] * 0.299 + img.at<Vec3b>(i,j)[1] * 0.587 + img.at<Vec3b>(i,j)[2] * 0.114);
        G.push_back(g);
    }
}

void SeamC::get_energy(int sx,int sy,int tx,int ty,bool dir)
{
    E.clear();
    vector<double > e;
    for(int i = sx;i < tx;++ i)
    {
        e.clear();
        for(int j = sy;j < ty;++ j)
            if(dir ? masky.at<uchar>(i,j): maskx.at<uchar>(i,j))
            {
                double ex = 0,ey = 0;
                for(int k = 0;k < 6;++ k) 
                {
                    int Xx = i + A[k],Xy = j + B[k];
                    int Yx = i + B[k],Yy = j + A[k];
                    // cout <<"solve_energy" << Xx <<" " << Xy <<" "<< Yx <<" "<< Yy <<" "<< endl;
                    if(Xx >= 0 && Xy >= 0 && Xx < n && Xy < m)
                        ex += W[k] * G[Xx][Xy];
                    
                    if(Yx >= 0 && Yy >= 0 && Yx < n && Yy < m)
                        ey += W[k] * G[Yx][Yy];
                }
                e.push_back(abs(ex) * 0.5 + abs(ey) * 0.5);
            }else e.push_back(1e5);
        E.push_back(e);
    }
    // cout <<"solve_energy" << endl;
}

double get(int i,int j,int dir)
{
    if(dir)
    {
        if(j == 0) return 0;
        if(i < n - 1 && i > 0) return min(C[i][j - 1],min(C[i + 1][j - 1],C[i - 1][j - 1]));
        if(i < n - 1 && i >= 0) return min(C[i][j - 1],C[i + 1][j - 1]);
        if(i < n && i > 0) return min(C[i][j - 1],C[i - 1][j - 1]);
        return C[i][j];
    }
    else
    {
        if(i == 0) return 0;
        if(j < m - 1 && j > 0) return min(C[i - 1][j],min(C[i - 1][j - 1],C[i - 1][j + 1]));
        if(j < m - 1 && j >= 0) return min(C[i - 1][j],C[i - 1][j + 1]);
        if(j < m && j > 0) return min(C[i - 1][j],C[i - 1][j - 1]);
        return C[i][j];
    }
}

void SeamC::get_cost(int sx,int sy,bool dir)
{
    // cout << E.size() << endl;
    // cout << E[0].size() << endl;
    
    C.clear();
    vector<double> c;
    for(int i = 0;i < E[0].size();++ i)
        c.push_back(0);
    for(int i = 0;i < E.size();++ i)
        C.push_back(c);
    // cout << "get cost sx,sy: "<< sx << " " << sy << endl;
    if(dir)
    {
        for(int j = 0; j < E[0].size(); ++ j)
            for(int i = 0; i < E.size(); ++ i)   
                // if(mask.at<uchar>(i + sx,j + sy))
                    C[i][j] = E[i][j] + get(i,j,dir);
    }
    else
    {
        for(int i = 0; i < E.size(); ++ i)
            for(int j = 0; j < E[0].size(); ++ j)   
                // if(mask.at<uchar>(i + sx,j + sy))
                    C[i][j] = E[i][j] + get(i,j,dir);
    }
                
    // cout << "get cost end" << endl;
}

Vec3b get_color(int x,int y,int dir)
{
    switch(dir)
    {
        case 1: return x == n - 1 ? img.at<Vec3b>(x,y) : img.at<Vec3b>(x,y) * 0.5 + img.at<Vec3b>(x+1,y) * 0.5;
        case 2: return y == m - 1 ? img.at<Vec3b>(x,y) : img.at<Vec3b>(x,y) * 0.5 + img.at<Vec3b>(x,y+1) * 0.5;
        case 3: return x == 0 ? img.at<Vec3b>(x,y) : img.at<Vec3b>(x,y) * 0.5 + img.at<Vec3b>(x-1,y) * 0.5;
        case 4: return y == 0 ? img.at<Vec3b>(x,y) : img.at<Vec3b>(x,y) * 0.5 + img.at<Vec3b>(x,y-1) * 0.5;
    }
}   

void SeamC::seam_insert(int sx,int sy,int tx,int ty,int dir)
{
    // cout << "seam insert dir :" << dir << endl;
    int n = tx - sx,m = ty - sy;
    // cout << n << ' ' << m << endl;
    // cout << sy << ' ' << ty << endl;

    vector<int> pos;
    if(dir&1)
    {

        double mx = 1e18,x = 0,w = 0;
        
        for(int i = 0; i < C.size(); ++ i)   
            if(mx > C[i][m - 1])
            {
                mx = C[i][m - 1];
                x = i; w = mx - E[i][m - 1];
            }
        pos.push_back(x);
        
        // cout <<fixed << setprecision(3) << "value: " << mx << endl;

        for(int j = m - 2; j >= 0;j --)
            for(int i = 0;i < n; ++ i)
                if(abs(w - C[i][j]) < 0.0001 && abs(i - x) <= 1 )
                {
                    x = i; w -= E[i][j];
                    pos.push_back(i);
                    break;
                }
        reverse(pos.begin(),pos.end());
        // cout << "size:" << pos.size() << ' ' << w<< endl;
        if(dir == 1)
            for(int j = 0 ;j < m;++ j)
                for(int i = 0 ;i < pos[j];++ i)
                {
                    U[i][j + sy] = U[i + 1][j + sy];
                    U[i][j + sy].x += 1;
                    img.at<Vec3b>(i,j + sy) = img.at<Vec3b>(i + 1,j + sy);
                    mask.at<uchar>(i,j + sy) = mask.at<uchar>(i + 1,j + sy);
                }
        else
            for(int j = 0 ;j < m;++ j)
                for(int i = n - 1 ;i > pos[j];-- i)
                {
                    U[i][j + sy] = U[i - 1][j + sy];
                    U[i][j + sy].x -= 1;
                    img.at<Vec3b>(i,j + sy) = img.at<Vec3b>(i - 1,j + sy);
                    mask.at<uchar>(i,j + sy) = mask.at<uchar>(i - 1,j + sy);
                }
        
        // cout << sy << endl;
        for(int j = 0 ;j < m;++ j)
        {
            // cout << pos[j] << ' ' << j << ' ' << j + sy << ' ' << img.cols << ' ' << masky.cols<< endl;
            masky.at<uchar>(pos[j],j + sy) = 0;
            img.at<Vec3b>(pos[j],j + sy) = get_color(pos[j]+sx,j+sy,dir);
        }

        if(cshow)
        {
            for(int j = 0 ;j < m;++ j)
            {
                if(dir == 1)
                for(int i = 0 ;i < pos[j];++ i)
                    cimg.at<Vec3b>(i,j + sy) = cimg.at<Vec3b>(i + 1,j + sy);
                else
                for(int i = n - 1 ;i > pos[j];-- i)
                    cimg.at<Vec3b>(i,j + sy) = cimg.at<Vec3b>(i - 1,j + sy);
            }
            for(int j = 0 ;j < m;++ j)
                cimg.at<Vec3b>(pos[j],j + sy) = Vec3b(0,0,255);
        }
    }else
    {
        double mx = 1e18,x = 0,w = 0;
        
        for(int i = 0; i < C[0].size(); ++ i)   
            if(mx > C[n - 1][i])
            {
                mx = C[n - 1][i];
                x = i; w = mx - E[n - 1][i];
            }
        pos.push_back(x);
        for(int i = n - 2;i >= 0; -- i)
            for(int j = 0; j < m ;j ++)
                if(abs(w - C[i][j]) < 0.1 && abs(j-x) <= 1)
                {
                    x = j; w -= E[i][j];
                    // cout << w << endl;
                    pos.push_back(j);
                    break;
                }
        reverse(pos.begin(),pos.end());
        // cout << pos.size() << endl;
        if(dir == 2)
            for(int i = 0 ;i < n;++ i)
                for(int j = 0;j < pos[i];++ j)
                {
                    U[i + sx][j] = U[i + sx][j + 1];
                    U[i + sx][j].y += 1;
                    img.at<Vec3b>(i + sx,j) = img.at<Vec3b>(i + sx,j + 1);
                    mask.at<uchar>(i + sx,j) = mask.at<uchar>(i + sx,j + 1);
                }
        else
            for(int i = 0;i < n;++ i)
                for(int j = m - 1;j > pos[i];-- j)
                {
                    U[i + sx][j] = U[i + sx][j - 1];
                    U[i + sx][j].y -= 1;
                    img.at<Vec3b>(i + sx,j) = img.at<Vec3b>(i + sx,j - 1);
                    mask.at<uchar>(i + sx,j) = mask.at<uchar>(i + sx,j - 1);
                }
            
        for(int i = 0 ;i < n;++ i)
        {
            maskx.at<uchar>(i + sx,pos[i]) = 0;
            img.at<Vec3b>(i + sx,pos[i]) = get_color(i+sx,pos[i],dir);
        }
        if(cshow)
        {
            for(int i = 0 ;i < n;++ i)
            {  
                if(dir == 2)
                    for(int j = 0;j < pos[i];++ j)
                        cimg.at<Vec3b>(i + sx,j) = cimg.at<Vec3b>(i + sx,j + 1);
                else
                    for(int j = m - 1;j > pos[i];-- j)
                        cimg.at<Vec3b>(i + sx,j) = cimg.at<Vec3b>(i + sx,j - 1);
            }
            for(int i = 0 ;i < n;++ i)
                cimg.at<Vec3b>(i + sx,pos[i]) = Vec3b(0,0,255);
        }
    }
}