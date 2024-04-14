#include "energy.h"

double* energy::line_img(Mat img,bool is_show)
{
    int n = img.rows,m = img.cols;

    int *n_line = new int;
    double *input_img = solve_img(img);
    double *out = lsd(n_line,input_img,m,n);
    if(is_show)
    {
        Mat line_img = Mat::zeros(n,m,CV_8UC1);
        line_img = ~line_img;
        int num = * n_line;
        for(int i = 0;i < num;++ i)
        {
            // cout << out[i*7]<< ' ' << out[i*7+1] << ' ' << out[i*7+2] << ' ' << out[i*7+3] << endl;
            line(line_img,Point(out[i*7],out[i*7+1]),Point(out[i*7+2],out[i*7+3]),0,1);
        }
        namedWindow("line",WINDOW_FREERATIO);
        imshow("line",line_img);
        waitKey(0);
    }
    return out;
}

double energy::get_gray(Mat&img,int i,int j) { return img.at<Vec3b>(i,j)[0] * 0.299 + img.at<Vec3b>(i,j)[1] * 0.587 + img.at<Vec3b>(i,j)[2] * 0.114; }

double* energy::solve_img(Mat img)
{
    int n = img.rows,m = img.cols;
    double*ans = new double[n * m];
    for(int j = 0;j < m; ++ j)
        for(int i = 0;i < n; ++ i)
        {
            // cout << n << ' ' << m <<' '<< i << ' ' << j <<' ' <<i +j*n <<endl;
            ans[i*m+j] = get_gray(img,i,j);
        }
    return ans;
}

MatrixXf energy::shape_energy(Mat img,vector<vector<pair<int,int> > > mesh)
{
    int n = img.rows,m = img.cols;
    MatrixXd shape_energy(400*8,20*20*8);
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
        {
            pair<int,int> p0 = mesh[i][j];
            pair<int,int> p1 = mesh[i][j+1];
            pair<int,int> p2 = mesh[i+1][j];
            pair<int,int> p3 = mesh[i+1][j+1];
            MatrixXd Aq(8,4);
            Aq << p0.second, - p0.first , 1 , 0,
                  p0.first,   p0.second , 0 , 1,
                  p1.second, - p1.first , 1 , 0,
                  p1.first,   p1.second , 0 , 1,
                  p2.second, - p2.first , 1 , 0,
                  p2.first,   p2.second , 0 , 1,
                  p3.second, - p3.first , 1 , 0,
                  p3.first,   p3.second , 0 , 1;

            MatrixXd es = Aq * (Aq.transpose() * Aq).inverse() * Aq.transpose() - MatrixXd::Indentity(8,8);
            
            int lx = i * 8,ly = j * 8;
            for(int i = 0;i < 8; ++ i)
                for(int j = 0;j < 8; ++ j)
                    shape_energy(i,j) = es(i,j);
        }
    cout << shape_energy << endl;
    return shape_energy;
}