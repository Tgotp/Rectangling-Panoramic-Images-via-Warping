

#include "mesh_placement.h"

vector<vector<Point > > Mesh_Placement(Mat img,Mat image,vector<vector<Point > > U,bool is_show)
{
    int n = img.rows,m = img.cols;
    // cout << n << ' ' << m << endl;
    // cout << image.rows << ' ' << image.cols << endl;
    vector<vector<Point > > pos;
    vector<Point > p;
    for(int i = 0;i <= 20; ++ i)
    {
        p.clear();
        int x = (n - 2) / 20.0 * i;
        if(x > n - 2) x = n - 2;
        for(int j = 0;j <= 20; ++ j)
        {
            int y = (m - 2) / 20.0 * j;
            if(y > m - 2) y = m - 2;
            p.push_back(Point(x,y));
        }
        pos.push_back(p);
    }
    Mat imgshow;
    img.copyTo(imgshow);
    if(is_show)
    {
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].x,pos[i][j].y);
        for(int i = 0;i <= 20; ++ i)
        {
            for(int j = 0;j <= 20;++ j)
            {
                cout << pos[i][j].x << ' '<<pos[i][j].y << ' ' << i << ' ' << j << endl;
                if(i != 0) line(imgshow,Point(pos[i][j].x,pos[i][j].y),Point(pos[i-1][j].x,pos[i-1][j].y),Scalar(0,255,0),2);
                if(j != 0) line(imgshow,Point(pos[i][j].x,pos[i][j].y),Point(pos[i][j-1].x,pos[i][j-1].y),Scalar(0,255,0),2);
            }
        }
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].x,pos[i][j].y);    
        namedWindow("mesh_placement",WINDOW_FREERATIO);
        imshow("mesh_placement",imgshow);

        waitKey(0);
    }

    image.copyTo(imgshow);
    for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
        pos[i][j].x += U[pos[i][j].x][pos[i][j].y].x,
        pos[i][j].y += U[pos[i][j].x][pos[i][j].y].y;
    if(is_show)
    {
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].x,pos[i][j].y);
            
        for(int i = 0;i <= 20; ++ i)
        {
            for(int j = 0;j <= 20;++ j)
            {
                cout << pos[i][j].y << ' '<< pos[i][j].x << ' ' << i << ' ' << j << endl;
                cout << U[pos[i][j].y][pos[i][j].x].x << ' '<< U[pos[i][j].y][pos[i][j].x].y<< ' ' << i << ' ' << j << endl;
                if(i != 0) line(imgshow,Point(pos[i][j].x,pos[i][j].y),Point(pos[i-1][j].x,pos[i-1][j].y),Scalar(0,0,255),1);
                if(j != 0) line(imgshow,Point(pos[i][j].x,pos[i][j].y),Point(pos[i][j-1].x,pos[i][j-1].y),Scalar(0,0,255),1);
                cout << "Draw" << endl;
            }
        }
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].x,pos[i][j].y);
        namedWindow("mesh_placement",WINDOW_FREERATIO);
        imshow("mesh_placement",imgshow);
        waitKey(0);
    }
    return pos;
}