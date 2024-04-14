

#include "mesh_placement.h"

vector<vector<pair<int,int> > > Mesh_Placement(Mat img,Mat image,vector<vector<pair<int,int> > > U,bool is_show)
{
    int n = img.rows,m = img.cols;
    // cout << n << ' ' << m << endl;
    // cout << image.rows << ' ' << image.cols << endl;
    vector<vector<pair<int,int> > > pos;
    vector<pair<int,int> > p;
    for(int i = 0;i <= 20; ++ i)
    {
        p.clear();
        int x = (n - 2) / 20.0 * i;
        if(x > n - 2) x = n - 2;
        for(int j = 0;j <= 20; ++ j)
        {
            int y = (m - 2) / 20.0 * j;
            if(y > m - 2) x = m - 2;
            p.push_back(make_pair(x,y));
        }
        pos.push_back(p);
    }
    Mat imgshow;
    img.copyTo(imgshow);
    if(is_show)
    {
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].first,pos[i][j].second);
        for(int i = 0;i <= 20; ++ i)
        {
            for(int j = 0;j <= 20;++ j)
            {
                // cout << pos[i][j].first << ' '<<pos[i][j].second << ' ' << i << ' ' << j << endl;
                if(i != 0) line(imgshow,Point(pos[i][j].first,pos[i][j].second),Point(pos[i-1][j].first,pos[i-1][j].second),Scalar(0,255,0),2);
                if(j != 0) line(imgshow,Point(pos[i][j].first,pos[i][j].second),Point(pos[i][j-1].first,pos[i][j-1].second),Scalar(0,255,0),2);
            }
        }
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].first,pos[i][j].second);    
        namedWindow("mesh_placement",WINDOW_FREERATIO);
        imshow("mesh_placement",imgshow);

        waitKey(0);
    }

    image.copyTo(imgshow);
    for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
        pos[i][j].first += U[pos[i][j].first][pos[i][j].second].first,
        pos[i][j].second += U[pos[i][j].first][pos[i][j].second].second;
    if(is_show)
    {
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].first,pos[i][j].second);
            
        for(int i = 0;i <= 20; ++ i)
        {
            for(int j = 0;j <= 20;++ j)
            {
                cout << pos[i][j].second << ' '<< pos[i][j].first << ' ' << i << ' ' << j << endl;
                cout << U[pos[i][j].second][pos[i][j].first].first << ' '<< U[pos[i][j].second][pos[i][j].first].second<< ' ' << i << ' ' << j << endl;
                if(i != 0) line(imgshow,Point(pos[i][j].first,pos[i][j].second),Point(pos[i-1][j].first,pos[i-1][j].second),Scalar(0,255,0),1);
                if(j != 0) line(imgshow,Point(pos[i][j].first,pos[i][j].second),Point(pos[i][j-1].first,pos[i][j-1].second),Scalar(0,255,0),1);
                cout << "Draw" << endl;
            }
        }
        for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
            swap(pos[i][j].first,pos[i][j].second);
        namedWindow("mesh_placement",WINDOW_FREERATIO);
        imshow("mesh_placement",imgshow);
    }
    return pos;
}