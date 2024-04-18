

#include "mesh_placement.h"

vector<vector<point > > Mesh_Placement(Mat img,Mat image,vector<vector<point > > U,vector<vector<point > >&P,bool is_show)
{
    int n = img.rows,m = img.cols;
    // cout << n << ' ' << m << endl;
    // cout << image.rows << ' ' << image.cols << endl;
    vector<vector<point > > pos;
    vector<point > p;
    for(int i = 0;i <= 20; ++ i)
    {
        p.clear();
        int x = (n - 3) / 20.0 * i;
        if(x > n - 3) x = n - 3;
        for(int j = 0;j <= 20; ++ j)
        {
            int y = (m - 3) / 20.0 * j;
            if(y > m - 3) y = m - 3;
            p.push_back(point(x,y));
        }
        pos.push_back(p);
    }
    P = pos;
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


void img_mesh(Mat&img,vector<vector<point > > pos)
{
    int n = img.rows,m = img.cols;
    for(int i = 0;i <= 20; ++ i) for(int j = 0;j <= 20;++ j)
        swap(pos[i][j].x,pos[i][j].y);
    for(int i = 0;i <= 20; ++ i)
    {
        for(int j = 0;j <= 20;++ j)
        {
            // cout << pos[i][j].x << ' '<<pos[i][j].y << ' ' << i << ' ' << j << endl;
            if(i != 0) line(img,Point(pos[i][j].x,pos[i][j].y),Point(pos[i-1][j].x,pos[i-1][j].y),Scalar(0,255,0),2);
            if(j != 0) line(img,Point(pos[i][j].x,pos[i][j].y),Point(pos[i][j-1].x,pos[i][j-1].y),Scalar(0,255,0),2);
        }
    } 
    namedWindow("mesh_placement",WINDOW_FREERATIO);
    imshow("mesh_placement",img);

    waitKey(0);

}

VectorXd get_mesh_point(vector<vector<point> >pos)
{
    VectorXd vec = VectorXd(21*21*2);
    for(int i = 0;i <= 20;++ i)
        for(int j = 0;j <= 20;++ j)
            vec((i*21+j)*2) = pos[i][j].x,
            vec((i*21+j)*2+1) = pos[i][j].y;
    return vec;
}

vector<vector<point > > vec_to_mesh(VectorXd V)
{
    vector<vector<point > > mesh;
    vector<point > m;
    for(int i = 0;i <= 20; ++ i)
    {
        m.clear();
        for(int j = 0;j <= 20; ++ j)
        {
            // if(V((i*21+j)*2) < -0.0000001 || V((i*21+j)*2+1)<-0.0000001) 
            // {
            //     printf("error qwq\n");
            //     cout << i << ' ' << j << endl;
            //     cout << V((i*21+j)*2) <<' ' <<V((i*21+j)*2+1) << endl;
            // }
            // if(V((i*21+j)*2) > 289 || V((i*21+j)*2+1) > 633) 
            //     printf("error qaq\n");
            // cout << V((i*21+j)*2)<< ' ' << V((i*21+j)*2+1) << endl;
            m.push_back(point(V((i*21+j)*2),V((i*21+j)*2+1)));
        }
        mesh.push_back(m);
    }
    return mesh;
}
