

#include "mesh_placement.h"

vector<vector<Point > > Mesh_Placement(Mat img,Mat image,vector<vector<Point > > U,vector<vector<Point > >&P,bool is_show)
{
    int n = img.rows,m = img.cols;
    // cout << n << ' ' << m << endl;
    // cout << image.rows << ' ' << image.cols << endl;
    vector<vector<Point > > pos;
    vector<Point > p;
    for(int i = 0;i <= 20; ++ i)
    {
        p.clear();
        int x = (n - 1) / 20.0 * i;
        if(x > n - 1) x = n - 1;
        for(int j = 0;j <= 20; ++ j)
        {
            int y = (m - 1) / 20.0 * j;
            if(y > m - 1) y = m - 1;
            p.push_back(Point(x,y));
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

VectorXd get_mesh_point(vector<vector<Point> >pos)
{
    VectorXd vec = VectorXd(21*21*2);
    for(int i = 0;i <= 20;++ i)
        for(int j = 0;j <= 20;++ j)
            vec((i*21+j)*2) = pos[i][j].x,
            vec((i*21+j)*2+1) = pos[i][j].y;
    return vec;
}

vector<vector<Point > > vec_to_mesh(VectorXd V)
{
    vector<vector<Point > > mesh;
    vector<Point > m;
    for(int i = 0;i <= 20; ++ i)
    {
        m.clear();
        for(int j = 0;j <= 20; ++ j)
        {
            cout << V((i*21+j)*2)<< ' ' << V((i*21+j)*2+1) << endl;
            m.push_back(Point(V((i*21+j)*2),V((i*21+j)*2+1)));
        }
        mesh.push_back(m);
    }
    return mesh;
}

SparseMatrix<double> merge_matrix(SparseMatrix<double> a,SparseMatrix<double> b)
{
    SparseMatrix<double> mat1(a.rows()+b.rows(),a.cols());
    for (int k=0; k<a.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(a,k); it; ++it)
            mat1.insert(it.row(), it.col()) = it.value();

    // 插入矩阵b的元素
    for (int k=0; k<b.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(b,k); it; ++it)
            mat1.insert(it.row() + a.rows(), it.col()) = it.value();
	mat1.makeCompressed();
    return mat1;
}



VectorXd merge_matrix(VectorXd a,VectorXd b)
{
    VectorXd mat1(a.rows()+b.rows());
    for(int i = 0;i < a.rows();++ i)
        mat1(i) = a(i);

    for(int i = 0;i < b.rows();++ i)
        mat1(i+a.rows()) = b(i);
    return mat1;
}