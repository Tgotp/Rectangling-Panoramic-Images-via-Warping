#include "energy.h"

SparseMatrix<double> energy::shape_energy(vector<vector<Point > > mesh)
{
    // double  = 0;
    SparseMatrix<double> shape_energy(20*20*8,21*21*2);
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
        {
            Point p0 = mesh[i][j];
            Point p1 = mesh[i][j+1];
            Point p2 = mesh[i+1][j];
            Point p3 = mesh[i+1][j+1];
            MatrixXd Aq(8,4);
            Aq << p0.y, - p0.x , 1 , 0,
                  p0.x,   p0.y , 0 , 1,
                  p1.y, - p1.x , 1 , 0,
                  p1.x,   p1.y , 0 , 1,
                  p2.y, - p2.x , 1 , 0,
                  p2.x,   p2.y , 0 , 1,
                  p3.y, - p3.x , 1 , 0,
                  p3.x,   p3.y , 0 , 1;
            // MatrixXd Vq(8,1);
            // Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
            //     V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;

            Aq = (Aq * (Aq.transpose() * Aq).inverse() * Aq.transpose() - MatrixXd::Identity(8,8));
            Aq = Aq.transpose()*Aq;
            
            // cout << i << ' ' << j << ' ' << (i*20+j)*8+7 << ' '<< 20*20*8<< endl;
            for(int x = 0;x < 8;++ x)
            {
                shape_energy.insert((i*20+j)*8+x,(i*21+j)*2) = Aq(x,0);
                shape_energy.insert((i*20+j)*8+x,(i*21+j)*2+1) = Aq(x,1);
                shape_energy.insert((i*20+j)*8+x,(i*21+j+1)*2) = Aq(x,2);
                shape_energy.insert((i*20+j)*8+x,(i*21+j+1)*2+1) = Aq(x,3);
                shape_energy.insert((i*20+j)*8+x,((i+1)*21+j)*2) = Aq(x,4);
                shape_energy.insert((i*20+j)*8+x,((i+1)*21+j)*2+1) = Aq(x,5);
                shape_energy.insert((i*20+j)*8+x,((i+1)*21+j+1)*2) = Aq(x,6);
                shape_energy.insert((i*20+j)*8+x,((i+1)*21+j+1)*2+1) = Aq(x,7);
            }
            
            // for(int x = 0;x < 8; ++ x)
            //     shape_energy += es(x,0);
        }
    // cout << shape_energy << endl;
    cout << "solved shape energy" << endl;
	shape_energy.makeCompressed();
    return shape_energy;
}

double energy::line_energy(vector<vector<vector<Line> > > mesh,vector<vector<Point > > V)
{
    return 0;
}

pair<SparseMatrix<double>,VectorXd > energy::bound_energy(vector<vector<Point > > mesh,double inf,int n,int m)
{
    // cout << "solve bound energy" << endl;
    pair<SparseMatrix<double>,VectorXd > bound_energy; 
    SparseMatrix<double> a(84,21*21*2);
    int cnt = 0;
    VectorXd b(84);
    for(int i = 0;i <= 20;++ i)
    {
        a.insert(cnt,i*2) = inf;
        b(cnt) = 0; // top
        ++ cnt;
    }
    for(int i = 0;i <= 20;++ i)
    {
        a.insert(cnt,i*21*2+1) = inf,
        b(cnt) = 0; // left
        ++ cnt;
    }
    for(int i = 0;i <= 20;++ i)
    {
        a.insert(cnt,(20*21+i)*2) = inf,  // bottom
        b(cnt) = inf * (n-1); 
        ++ cnt;
    }
    for(int i = 0;i <= 20;++ i)
    {
        a.insert(cnt,(i*21+20)*2+1) = inf,  // right
        b(cnt) = inf * (m-1); 
        ++ cnt;
    }
    bound_energy = make_pair(a,b);
    // cout << cnt << endl;
    cout << "solved bound energy" << endl;
    return bound_energy;
}

double inverse_bilinear_interpolation(vector<Point> P[2],Line X)
{
    MatrixXd A = MatrixXd::Zero(8,8);
    VectorXd B(8,1);
    for(int i = 0;i < 4; ++ i)
    {
        A(i*2,0) = P[0][i].x;
        A(i*2,1) = P[0][i].y;
        A(i*2,2) = P[0][i].x * P[0][i].y;
        A(i*2,3) = 1;
        
        A(i*2+1,4) = P[0][i].x;
        A(i*2+1,5) = P[0][i].y;
        A(i*2+1,6) = P[0][i].x * P[0][i].y;
        A(i*2+1,7) = 1;

        B(i*2) = P[1][i].x;
        B(i*2+1) = P[1][i].y;
    }
    // cout << "----A------"<<endl;
    // cout << A << endl;
    // cout << "----B------"<<endl;
    // cout << B << endl;
    MatrixXd F = A.colPivHouseholderQr().solve(B);
    MatrixXd x(2,8),y(2,8);
    x << X.x.x,X.x.y,X.x.x*X.x.y,1,0,0,0,0,
        0,0,0,0,X.x.x,X.x.y,X.x.x*X.x.y,1;
    y << X.y.x,X.y.y,X.y.x*X.y.y,1,0,0,0,0,
        0,0,0,0,X.y.x,X.y.y,X.y.x*X.y.y,1;
    MatrixXd a = x * F;
    MatrixXd b = y * F;
    double q = b(0) - a(0);
    double e = b(1) - a(1);
    // cout << X.x << ' ' << X.y << endl;
    // cout <<"----a----"<< endl << a << endl;
    // cout <<"----b----"<< endl << b << endl;
    // cout << acos(q / sqrt(q*q + e*e)) << endl;
    return acos(q / sqrt(q*q + e*e));
}

double* energy::Line_rotate_count(int *num,vector<vector<vector<Line > > > mesh_line,vector<vector<Point> >mesh,vector<vector<Point> >V,int number)
{
    double *bins;
    bins = new double[number];
    vector<Point> p[2];
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
            for(auto k:mesh_line[i][j])
            {
                p[0].clear();p[1].clear();
                p[0].push_back(mesh[i][j]); p[0].push_back(mesh[i][j+1]);
                p[0].push_back(mesh[i+1][j+1]); p[0].push_back(mesh[i+1][j]);
                p[1].push_back(V[i][j]); p[1].push_back(V[i][j+1]);
                p[1].push_back(V[i+1][j+1]); p[1].push_back(V[i+1][j]);
                k.theta = inverse_bilinear_interpolation(p,k);
                k.pos = floor(k.theta/PI * 50);
                k.pos = min(49,k.pos);
                // cout <<"theta and pos : "<< k.theta << ' ' << k.pos << endl;
                bins[k.pos] += k.theta;
                num[k.pos] ++;
            }
    for(int i = 0;i < number; ++ i)
        if(num[i]) bins[i] /= num[i];
    // cout << "Line rotate count solved out" << endl;
    return bins;
}

