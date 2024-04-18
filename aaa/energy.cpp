#include "energy.h"


SparseMatrix<double> energy::shape_energy(vector<vector<point > > mesh)
{
    // double  = 0;
    SparseMatrix<double> shape_energy(20*20*8,21*21*2);
    for(int i = 0;i < 20;++ i)
        for(int j = 0;j < 20;++ j)
        {
            point p0 = mesh[i][j];
            point p1 = mesh[i][j+1];
            point p2 = mesh[i+1][j];
            point p3 = mesh[i+1][j+1];
            MatrixXd Aq(8,4);
            Aq << p0.x, - p0.y , 1 , 0,
                  p0.y,   p0.x , 0 , 1,
                  p1.x, - p1.y , 1 , 0,
                  p1.y,   p1.x , 0 , 1,
                  p2.x, - p2.y , 1 , 0,
                  p2.y,   p2.x , 0 , 1,
                  p3.x, - p3.y , 1 , 0,
                  p3.y,   p3.x , 0 , 1;
            // MatrixXd Vq(8,1);
            // Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
            //     V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;

            Aq = (Aq * (Aq.transpose() * Aq).inverse() * Aq.transpose() - MatrixXd::Identity(8,8));
            // Aq = Aq.transpose()*Aq;
            
            // if(i == 0 && j == 0) output(Aq,"Aq");
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
            // int t = (i*20+j)*8;
            // for(int x = 0;x < 8;++ x)
            //     for(int y = 0;y < 8;++ y)
            //         shape_energy.insert(t+x,y) = Aq(x,y);
        }
    // cout << shape_energy << endl;
    // cout << "solved shape energy" << endl;
	shape_energy.makeCompressed();
    return shape_energy;
}

pair<SparseMatrix<double>,VectorXd > energy::bound_energy(vector<vector<point > > mesh,double inf,int n,int m)
{
    // cout << "solve bound energy" << endl;
    pair<SparseMatrix<double>,VectorXd > bound_energy; 
    SparseMatrix<double> a(80,21*21*2);
    int cnt = 0;
    VectorXd b(80);
    for(int i = 0;i <= 20;++ i)
    {
        a.insert(cnt,i*2) = inf;
        b(cnt) = 0; // top
        ++ cnt;
    }
    for(int i = 1;i < 20;++ i)
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
    for(int i = 1;i <= 19;++ i)
    {
        a.insert(cnt,(i*21+20)*2+1) = inf,  // right
        b(cnt) = inf * (m-1); 
        ++ cnt;
    }
    bound_energy = make_pair(a,b);
    // cout << cnt << endl;
    // cout << "solved bound energy" << endl;
    return bound_energy;
}

MatrixXd inverse_bilinear_interpolation_function(vector<point> P[2],Line X)
{
    MatrixXd A = MatrixXd::Zero(8,8);
    VectorXd B(8,1);
    for(int i = 0;i < 4; ++ i)
    {
        A(i*2,0) = P[0][i].y;
        A(i*2,1) = P[0][i].x;
        A(i*2,2) = P[0][i].x * P[0][i].y;
        A(i*2,3) = 1;
        
        A(i*2+1,4) = P[0][i].y;
        A(i*2+1,5) = P[0][i].x;
        A(i*2+1,6) = P[0][i].x * P[0][i].y;
        A(i*2+1,7) = 1;

        B(i*2) = P[1][i].y;
        B(i*2+1) = P[1][i].x;
    }
    // cout << "----A------"<<endl;
    // cout << A << endl;
    // cout << "----B------"<<endl;
    // cout << B << endl;
    MatrixXd F = A.colPivHouseholderQr().solve(B);
    return F;
}

// MatrixXd zzz,zzzz;
point inverse_bilinear_interpolation(vector<point> P[2],Line X)
{
    // cout <<"invers V :  " << P[0] << ' ' << P[1] << endl;
    MatrixXd F = inverse_bilinear_interpolation_function(P,X);
    MatrixXd x(2,8),y(2,8),z(2,8);
    x << X.x.x,X.x.y,X.x.x*X.x.y,1,0,0,0,0,
        0,0,0,0,X.x.x,X.x.y,X.x.x*X.x.y,1;
    y << X.y.x,X.y.y,X.y.x*X.y.y,1,0,0,0,0,
        0,0,0,0,X.y.x,X.y.y,X.y.x*X.y.y,1;
    MatrixXd a = x * F;
    MatrixXd b = y * F;
    // zzz = a;zzzz = b;
    double q = a(0) - b(0);
    double e = a(1) - b(1);
    // cout << X.x << ' ' << X.y << endl;
    // cout <<"----a----"<< endl << a << endl;
    // cout <<"----b----"<< endl << b << endl;
    // cout << q << ' ' << e << endl;
    // cout << X.x.x - X.y.x << ' ' << X.x.y - X.y.y << endl;
    return point(q,e);
    // if(q == 0 && e == 0) return 4;
    // cout <<"theta :" << q/ sqrt(q*q + e*e) << " " << acos((q)/ sqrt(q*q + e*e))<< endl;
    // cout << acos(q / sqrt(q*q + e*e)) << endl;
    // return acos((q)/ sqrt(q*q + e*e));
}

void energy::Line_rotate_count(double *bins,int *num,vector<vector<vector<Line > > >&mesh_line,vector<vector<point> >mesh,vector<vector<point> >V,int number)
{
    // cout << "solve Line rotate count" << endl;
    for(int i = 0;i < number;++ i) bins[i] = 0;
    vector<point> p[2];
    for(int i = 0;i < 20;++ i)
    {
        for(int j = 0;j < 20;++ j)
        {
            for(auto&k: mesh_line[i][j])
            {
                p[0].clear();p[1].clear();
                p[0].push_back(mesh[i][j]); p[0].push_back(mesh[i][j+1]);
                p[0].push_back(mesh[i+1][j+1]); p[0].push_back(mesh[i+1][j]);
                p[1].push_back(V[i][j]); p[1].push_back(V[i][j+1]);
                p[1].push_back(V[i+1][j+1]); p[1].push_back(V[i+1][j]); 
                point t = inverse_bilinear_interpolation(p,k);
                // cout << t <<' ' << k.x - k.y << endl;
                if(t.x == 0 && t.x == 0) { k.solve = 0; continue; }
                else k.solve = 1;

                double theta = t.x / (sqrt(t.x * t.x + t.y * t.y));
                if(theta > 1) theta = 1;
                theta = acos ( theta );

                bins[k.pos] += theta - k.theta;
                // num[k.pos] ++;

                // Mat img = Mat::zeros(700,1200,CV_8UC3);
                // line(img,Point(k.x.x,k.x.y),Point(k.y.x,k.y.y),Scalar(0,255,0));
                // line(img,Point(mesh[i][j].y,mesh[i][j].x),Point(mesh[i][j+1].y,mesh[i][j+1].x),Scalar(0,255,0));
                // line(img,Point(mesh[i][j].y,mesh[i][j].x),Point(mesh[i+1][j].y,mesh[i+1][j].x),Scalar(0,255,0));
                // line(img,Point(mesh[i+1][j+1].y,mesh[i+1][j+1].x),Point(mesh[i][j+1].y,mesh[i][j+1].x),Scalar(0,255,0));
                // line(img,Point(mesh[i+1][j].y,mesh[i+1][j].x),Point(mesh[i+1][j+1].y,mesh[i+1][j+1].x),Scalar(0,255,0));
                
                // line(img,Point(zzz(0),zzz(1)),Point(zzzz(0),zzzz(1)),Scalar(0,255,0));
                // line(img,Point(V[i][j].y,V[i][j].x),Point(V[i][j+1].y,V[i][j+1].x),Scalar(0,0,255));
                // line(img,Point(V[i][j].y,V[i][j].x),Point(V[i+1][j].y,V[i+1][j].x),Scalar(0,0,255));
                // line(img,Point(V[i+1][j+1].y,V[i+1][j+1].x),Point(V[i][j+1].y,V[i][j+1].x),Scalar(0,0,255));
                // line(img,Point(V[i+1][j].y,V[i+1][j].x),Point(V[i+1][j+1].y,V[i+1][j+1].x),Scalar(0,0,255));
                // cout << "X : \n" << k.x - k.y << "Trans X: \n" << t << endl;
                // namedWindow("trans k: ",WINDOW_FREERATIO);
                // imshow("trans k: ",img);
                // if(waitKey(0) == 'q') exit(0);
            }
        }
    }
    for(int i = 0;i < number; ++ i)
        if(num[i]) bins[i] /= num[i];
    // cout << "-----bins----" << endl;
    // for(int i = 0;i < number; ++ i)
    //     cout << bins[i] << ' '<< endl;
    // cout << "Line rotate count solved out" << endl;
}

MatrixXd get_bilinear_mat(vector<point> p,Line L)
{
    // cout << "solve get bilinear mat" << endl;
    MatrixXd e1(4,8);
    e1 << 
        L.x.x,L.x.y,L.x.x*L.x.y,1,0,0,0,0,
        0,0,0,0,L.x.x,L.x.y,L.x.x*L.x.y,1,
        L.y.x,L.y.y,L.y.x*L.y.y,1,0,0,0,0,
        0,0,0,0,L.y.x,L.y.y,L.y.x*L.y.y,1;
    
    MatrixXd kx = MatrixXd::Zero(8,8);
    for(int t = 0;t < 4; ++ t)
    {
        kx(t*2,0) = p[t].y;
        kx(t*2,1) = p[t].x;
        kx(t*2,2) = p[t].x * p[t].y;
        kx(t*2,3) = 1;

        kx(t*2+1,4) = p[t].y;
        kx(t*2+1,5) = p[t].x;
        kx(t*2+1,6) = p[t].x * p[t].y;
        kx(t*2+1,7) = 1;
    }
    // cout << "line energy Kx: " << kx.rows() << ' ' << kx.cols() << endl;
    // output(kx,"kx");
    // output(kx.inverse(),"inv kx");
    MatrixXd K = e1*kx.inverse();

    // cout << "solved get bilinear mat" << endl;
    return K;
}

bool check(MatrixXd X,int m,int n)
{
    double x = X(1),y = X(0);
    return x < n && y < m && x >= 0 && y >= 0;
}

SparseMatrix<double> energy::line_energy(vector<vector<vector<Line> > > mesh_line,vector<vector<point > > mesh,vector<vector<point > > V,double *bins,int num,int n,int m)
{
    // cout << "line energy solve" << endl;
    // cout << num << endl;
    SparseMatrix<double> line_energy(num * 2,21 * 21 * 2);
    vector<point> p;
    int cnt = 0;
    double cost = 0;
    for(int i = 0;i < 20;++ i)
    {
        for(int j = 0;j < 20;++ j)
        {
            for(auto k : mesh_line[i][j]) 
            {
                if(k.solve == 0) continue; 
                // cout << "line energy i,j: " << i << ' ' << j << endl;
                MatrixXd R(2,2);
                R << cos(bins[k.pos]),-sin(bins[k.pos]),
                        sin(bins[k.pos]),cos(bins[k.pos]);
                MatrixXd e(2,1);
                e << k.x.x - k.y.x,
                        k.x.y - k.y.y;
                    
                MatrixXd C = R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose() - MatrixXd::Identity(2,2);
                
                p.clear();
                p.push_back(mesh[i][j]); p.push_back(mesh[i][j+1]);
                p.push_back(mesh[i+1][j]); p.push_back(mesh[i+1][j+1]);

                MatrixXd K = get_bilinear_mat(p,k);
                MatrixXd Vq(8,1);
                Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
                    V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;
                if(!check(K*Vq,m,n)) 
                {
                    // MatrixXd fail = K*Vq; 
                    // cout << m<< ' ' << n << endl;
                    // cout << fail(1) << ' ' << fail(0) << endl;
                    continue;
                }

                MatrixXd dec(2,4); 
                dec << 1,0,-1,0,
                        0,1,0,-1;
                // MatrixXd K = dec * (K1 - K2);
                K = dec * K;
                
                MatrixXd out = C * K;
                // output(out,"out");
                for(int x = 0;x < 2; ++ x)
                {
                    line_energy.insert(cnt,(i*21+j)*2) = out(x,1);
                    line_energy.insert(cnt,(i*21+j)*2+1) = out(x,0);
                    line_energy.insert(cnt,(i*21+j+1)*2) = out(x,3);
                    line_energy.insert(cnt,(i*21+j+1)*2+1) = out(x,2);
                    line_energy.insert(cnt,((i+1)*21+j)*2) = out(x,5);
                    line_energy.insert(cnt,((i+1)*21+j)*2+1) = out(x,4);
                    line_energy.insert(cnt,((i+1)*21+j+1)*2) = out(x,7);
                    line_energy.insert(cnt,((i+1)*21+j+1)*2+1) = out(x,6);
                    ++ cnt;
                }
                // exit(0);
            }
        }
    }
    SparseMatrix<double> final_line_energy(cnt,line_energy.cols());
    for (int k=0; k< line_energy.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(line_energy,k); it; ++it)
            final_line_energy.insert(it.row(), it.col()) = it.value();
    // cout << "solved line energy " << endl;
    final_line_energy.makeCompressed();
    return final_line_energy;
}

// SparseMatrix<double> energy::line_energy(vector<vector<vector<Line> > > mesh_line,vector<vector<point > > mesh,vector<vector<point > > V,double *bins,int num)
// {
//     // cout << "line energy solve" << endl;
//     // cout << num << endl;
//     SparseMatrix<double> line_energy(num * 2,21 * 21 * 2);
//     vector<point> p;
//     int cnt = 0;
//     double cost = 0;
//     for(int i = 0;i < 20;++ i)
//     {
//         for(int j = 0;j < 20;++ j)
//         {
//             for(auto k : mesh_line[i][j]) 
//             {
//                 if(k.solve == 0) continue; 
//                 // cout << "line energy i,j: " << i << ' ' << j << endl;
//                 MatrixXd R(2,2);
//                 R << cos(bins[k.pos]),-sin(bins[k.pos]),
//                         sin(bins[k.pos]),cos(bins[k.pos]);
//                 // R << cos(k.theta),sin(k.theta),
//                 //         -sin(k.theta),cos(k.theta);
//                 // R << 1,0,
//                 //     0,1;
//                 // cout << "line energy R: " << R.rows() << ' ' << R.cols() << endl;
//                 // cout <<"belong bins and theta: " << k.pos << ' ' << bins[k.pos] << ' '<< acos((k.x.x - k.y.x)/sqrt((k.x.y - k.y.y)*((k.x.y - k.y.y))+((k.x.x - k.y.x))*((k.x.x - k.y.x))))<< endl;
//                 // output(R,"R");
//                 // cout << "line e: \n" << k.x.x << ' ' << k.x.y << endl << k.y.x<< ' ' << k.y.y<< endl;
//                 MatrixXd e(2,1);
//                 e << k.x.x - k.y.x,
//                         k.x.y - k.y.y;
                    
//                 // cout << "line energy e: " << e.rows() << ' ' << e.cols() << endl;
//                 // output(e,"e");
//                 // output(R*e,"e'");
//                 // output(R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose(),"e''");
//                 MatrixXd C = R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose() - MatrixXd::Identity(2,2);
//                 // cout << "line energy C: " << C.rows() << ' ' << C.cols() << endl;
//                 // output(C,"C");
//                 // output(C*e,"Ce");
                
//                 p.clear();
//                 p.push_back(mesh[i][j]); p.push_back(mesh[i][j+1]);
//                 p.push_back(mesh[i+1][j]); p.push_back(mesh[i+1][j+1]);

//                 MatrixXd K = get_bilinear_mat(p,k);
//                 // p.clear();
//                 // p.push_back(V[i][j]); p.push_back(V[i][j+1]);
//                 // p.push_back(V[i+1][j]); p.push_back(V[i+1][j+1]);
//                 // MatrixXd K2 = get_bilinear_mat(p,k);

//                 // output(K1,"K1");
//                 // output(K2,"K2");
//                 MatrixXd Vq(8,1);
//                 // Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
//                 //     V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;
//                 Vq << mesh[i][j].y,mesh[i][j].x,mesh[i][j+1].y,mesh[i][j+1].x,
//                     mesh[i+1][j].y,mesh[i+1][j].x,mesh[i+1][j+1].y,mesh[i+1][j+1].x;
                    
//                 // MatrixXd F = kx.colPivHouseholderQr().solve(Vq);
//                 // cout << "K * VQ : "<< endl << K * Vq << endl;
//                 // cout << k.x.x << ' ' << k.x.y << endl;
//                 // cout << k.y.x << ' ' << k.y.y << endl;
//                 // cout << "Linear:"<< endl << x * F << endl;
//                 // exit(0);
//                 // VectorXd Z(8),B(2); 
//                 // Z << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x; 
//                 // cout << "K1 * Vq" << endl;
//                 // cout << K * Vq << endl;
//                 // cout << "-----S------" <<endl;
//                 // B << k.x.x , k.x.y;
//                 // cout << (K1 * Z - B).norm()<<endl;
//                 // if()
//                 // if((K * Vq - e1 * F).norm() > 0.1)
//                 //     continue;
//                 // show_new_line(K,V,i,j);
//                 MatrixXd dec(2,4); 
//                 dec << 1,0,-1,0,
//                         0,1,0,-1;
                        
//                 K = dec * K;
                
//                 // output(K*Vq,"Ksss");
//                 // output(F,"F");

//                 MatrixXd out = C * K;
//                 // [2,4] * [4,8] * [8,8] * [8,1] = [2,1]
//                 // MatrixXd Z = out * Vq;
//                 // output(Z,"Zzzz");
//                 // MatrixXd D = R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose()*e;
//                 // output(D,"ddd");
//                 // cout << endl << endl;

//                 // cout << "----cost-----"<<endl << out<<endl;
//                 // cost += out(0,0) * Vq(1) + out(1,1) * Vq(0) + out(2,2) * Vq(3) + out(3,3) * Vq(2) + 
//                 //         out(4,4) * Vq(5) + out(5,5) * Vq(4) + out(6,6) * Vq(7) + out(7,7) * Vq(6);
//                 // cout << cost << endl;
                
//                 // cout << "line energy out: " << out.rows() << ' ' << out.cols() << endl;
//                 // output(out,"out");
//                 // cout << cnt << ' '<<num << endl;
//                 for(int x = 0;x < 2; ++ x)
//                 {
//                     line_energy.insert(cnt,(i*21+j)*2) = out(x,1);
//                     line_energy.insert(cnt,(i*21+j)*2+1) = out(x,0);
//                     line_energy.insert(cnt,(i*21+j+1)*2) = out(x,3);
//                     line_energy.insert(cnt,(i*21+j+1)*2+1) = out(x,2);
//                     line_energy.insert(cnt,((i+1)*21+j)*2) = out(x,5);
//                     line_energy.insert(cnt,((i+1)*21+j)*2+1) = out(x,4);
//                     line_energy.insert(cnt,((i+1)*21+j+1)*2) = out(x,7);
//                     line_energy.insert(cnt,((i+1)*21+j+1)*2+1) = out(x,6);
//                     ++ cnt;
//                 }
//                 // exit(0);
//             }
//         }
//     }
//     // cout << "solved line energy " << endl;
//     line_energy.makeCompressed();
//     return line_energy;
// }
