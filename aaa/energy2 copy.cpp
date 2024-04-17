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
    cout << "solved shape energy" << endl;
	shape_energy.makeCompressed();
    return shape_energy;
}

pair<SparseMatrix<double>,VectorXd > energy::bound_energy(vector<vector<Point > > mesh,double inf,int n,int m)
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
    cout << "solved bound energy" << endl;
    return bound_energy;
}

MatrixXd inverse_bilinear_interpolation_function(vector<Point> P[2],Line X)
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

Line inverse_bilinear_interpolation(vector<Point> P[2],Line X,Point st)
{
    // cout <<"invers V :  " << P[0] << ' ' << P[1] << endl;
    MatrixXd F = inverse_bilinear_interpolation_function(P,X);
    MatrixXd x(2,8),y(2,8),z(2,8);
    x << X.x.x,X.x.y,X.x.x*X.x.y,1,0,0,0,0,
        0,0,0,0,X.x.x,X.x.y,X.x.x*X.x.y,1;
    y << X.y.x,X.y.y,X.y.x*X.y.y,1,0,0,0,0,
        0,0,0,0,X.y.x,X.y.y,X.y.x*X.y.y,1;
    z << st.y,st.x,st.x*st.y,1,0,0,0,0,
        0,0,0,0,st.y,st.x,st.x*st.y,1;
    MatrixXd a = x * F;
    MatrixXd b = y * F;
    MatrixXd c = z * F;
    // double q = a(0) - b(0);
    // double e = a(1) - b(1);
    return Line(Point(a(0)-c(0),a(1)-c(1)),Point(b(0)-c(0),b(1)-c(1)));
    // if(q == 0 && e == 0) return 4;
    // // cout << X.x << ' ' << X.y << ' ' << q << ' ' << e << endl;
    // // cout <<"----a----"<< endl << a << endl;
    // // cout <<"----b----"<< endl << b << endl;
    // // cout <<"theta :" << q/ sqrt(q*q + e*e) << " " << acos((q)/ sqrt(q*q + e*e))<< endl;
    // // cout << acos(q / sqrt(q*q + e*e)) << endl;
    // return acos((q)/ sqrt(q*q + e*e));
}
vector<vector<vector<Line> > > qwq;
void energy::Line_rotate_count(double *bins,int *num,vector<vector<vector<Line > > >&mesh_line,vector<vector<Point> >mesh,vector<vector<Point> >V,int number)
{
    qwq.clear();
    vector<Line> q2q;
    vector<vector<Line> > qaq;
    for(int i = 0;i < number;++ i)
        bins[i] = 0,num[i] = 0;
    vector<Point> p[2];
    for(int i = 0;i < 20;++ i)
    {
        qaq.clear();
        for(int j = 0;j < 20;++ j)
        {
            q2q.clear();
            for(auto&k: mesh_line[i][j])
            {
                p[0].clear();p[1].clear();
                p[0].push_back(mesh[i][j]); p[0].push_back(mesh[i][j+1]);
                p[0].push_back(mesh[i+1][j+1]); p[0].push_back(mesh[i+1][j]);
                p[1].push_back(V[i][j]); p[1].push_back(V[i][j+1]);
                p[1].push_back(V[i+1][j+1]); p[1].push_back(V[i+1][j]); 
                Line t = inverse_bilinear_interpolation(p,k,mesh[i][j]);
                q2q.push_back(t);
                // cout <<"line rotate theta : " << theta << ' ' << acos((k.x.x - k.y.x)/sqrt((k.x.y - k.y.y)*((k.x.y - k.y.y))+((k.x.x - k.y.x))*((k.x.x - k.y.x)))) << endl;
                double theta;

                double q = t.x.x - t.y.x;
                double e = t.x.y - t.y.y;
                if(q == 0 && e == 0) theta = 4;
                // cout << X.x << ' ' << X.y << ' ' << q << ' ' << e << endl;
                // cout <<"----a----"<< endl << a << endl;
                // cout <<"----b----"<< endl << b << endl;
                // cout <<"theta :" << q/ sqrt(q*q + e*e) << " " << acos((q)/ sqrt(q*q + e*e))<< endl;
                // cout << acos(q / sqrt(q*q + e*e)) << endl;
                else theta = acos((q)/ sqrt(q*q + e*e));

                if(theta == 4)
                {
                    k.pos = 50;
                    continue;
                }
                k.pos = floor(theta/PI * 50);
                k.pos = min(49,k.pos);
                // cout <<"theta and pos : "<< k.theta << ' ' << k.pos << endl;
                bins[k.pos] += theta;
                num[k.pos] ++;
            }
            qaq.push_back(q2q);
        }
        qwq.push_back(qaq);
    }
    for(int i = 0;i < number; ++ i)
        if(num[i]) bins[i] /= num[i];
    // cout << "-----bins----" << endl;
    // for(int i = 0;i < number; ++ i)
    //     cout << bins[i] << ' '<< endl;
    // cout << "Line rotate count solved out" << endl;
}

MatrixXd get_bilinear_mat(vector<Point> p,Line L)
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

    MatrixXd dec(2,4); 
    dec << -1,0,1,0,
            0,-1,0,1;
                        
    K = dec * K;
    // cout << "solved get bilinear mat" << endl;
    return K;
}

SparseMatrix<double> energy::line_energy(vector<vector<vector<Line> > > mesh_line,vector<vector<Point > > mesh,vector<vector<Point > > V,double *bins,int num)
{
    cout << "line energy solve" << endl;
    // cout << num << endl;
    SparseMatrix<double> line_energy(num * 2,21 * 21 * 2);
    vector<Point> p;
    int cnt = 0; int aaa = 0; 
    double cost = 0;
    for(int i = 0;i < 20;++ i)
    {
        for(int j = 0;j < 20;++ j)
        {
            aaa = -1;
            for(auto k : mesh_line[i][j]) 
            {
                ++ aaa;
                if(k.pos == 50) continue; 
                // cout << "line energy i,j: " << i << ' ' << j << endl;
                MatrixXd R(2,2);
                R << cos(bins[k.pos]),-sin(bins[k.pos]),
                        sin(bins[k.pos]),cos(bins[k.pos]);
                // R << cos(bins[k.pos]),sin(bins[k.pos]),
                //         -sin(bins[k.pos]),cos(bins[k.pos]);
                // cout << "line energy R: " << R.rows() << ' ' << R.cols() << endl;
                // cout <<"belong bins and theta: " << k.pos << ' ' << bins[k.pos] << ' '<< acos((k.x.x - k.y.x)/sqrt((k.x.y - k.y.y)*((k.x.y - k.y.y))+((k.x.x - k.y.x))*((k.x.x - k.y.x))))<< endl;
                // output(R,"R");
                // cout << "line e: \n" << k.x.x << ' ' << k.x.y << endl << k.y.x<< ' ' << k.y.y<< endl;
                MatrixXd e(2,1);
                e << k.x.x - k.y.x,
                        k.x.y - k.y.y;
                    
                // cout << "line energy e: " << e.rows() << ' ' << e.cols() << endl;
                // output(e,"e");
                // output(R*e,"e'");
                // output(R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose(),"e''");
                MatrixXd C = R*e*(e.transpose()*e).inverse()*e.transpose()*R.transpose() - MatrixXd::Identity(2,2);
                // cout << "line energy C: " << C.rows() << ' ' << C.cols() << endl;
                // output(C,"C");
                // output(C*e,"Ce");
                
                p.clear();
                p.push_back(mesh[i][j]); p.push_back(mesh[i][j+1]);
                p.push_back(mesh[i+1][j]); p.push_back(mesh[i+1][j+1]);

                MatrixXd K1 = get_bilinear_mat(p,Line(mesh[i][j],k.x));
                MatrixXd K2 = get_bilinear_mat(p,Line(mesh[i][j],k.y));

                // output(K1,"K1");
                // output(K2,"K2");
                // MatrixXd Vq(8,1);
                // Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
                //     V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;
                    
                // MatrixXd F = kx.colPivHouseholderQr().solve(Vq);
                // cout << "K*VQ : "<< endl << K * Vq << endl;
                // cout << "Linear:"<< endl << x * F << endl;
                // exit(0);
                VectorXd Z(8),B(2); 
                Z << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x; 
                cout << "K1 * Z" << endl;
                cout << K1 * Z << endl;
                cout << "-----S------" <<endl;
                cout << qwq[i][j][aaa].x.x - mesh[i][j].y << ' ' << qwq[i][j][aaa].x.y - mesh[i][j].x << endl;
                cout << qwq[i][j][aaa].y.x - mesh[i][j].y << ' ' << qwq[i][j][aaa].y.y - mesh[i][j].x << endl;
                cout << qwq[i][j][aaa].x.x << ' ' << qwq[i][j][aaa].x.y << endl;
                cout << qwq[i][j][aaa].y.x << ' ' << qwq[i][j][aaa].y.y << endl;
                B << k.x.x , k.x.y;
                // cout << (K1 * Z - B).norm()<<endl;
                // if()
                // if((K * Vq - e1 * F).norm() > 0.1)
                //     continue;
                // show_new_line(K,V,i,j);
                // MatrixXd dec(2,4); 
                // dec << 1,0,-1,0,
                //         0,1,0,-1;
                        
                // K = dec * K;
                
                // output(K,"K");
                // output(K,"K");
                // output(F,"F");

                MatrixXd out = C * (K1 - K2);
                // [2,4] * [4,8] * [8,8] * [8,1] = [2,1]

                // cout << "----cost-----"<<endl << out<<endl;
                // cost += out(0,0) * Vq(1) + out(1,1) * Vq(0) + out(2,2) * Vq(3) + out(3,3) * Vq(2) + 
                //         out(4,4) * Vq(5) + out(5,5) * Vq(4) + out(6,6) * Vq(7) + out(7,7) * Vq(6);
                // cout << cost << endl;
                
                // cout << "line energy out: " << out.rows() << ' ' << out.cols() << endl;
                // output(out,"out");
                // cout << cnt << ' '<<num << endl;
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
                exit(0);
            }
        }
    }
    cout << "solved line energy " << endl;
    line_energy.makeCompressed();
    return line_energy;
}
