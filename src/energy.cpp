#include "energy.h"


double energy::shape_energy(vector<vector<Point > > mesh,vector<vector<Point > > V)
{
    double shape_energy = 0;
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
            MatrixXd Vq(8,1);
            Vq << V[i][j].y,V[i][j].x,V[i][j+1].y,V[i][j+1].x,
                V[i+1][j].y,V[i+1][j].x,V[i+1][j+1].y,V[i+1][j+1].x;

            MatrixXd es = (Aq * (Aq.transpose() * Aq).inverse() * Aq.transpose() - MatrixXd::Identity(8,8)) * Vq;
            
            for(int x = 0;x < 8; ++ x)
                shape_energy += es(x,0);
        }
    // cout << shape_energy << endl;
	// shape_energy.makeCompressed();
    return shape_energy;
}

double energy::line_energy(vector<vector<vector<Line> > > mesh,vector<vector<Point > > V)
{

}