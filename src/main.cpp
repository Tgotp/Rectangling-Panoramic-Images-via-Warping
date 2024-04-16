#include <cstdio>
#include <iostream> 

#include <iomanip>
#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "seam_carving.h"
#include "mesh_placement.h"
#include "energy.h"
#include "fLine.h"
#include <Eigen/Eigen>
#include "GLout.h"
#include "tool.h"
using namespace cv;
using namespace std;


int main()
{
    double st = clock();
    Mat img = imread("input/3_input.jpg");
    // namedWindow("origin image",WINDOW_FREERATIO);
    // imshow("origin image", img);
    SeamC::init(img,0); //op
    img = SeamC::get_rec_img();

    double st1 = clock();
    cout << fixed << setprecision(3)<< "seam carving cost time: " << (st1 - st) / CLOCKS_PER_SEC << endl; 

    vector<vector<Point > > U = SeamC::get_U(),pos;
    // waitKey(0);
    vector<vector<Point > > mesh = Mesh_Placement(SeamC::get_seam_carving(),img,U,pos,0); // 网格坐标 op
    double st2 = clock();
    cout << fixed << setprecision(3)<< "mesh cost time: " << (st2 - st1) / CLOCKS_PER_SEC << endl; 

    int *n_line = new int,*Nl = new int; // 线段数量
    double* line = fline::line_img(img,n_line,0); // 线段 and op
    vector<vector<vector<Line > > > mesh_line = fline::init_line(mesh,line,*n_line,Nl); // 每个网格中的线段
    // fline::check_line(img,mesh_line);// op

    int M = 50; // 桶的数量
    int*num = new int[M];
    double*bins = energy::Line_rotate_count(num,mesh_line,mesh,pos,M);

    double Nq = 20 * 20,lambdab = 1e5;
    // SparseMatrix<double> ES = (1.0 / Nq) * energy::shape_energy(mesh);
    SparseMatrix<double> ES = energy::shape_energy(mesh);
    pair<SparseMatrix<double>,MatrixXd > EB = energy::bound_energy(mesh,lambdab,img.rows,img.cols);

    // VectorXd V = get_mesh_point(pos);

    for(int iter = 0;iter < 1; ++ iter)
    {
        cout << 8*(*Nl) << endl;
        SparseMatrix<double> EL = 100 * (1.0 / *Nl) * energy::line_energy(mesh_line,mesh,bins,*Nl);
        cout << "line energy : "<< EL.rows() << " " << 8 * (*Nl) << endl;

        VectorXd B = VectorXd::Zero(20 * 20 * 8 + 4 * (*Nl));
        // VectorXd B = VectorXd::Zero(20 * 20 * 8);
        // cout << "qwq:" << endl << EB.second.rows() << ' ' << B.rows() << endl;
        SparseMatrix<double> A = merge_matrix(merge_matrix(ES,EL),EB.first);
        // SparseMatrix<double> A = EB.first;A.makeCompressed();
        // SparseMatrix<double> A = merge_matrix(ES,EB.first);A.makeCompressed();
        // output(A,"Matrix A");
        B = merge_matrix(B,EB.second);
        // B = EB.second;
        
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        solver.compute(A);
        VectorXd V = solver.solve(B);
        SparseMatrix<double> Vq(21*21*2,1),Wq(B.rows(),1);
        for(int i = 0;i < V.rows();++ i) Vq.insert(i,0) = V(i);
        for(int i = 0;i < B.rows();++ i) Wq.insert(i,0) = B(i);
        SparseMatrix<double> energy_cost;
        energy_cost = (A * Vq) - Wq;
        energy_cost = energy_cost.transpose()*energy_cost;
        cout << "enery cost: " << endl << energy_cost.coeffRef(0,0) << endl;

        pos = vec_to_mesh(V);
        double*bins = energy::Line_rotate_count(num,mesh_line,mesh,pos,M);
    }
    double st3 = clock();
    cout << fixed << setprecision(3)<< "energy cost time: " << (st3 - st2) / CLOCKS_PER_SEC << endl; 

    GLout(img,mesh,pos);
    waitKey(0);
    return 0;
}