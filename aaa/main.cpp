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

void init(int*op)
{
    op[0] = 3; // pic id
    op[1] = 0; // show seam carving
    op[2] = 1; // show mesh line
    op[3] = 0; // show segment
    op[4] = 0; // segment length limit
    op[5] = 0; // show mesh segment
    op[6] = 1; // show energy cost
    op[7] = 0; // show shape cost
    op[8] = 0; // show line cost
    op[9] = 0; // show final mesh and line
}

int main()
{
    int*op = new int(10);
    init(op);
    double st = clock();
    // Mat img = imread("input/container.jpg");
    Mat img = imread("input/"+ to_string(op[0])+ "_input.jpg"),line_img;
    // namedWindow("origin image",WINDOW_FREERATIO);
    // imshow("origin image", img);
    SeamC::init(img,op[1]); //op
    img = SeamC::get_rec_img();

    double st1 = clock();
    cout << fixed << setprecision(3)<< "seam carving cost time: " << (st1 - st) / CLOCKS_PER_SEC << endl; 

    vector<vector<point > > U = SeamC::get_U(),pos;
    // waitKey(0);
    vector<vector<point > > mesh = Mesh_Placement(SeamC::get_seam_carving(),img,U,pos,op[2]); // 网格坐标 op
    double st2 = clock();
    cout << fixed << setprecision(3)<< "mesh cost time: " << (st2 - st1) / CLOCKS_PER_SEC << endl; 

    int *n_line = new int,*Nl = new int; // 线段数量
    double* line = fline::line_img(img,n_line,line_img,op[3]); // 线段 and op
    vector<vector<vector<Line > > > mesh_line = fline::init_line(mesh,line,*n_line,Nl,op[4]); // 每个网格中的线段
    
    if(op[5]) fline::check_line(img,mesh_line);// op

    int M = 50; // 桶的数量
    int*num = new int[M];
    double *bins;
    bins = new double[M];
    fline::init_bins(num,mesh_line,M);
    energy::Line_rotate_count(bins,num,mesh_line,mesh,pos,M);

    double Nq = 20 * 20,lambdab = 1e6,lambdal=100;
    SparseMatrix<double> ES = (1.0 / Nq) * energy::shape_energy(mesh);
    // SparseMatrix<double> ES = energy::shape_energy(mesh);
    pair<SparseMatrix<double>,MatrixXd > EB = energy::bound_energy(mesh,lambdab,img.rows,img.cols);

    // VectorXd V = get_mesh_point(pos);
    double st3 = clock();
    cout << fixed << setprecision(3)<< "shape and bound cost time: " << (st3 - st2) / CLOCKS_PER_SEC << endl; 

    double lstcost = 1000000,cost;

    for(int iter = 0;iter < 10; ++ iter)
    {
        SparseMatrix<double> EL = lambdal*(1.0 / *Nl) * energy::line_energy(mesh_line,mesh,pos,bins,*Nl,img.rows,img.cols);
        // cout << "line energy : "<< EL.rows() << " " << 8 * (*Nl) << endl;

        VectorXd B = VectorXd::Zero(ES.rows() + EL.rows());
        // VectorXd B = VectorXd::Zero(8 * (*Nl));
        // VectorXd B = VectorXd::Zero(20 * 20 * 8);
        // cout << "qwq:" << endl << EB.second.rows() << ' ' << B.rows() << endl;
        SparseMatrix<double> A = merge_matrix(merge_matrix(ES,EL),EB.first);
        // SparseMatrix<double> A = merge_matrix(EL,EB.first);
        // SparseMatrix<double> A = ES;
        // SparseMatrix<double> A = EB.first;A.makeCompressed();
        // SparseMatrix<double> A = merge_matrix(ES,EB.first);A.makeCompressed();
        // output(A,"Matrix A");
        B = merge_matrix(B,EB.second);
        // B = EB.second;
        
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        solver.compute(A);
        VectorXd V = solver.solve(B);
        // V = get_mesh_point(mesh);
        SparseMatrix<double> Vq(21*21*2,1),Wq(B.rows(),1);
        for(int i = 0;i < V.rows();++ i) Vq.insert(i,0) = V(i);
        for(int i = 0;i < B.rows();++ i) Wq.insert(i,0) = B(i);
        
        SparseMatrix<double> energy_cost;
        energy_cost = (A * Vq) - Wq;
        energy_cost = energy_cost.transpose()*energy_cost;
        cost = energy_cost.coeffRef(0,0);
        if(op[6])
        {
            cout << "energy cost: " << endl << energy_cost.coeffRef(0,0) << endl;
        }
        if(op[7])
        {
            SparseMatrix<double> shape_cost;
            shape_cost = ES * Vq;
            shape_cost = shape_cost.transpose() * shape_cost;
            cout << "shape cost: " << endl << shape_cost.coeffRef(0,0) << endl;
        }
        if(op[8])
        {
            SparseMatrix<double> line_cost;
            line_cost = EL * Vq;
            line_cost = line_cost.transpose() * line_cost;
            cout << "line cost: " << endl << line_cost.coeffRef(0,0) << endl;
        }

        // cout << "vec_to_mesh" << endl;
        pos = vec_to_mesh(V);
        if(abs(cost - lstcost) < 0.01) break;
        else lstcost = cost;
        // cout << "recompute Line" << endl;
        energy::Line_rotate_count(bins,num,mesh_line,mesh,pos,M);
    }
    double st4 = clock();
    cout << fixed << setprecision(3)<< "energy cost time: " << (st4 - st3) / CLOCKS_PER_SEC << endl; 

    if(op[9])
    {
        img_mesh(line_img,mesh);
        GLout(line_img,mesh,pos);
    }
    cout << fixed << setprecision(3)<< "total cost time: " << (st4 - st) / CLOCKS_PER_SEC << endl; 
    GLout(img,mesh,pos);
    waitKey(0);
    return 0;
}