#include "tool.h"

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
void output(MatrixXd x,string name)
{
    cout << "--------"<< name << "--------"<< endl;
    cout << "rows and cols : "<< x.rows() << ' ' << x.cols() << endl;
    for(int i = 0;i < x.rows(); ++ i)
    {
        for(int j = 0;j < x.cols(); ++ j)
            cout << fixed << setprecision(16) << x(i,j)<< ' ';
        cout << endl;
    }
}

void output(SparseMatrix<double> x,string name)
{
    cout << name << endl;
    cout << "rows and cols :"<< x.rows() << ' ' << x.cols() << endl;
    for(int i = 0;i < x.rows(); ++ i)
    {
        for(int j = 0;j < x.cols(); ++ j)
            cout << fixed << setprecision(16) << x.coeffRef(i,j)<< ' ';
        cout << endl;
    }
}

void show_new_line(MatrixXd K,vector<vector<point > > V,int i,int j)
{
    MatrixXd Vq(8,1);
    Vq(0) = V[i][j].x,Vq(1) = V[i][j].y;
    Vq(2) = V[i][j+1].x,Vq(3) = V[i][j+1].y;
    Vq(4) = V[i+1][j].x,Vq(5) = V[i+1][j].y;
    Vq(6) = V[i+1][j+1].x,Vq(7) = V[i+1][j+1].y;
    
    MatrixXd new_V;
    new_V = K * Vq;
    cout << "new e: " << new_V.rows() << ' ' << new_V.cols() << endl
                        << new_V(0) << ' ' << new_V(1) << endl 
                        << new_V(2) << ' ' << new_V(3)<< endl
                        << new_V(0) - new_V(1) << endl 
                        << new_V(2) - new_V(3)<< endl;
}