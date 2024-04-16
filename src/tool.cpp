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
    cout << name << endl;
    cout << 'rows and cols : '<< x.rows() << ' ' << x.cols() << endl;
    for(int i = 0;i < x.rows(); ++ i)
    {
        for(int j = 0;j < x.cols(); ++ j)
            cout << x(i,j) * 10000<< ' ';
        cout << endl;
    }
}

void output(SparseMatrix<double> x,string name)
{
    cout << name << endl;
    cout << 'rows and cols : '<< x.rows() << ' ' << x.cols() << endl;
    for(int i = 0;i < x.rows(); ++ i)
    {
        for(int j = 0;j < x.cols(); ++ j)
            cout << x.coeffRef(i,j) * 10000<< ' ';
        cout << endl;
    }
}