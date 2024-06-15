#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main() {
    int n= 3;
    // Equation Ax = nBx
    MatrixXf A(n, n);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    MatrixXf B(n, n);
    B << 9, 8, 7,
         6, 5, 4,
         3, 2, 1;
    GeneralizedEigenSolver<MatrixXf> ges;
    ges.compute(A, B);
    MatrixXf eigenvalues = ges.eigenvalues().real();
    MatrixXf eigenvectors = ges.eigenvectors().real();
    cout << "Eigenvalues: \n" << eigenvalues << endl;
    cout << "Eigenvectors: \n" << eigenvectors << endl;
    return 0;
}
