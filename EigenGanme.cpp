#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
#define MAX_SIZE 10000
int n;
void setValueOfn(){
    cout<<"Enter order of matrix:";
    cin>>n;
    if (n > MAX_SIZE) {
        cout << "Error: Order exceeds maximum size. Setting order to maximum size.\n";
        n = MAX_SIZE;
    }
}
MatrixXf makeSymmetric(const MatrixXf& A) {
    return 0.5 * (A + A.transpose());
}
int main() {
     setValueOfn();
    // Equation Ax = nBx
    MatrixXf randomMatrixA , randomMatrixB;
    randomMatrixA.setRandom(n,n);
    randomMatrixA = randomMatrixA.array().abs(); // Ensure positive values
    randomMatrixB.setRandom(n,n);
    randomMatrixB = randomMatrixB.array().abs(); // Ensure positive values
    MatrixXf A = makeSymmetric(randomMatrixA);
    MatrixXf B = makeSymmetric(randomMatrixB);
          // Eigenvalue and Eigenvector Computation: The GeneralizedEigenSolver class in Eigen is used
          // to solve the generalized eigenvalue problem Ax=λBx. The .compute(A, B) method computes the 
          // eigenvalues and eigenvectors.
    GeneralizedEigenSolver<MatrixXf> ges;
    ges.compute(A, B);
    MatrixXf eigenvalues = ges.eigenvalues().real();
    MatrixXf eigenvectors = ges.eigenvectors().real();
    cout << "Eigenvalues: \n" << eigenvalues << endl;
    cout << "Eigenvectors: \n" << eigenvectors << endl;
    return 0;
}
