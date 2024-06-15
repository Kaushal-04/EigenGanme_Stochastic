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
MatrixXf normalizeVector(const MatrixXf& vect) {
    MatrixXf normalizedVector = vect;
    float maxElement = vect.cwiseAbs().maxCoeff();
    for (int i = 0; i < vect.size(); ++i) {
        normalizedVector(i) /= maxElement;
    }
    return normalizedVector;
}
void insertColumnMatrix(MatrixXf& originalMatrix, const MatrixXf& columnMatrix) {
    if (columnMatrix.rows() != originalMatrix.rows()) {
        cerr << "Dimention of Column & Rows not match" <<endl;
        return;
    }
    originalMatrix.conservativeResize(originalMatrix.rows(), originalMatrix.cols() + columnMatrix.cols());
    originalMatrix.rightCols(columnMatrix.cols()) = columnMatrix;
}
void deleteColumn(MatrixXf& matrix, int columnIndex) {
    if (columnIndex < 0 || columnIndex >= matrix.cols()) {
        cerr << "Error: Column index out of bounds." << endl;
        return;
    }
    for (int i = columnIndex; i < matrix.cols() - 1; ++i) {
        matrix.col(i) = matrix.col(i + 1);
    }
    matrix.conservativeResize(matrix.rows(), matrix.cols() - 1);
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

    
    MatrixXf EigenVector(n,1);
    MatrixXf temp;
    for(int i=0; i<eigenvectors.cols(); i++){
        temp = normalizeVector(eigenvectors.col(i));
        insertColumnMatrix(EigenVector , temp);
        deleteColumn(temp , 0);
    }
    deleteColumn(EigenVector , 0);
    cout<<"Normalized Eigen Vector\n"<<EigenVector<<endl;

    return 0;
}
