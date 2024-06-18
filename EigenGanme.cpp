#include <iostream>
#include <Eigen/Dense>
#include <algorithm>
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
int calculateRank(const MatrixXf& matrix) {
    // Perform full pivoting QR decomposition
    FullPivHouseholderQR<MatrixXf> qr(matrix);
    // Rank is the number of non-zero pivots
    return qr.rank();
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
    MatrixXf w = B *EigenVector;
    int T=1; //3 for just check
    int M=2; //check it
    vector<MatrixXf> matA;
    MatrixXf tempMat;
    for(int i=0; i<M; i++){
        tempMat.setRandom(n,n);
        tempMat = tempMat.array().abs(); // Ensure positive values
        matA.push_back(tempMat);
    }
    vector<MatrixXf> matB;
    for(int i=0; i<M; i++){
        tempMat.setRandom(n,n);
        tempMat = tempMat.array().abs(); // Ensure positive values
        matB.push_back(tempMat);
    }
    MatrixXf Atm = MatrixXf::Zero(n,n);
    MatrixXf Btm = MatrixXf::Zero(n,n);
    int colNum;
    MatrixXf vi , vj;
    for(int j=0; j<T; j++){
        for(int i=0; i<n; i++){
            for(int m=0; m<M; m++){
                Atm += matA[m];
                Atm /= m+1;
                Btm += matB[m];
                Btm /= m+1;
                colNum = rand()%EigenVector.cols();
                vi = EigenVector.col(colNum);
                colNum = rand()%EigenVector.cols();
                vj = EigenVector.col(colNum);

                float result = (vj.transpose() * B * vj)(0, 0);
                if(result < 0) //to fin nan result
                    result = result * (-1);
                float rank = calculateRank(A);
                float maxi = max(result , rank);
                float sqrtResult = sqrt(maxi);
                sqrtResult = 1.0 /sqrtResult;
                MatrixXf yj;
                yj = vj / sqrtResult;
            }
            cout<<"Atm\n"<<Atm<<endl;
            cout<<"Btm\n"<<Btm<<endl;         
        }
    }

    return 0;
}
