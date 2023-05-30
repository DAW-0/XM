
/*****************************************************************************************************************************/

/*
   LibMatrix.mqh
   Copyright © 2009 Logunov Eugene (lea)
   rootrat@mail.ru
*/

#property copyright "Copyright © 2009 Logunov Eugene (lea)"
#property link "rootrat@mail.ru"

/*****************************************************************************************************************************/

// useful mathematical constants
#define ZERO 0.0
#define ONE 1.0
#define MINUS_ONE -1.0
#define DEFAULT_TOLERANCE 0.000000001

/*****************************************************************************************************************************/

#import "LibMatrix.ex4"

/*****************************************************************************************************************************/

   /*
      abstract:
         linear interpolation between rangeLowLimit and rangeHighLimit
      restrictions:
         balance should be in range [0;1] (interpolation)
         balance is allowed to be out of range for linear extrapolation purposes olny
   */
   double MathLerp(double rangeLowLimit, double rangeHighLimit, double balance);

   /*
      abstract:
         returns a random number in range [rangeLowLimit;rangeHighLimit]
      restrictions:
         none
   */
   double MathInRangeRandom(double rangeLowLimit, double rangeHighLimit);

   /*
      abstract:
         checks two real numbers for equality with given tolerance
      restrictions:
         tolerance should be greater or equal zero
   */
   bool MathDoublesEqual(double value1, double value2, double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         converts 2d matrix indicies to 1d array offset
         should be used to access matrix elements via indicies
      restrictions:
         row should be in [0;numRows-1]
         col should be in [0;numCols-1]
   */
   int MatrIndiciesToOffset(int row, int col, int numRows, int numCols);

   /*
      abstract:
         copies src matrix to dst
      restrictions:
         none
   */
   void MatrCopy(double& src[], double& dst[]);

   /*
      abstract:
         resizes 1d array so that matrix of given size could be stored in it
      restrictions:
         numRows and numCols should be positive
   */
   void MatrSetSize(double& matr[], int numRows, int numCols);

   /*
      abstract:
         resizes given matrix without data loss
         new elements are set to zero
         old elements may be destroyed if one of dimensions is resized to a smaller value
      restrictions:
         none
   */
   void MatrResize(double& matr[], int numRowsOld, int numColsOld, int numRowsNew, int numColsNew);

/*****************************************************************************************************************************/

   /*
      abstract:
         swaps two rows (row1 and row2) of matrix
      restrictions:
         row1 and row2 should be in range [0;numRows-1]
   */
   void MatrSwapRows(double& matr[], int numRows, int numCols, int row1, int row2);

   /*
      abstract:
         swaps two columns (col1 and col2) of matrix
      restrictions:
         col1 and col2 should be in range [0;numCols-1]
   */
   void MatrSwapCols(double& matr[], int numRows, int numCols, int col1, int col2);

   /*
      abstract:
         copy row (replace one of rows with another existing)
      restrictions:
         rowToCopy and rowToReplace should be in range [0;numRows-1]
   */
   void MatrCopyRow(double& matr[], int numRows, int numCols, int rowToCopy, int rowToReplace);

   /*
      abstract:
         copy column (replace one of cols with another existing)
      restrictions:
         colToCopy and colToReplace should be in range [0;numCols-1]
   */
   void MatrCopyCol(double& matr[], int numRows, int numCols, int colToCopy, int colToReplace);

   /*
      abstract:
         check whether row is zero with given tolerance
      restrictions:
         row should be in range [0;numRows-1]
   */
   bool MatrRowIsZero(double& matr[], int numRows, int numCols, int row, double tolerance);

   /*
      abstract:
         check whether column is zero with given tolerance
      restrictions:
         col should be in range [0;numCols-1]
   */
   bool MatrColIsZero(double& matr[], int numRows, int numCols, int col, double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         check whether matrix is square
      restrictions:
         none
   */
   bool MatrIsSquare(int numRows, int numCols);

   /*
      abstract:
         check whether matrix element belongs to main diagonal
      restrictions:
         none
   */
   bool MatrIsElemOnMainDiagonal(int row, int col);

   /*
      abstract:
         check whether two matricies are compatible for addition
      restrictions:
         none
   */
   bool MatrCompatiblityCheckAdd(int numRows1, int numCols1, int numRows2, int numCols2);

   /*
      abstract:
         check whether two matricies are compatible for multiplication
      restrictions:
         none
   */
   bool MatrCompatiblityCheckMul(int numRows1, int numCols1, int numRows2, int numCols2);

/*****************************************************************************************************************************/

   /*
      abstract:
         load zero matrix
      restrictions:
         none
   */
   void MatrLoadZero(double& matr[], int numRows, int numCols);

   /*
      abstract:
         load identity matrix
      restrictions:
         matrix should be square
   */
   void MatrLoadIdentity(double& matr[], int numRows, int numCols);

   /*
      abstract:
         load matrix with random numbers in range [rangeLowLimit;rangeHighLimit]
      restrictions:
         none
   */
   void MatrLoadInRangeRandom(double& matr[], int numRows, int numCols, double rangeLowLimit, double rangeHighLimit);

/*****************************************************************************************************************************/

   /*
      abstract:
         check whether matrix is zero using given comparison tolerance
      restrictions:
         tolerance should be greater or equal zero
   */
   bool MatrIsZero(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         check whether matrix is diagonal using given comparison tolerance
         matrix is assumed to be diagonal if all elements outside of main diagonal are zero
      restrictions:
         matrix should be square
         tolerance should be greater or equal zero
   */
   bool MatrIsDiagonal(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         check whether matrix is identity using given comparison tolerance
      restrictions:
         matrix should be square
         tolerance should be greater or equal zero
   */
   bool MatrIsIdentity(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         check whether matrix is symmetric (a[i][j]=a[j][i]) using given comparison tolerance
      restrictions:
         matrix should be square
         tolerance should be greater or equal zero
   */
   bool MatrIsSymmetric(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         check whether matrix is antisymmetric (a[i][j]=-a[j][i]; a[k][k]=0) using given comparison tolerance
      restrictions:
         matrix should be square
         tolerance should be greater or equal zero
   */
   bool MatrIsAntisymmetric(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         check whether matricies are equal
      restrictions:
         matricies should have same size
         tolerance should be greater or equal zero
   */
   bool MatrEqual(double& matr1[], double& matr2[], int numRows, int numCols, double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         adds scalar value to each element of matrix
      restrictions:
         none
   */
   void MatrAddScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]);

   /*
      abstract:
         subtracts scalar value from each element of matrix
      restrictions:
         none
   */
   void MatrSubScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]);

   /*
      abstract:
         multiplies each element of matrix by scalar value
      restrictions:
         none
   */
   void MatrMulByScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]);

   /*
      abstract:
         divides each element of matrix by scalar value
      restrictions:
         scalar should be non-zero
   */
   void MatrDivByScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]);

/*****************************************************************************************************************************/

   /*
      abstract:
         adds second matrix to first
      restrictions:
         matricies should have same size
   */
   void MatrAddMatr(double& matr1[], double& matr2[], int numRows, int numCols, double& result[]);

   /*
      abstract:
         subtracts second matrix from first
      restrictions:
         matricies should have same size
   */
   void MatrSubMatr(double& matr1[], double& matr2[], int numRows, int numCols, double& result[]);

   /*
      abstract:
         multiplies first matrix by second
      restrictions:
         matricies should be compatible for multiplication
   */
   void MatrMulMatr(
         double& matr1[], int numRows1, int numCols1,
         double& matr2[], int numRows2, int numCols2,
         double& result[], int& numRowsRes, int& numColsRes);

   /*
      abstract:
         computes trace (sum of diagonal elements) of given matrix
         source matrix is not modified
      restrictions:
         matrix should be square
   */
   double MatrTrace(double& matr[], int numRows, int numCols);

   /*
      abstract:
         performs matrix transpose
      restrictions:
         none
   */
   void MatrTranspose(double& matr[], int numRows, int numCols, double& result[], int& numRowsRes, int& numColsRes);

/*****************************************************************************************************************************/

   /*
      abstract:
         performs Gaussian triangulation of matrix
         returns size of subset of linear-independent rows (matrix rank)
      restrictions:
         tolerance should be greater or equal zero
   */
   int MatrGaussianElimination(double& matr[], int numRows, int numCols, double& result[], double tolerance);

   /*
      abstract:
         an implementation of Gauss-Jordan reduction for solving systems of simultaneous equations.
         usable for batch solving of systems with same matrix and different right-hand sides (RHS).
         right-hand sides should be represented be matrix, where each right-hand side is a column.
         output matrix (roots) will contain solutions in each column according to right-hand sides 
         matrix. if batch solution was found - method returns true; otherwise - false.
      restrictions:
         numRows should be equal numRowsRHS
   */
   bool MatrGJBatchSolve(
         double& matr[], int numRows, int numCols,
         double& rhs[], int numRowsRHS, int numColsRHS,
         double& roots[], int& numRowsRoots, int& numColsRoots,
         double tolerance
         );

/*****************************************************************************************************************************/

   /*
      abstract:
         extracts minor of given element
      restrictions:
         both numRows and numCols should be greater or equal 2
   */
   void MatrMinor(
         double& matr[], int numRows, int numCols,
         int rowToExclude, int colToExclude,
         double& result[], int& numRowsRes, int& numColsRes);

   /*
      abstract:
         computes algebraic complement of given element (a determinant of its minor with sign)
         Gaussian elimination with given pivoting tolerance is used to compute determinant
      restrictions:
         tolerance should be greater or equal zero
         matrix should be square
         both numRows and numCols should be greater or equal 2
   */
   double MatrAlgebraicComplement(double& matr[], int numRows, int numCols, int elemRow, int elemCol, double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         computes matrix inverse (complexity O(N^5))
         returns true if matrix inversion is complete
      restrictions:
         tolerance should be greater or equal zero
         matrix should be square
         both numRows and numCols should be greater or equal 2
         matrix should have non-zero determinant
   */
   bool MatrInvertUsingMinors(double& matr[], int numRows, int numCols, double& result[], double tolerance);

   /*
      abstract:
         computes matrix inverse (complexity O(N^3) - using Gauss-Jordan reduction)
         returns true if matrix inversion is complete
      restrictions:
         tolerance should be greater or equal zero
         matrix should be square
         both numRows and numCols should be greater or equal 2
         matrix should have non-zero determinant
   */
   bool MatrInvertUsingGJ(double& matr[], int numRows, int numCols, double& result[], double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         computes matrix determinant using Gaussian elimination with given pivoting tolerance
         decomposition throung minors is much slower
         source matrix is not modified
      restrictions:
         matrix should be square
         tolerance should be greater or equal zero
   */
   double MatrDet(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         computes determinant of triangulated matrix
         source matrix is not modified
      restrictions:
         matrix should be square
   */
   double MatrDetTriang(double& matr[], int numRows, int numCols);

/*****************************************************************************************************************************/

   /*
      abstract:
         computes matrix rank using Gaussian elimination with given pivoting tolerance
         source matrix is not modified
      restrictions:
         tolerance should be greater or equal zero
   */
   double MatrRank(double& matr[], int numRows, int numCols, double tolerance);

   /*
      abstract:
         computes rank of triangulated matrix
         source matrix is not modified
         tolerance is used when searching for zero row
      restrictions:
         tolerance should be greater or equal zero
   */
   double MatrRankTriang(double& matr[], int numRows, int numCols, double tolerance);

/*****************************************************************************************************************************/

   /*
      abstract:
         computes connected matrix (can be used for matrix inversion with complexity O(n^2)*Odet=O(N^5))
      restrictions:
         tolerance should be greater or equal zero
         matrix should be square
         both numRows and numCols should be greater or equal 2
   */
   void MatrComputeConnectedMatr(double& matr[], int numRows, int numCols, double& result[], double tolerance);

   /*
      abstract:
         performs per-element linear interpolation of two matricies
      restrictions:
         matricies should have same size
   */
   void MatrLerpMatr(double& matr1[], double& matr2[], int numRows, int numCols, double balance, double& result[]);

/*****************************************************************************************************************************/

   /*
      abstract:
         prints matrix to log
         source matrix is not modified
      restrictions:
         none
   */
   void MatrPrint(double& matr[], int numRows, int numCols);

   /*
      abstract:
         saves matrix to file
      restrictions:
         file should be opened with flags FILE_BIN|FILE_WRITE
   */
   void FileWriteMatr(double& matr[], int numRows, int numCols, int handle);

   /*
      abstract:
         saves matrix to file
      restrictions:
         file should be opened with flags FILE_BIN|FILE_READ
   */
   void FileReadMatr(double& matr[], int& numRows, int& numCols, int handle);

/*****************************************************************************************************************************/

#import

/*****************************************************************************************************************************/

