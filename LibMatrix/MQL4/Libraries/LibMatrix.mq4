
/*****************************************************************************************************************************/

/*
   LibMatrix.mq4
   Copyright © 2009 Logunov Eugene (lea)
   rootrat@mail.ru
*/

#property copyright "Copyright © 2009 Logunov Eugene (lea)"
#property link "rootrat@mail.ru"
#property library

/*****************************************************************************************************************************/

// useful mathematical constants
#define ZERO 0.0
#define ONE 1.0
#define MINUS_ONE -1.0
#define DEFAULT_TOLERANCE 0.000000001
#define RAND_MAX 32767.0

/*****************************************************************************************************************************/

/*
   abstract:
      linear interpolation between rangeLowLimit and rangeHighLimit
   restrictions:
      balance should be in range [0;1] (interpolation)
      balance is allowed to be out of range for linear extrapolation purposes olny
*/
double MathLerp(double rangeLowLimit, double rangeHighLimit, double balance) {
   // x=x0+k*(x1-x0) is equivalent of x=(1-k)*x0+k*x1, but faster
   return (rangeLowLimit + balance * (rangeHighLimit - rangeLowLimit));
}

/*
   abstract:
      returns a random number in range [rangeLowLimit;rangeHighLimit]
   restrictions:
      none
*/
double MathInRangeRandom(double rangeLowLimit, double rangeHighLimit) {
   // this code performs linear mapping of ranges in a such way:
   // [0;RAND_MAX] --> [0;1] --> [rangeLowLimit;rangeHighLimit]
   return (MathLerp(rangeLowLimit, rangeHighLimit, MathRand() / RAND_MAX));
}

/*
   abstract:
      checks two real numbers for equality with given tolerance
   restrictions:
      tolerance should be greater or equal zero
*/
bool MathDoublesEqual(double value1, double value2, double tolerance) {
   return (MathAbs(value1 - value2) <= tolerance);
}

/*****************************************************************************************************************************/

/*
   abstract:
      converts 2d matrix indicies to 1d array offset
      should be used to access matrix elements via indicies
   restrictions:
      row should be in [0;numRows-1]
      col should be in [0;numCols-1]
*/
int MatrIndiciesToOffset(int row, int col, int numRows, int numCols) {
   // matricies are packed by columns (offset is a sum of row start offset and column number)
   return ((row * numCols) + col);
}

/*
   abstract:
      copies src matrix to dst
   restrictions:
      none
*/
void MatrCopy(double& src[], double& dst[]) {
   // set dst size equal src size
   ArrayResize(dst, ArraySize(src));
   // copy all data from src to dst
   ArrayCopy(dst, src, 0, 0, WHOLE_ARRAY);
}

/*
   abstract:
      resizes 1d array so that matrix of given size could be stored in it
   restrictions:
      numRows and numCols should be positive
*/
void MatrSetSize(double& matr[], int numRows, int numCols) {
   ArrayResize(matr, numRows * numCols);
}

/*
   abstract:
      resizes given matrix without data loss
      new elements are set to zero
      old elements may be destroyed if one of dimensions is resized to a smaller value
   restrictions:
      none
*/
void MatrResize(double& matr[], int numRowsOld, int numColsOld, int numRowsNew, int numColsNew) {
   // create matrix data backup
   double tmp[];
   MatrCopy(matr, tmp);
   // resize matrix
   MatrSetSize(matr, numRowsNew, numColsNew);
   // initialize elements of resized matrix
   for (int row = 0; row < numRowsNew; row++) {
      for (int col = 0; col < numColsNew; col++) {
         // by default element is initialized with zero
         double elemValue = ZERO;
         // if non-resized matrix contained element with such indicies
         if ((row < numRowsOld) && (col < numColsOld)) {
            // copy its value from data backup
            elemValue = tmp[MatrIndiciesToOffset(row, col, numRowsOld, numColsOld)];
         }
         matr[MatrIndiciesToOffset(row, col, numRowsNew, numColsNew)] = elemValue;
      }
   }
}

/*****************************************************************************************************************************/

/*
   abstract:
      swaps two rows (row1 and row2) of matrix
   restrictions:
      row1 and row2 should be in range [0;numRows-1]
*/
void MatrSwapRows(double& matr[], int numRows, int numCols, int row1, int row2) {
   // if row is swapped with itself - there is nothing to do
   if (row1 == row2) {
      return;
   }
   // cycle through all elements of rows
   for (int col = 0; col < numCols; col++) {
      // precalculate offsets of swapped elements ([row1;col] and [row2;col])
      int offsetElem1 = MatrIndiciesToOffset(row1, col, numRows, numCols);
      int offsetElem2 = MatrIndiciesToOffset(row2, col, numRows, numCols);
      // perform swapping
      double tmp = matr[offsetElem1];
      matr[offsetElem1] = matr[offsetElem2];
      matr[offsetElem2] = tmp;
   }
}

/*
   abstract:
      swaps two columns (col1 and col2) of matrix
   restrictions:
      col1 and col2 should be in range [0;numCols-1]
*/
void MatrSwapCols(double& matr[], int numRows, int numCols, int col1, int col2) {
   // if column is swapped with itself - there is nothing to do
   if (col1 == col2) {
      return;
   }
   // cycle through all elements of columns
   for (int row = 0; row < numRows; row++) {
      // precalculate offsets of swapped elements ([row;col1] and [row;col2])
      int offsetElem1 = MatrIndiciesToOffset(row, col1, numRows, numCols);
      int offsetElem2 = MatrIndiciesToOffset(row, col2, numRows, numCols);
      // perform swapping
      double tmp = matr[offsetElem1];
      matr[offsetElem1] = matr[offsetElem2];
      matr[offsetElem2] = tmp;
   }
}

/*
   abstract:
      copy row (replace one of rows with another existing)
   restrictions:
      rowToCopy and rowToReplace should be in range [0;numRows-1]
*/
void MatrCopyRow(double& matr[], int numRows, int numCols, int rowToCopy, int rowToReplace) {
   // if row is copied over itself - there is nothing to do
   if (rowToCopy == rowToReplace) {
      return;
   }
   // cycle through all elements of rows
   for (int col = 0; col < numCols; col++) {
      // precalculate offsets of copied and replaced element
      int offsetCopy = MatrIndiciesToOffset(rowToCopy, col, numRows, numCols);
      int offsetReplace = MatrIndiciesToOffset(rowToReplace, col, numRows, numCols);
      // perform replacing
      matr[offsetReplace] = matr[offsetCopy];
   }
}

/*
   abstract:
      copy column (replace one of cols with another existing)
   restrictions:
      colToCopy and colToReplace should be in range [0;numCols-1]
*/
void MatrCopyCol(double& matr[], int numRows, int numCols, int colToCopy, int colToReplace) {
   // if column is copied over itself - there is nothing to do
   if (colToCopy == colToReplace) {
      return;
   }
   // cycle through all elements of columns
   for (int row = 0; row < numRows; row++) {
      // precalculate offsets of copied and replaced element
      int offsetCopy = MatrIndiciesToOffset(row, colToCopy, numRows, numCols);
      int offsetReplace = MatrIndiciesToOffset(row, colToReplace, numRows, numCols);
      // perform replacing
      matr[offsetReplace] = matr[offsetCopy];
   }
}

/*
   abstract:
      check whether row is zero with given tolerance
   restrictions:
      row should be in range [0;numRows-1]
*/
bool MatrRowIsZero(double& matr[], int numRows, int numCols, int row, double tolerance) {
   for (int col = 0; col < numCols; col++) {
      if (!MathDoublesEqual(matr[MatrIndiciesToOffset(row, col, numRows, numCols)], ZERO, tolerance)) {
         return (false);
      }
   }
   return (true);
}

/*
   abstract:
      check whether column is zero with given tolerance
   restrictions:
      col should be in range [0;numCols-1]
*/
bool MatrColIsZero(double& matr[], int numRows, int numCols, int col, double tolerance) {
   for (int row = 0; row < numRows; row++) {
      if (!MathDoublesEqual(matr[MatrIndiciesToOffset(row, col, numRows, numCols)], ZERO, tolerance)) {
         return (false);
      }
   }
   return (true);
}

/*****************************************************************************************************************************/

/*
   abstract:
      check whether matrix is square
   restrictions:
      none
*/
bool MatrIsSquare(int numRows, int numCols) {
   return (numRows == numCols);
}

/*
   abstract:
      check whether matrix element belongs to main diagonal
   restrictions:
      none
*/
bool MatrIsElemOnMainDiagonal(int row, int col) {
   return (row == col);
}

/*
   abstract:
      check whether two matricies are compatible for addition
   restrictions:
      none
*/
bool MatrCompatiblityCheckAdd(int numRows1, int numCols1, int numRows2, int numCols2) {
   return ((numRows1 == numRows2) && (numCols1 == numCols2));
}

/*
   abstract:
      check whether two matricies are compatible for multiplication
   restrictions:
      none
*/
bool MatrCompatiblityCheckMul(int numRows1, int numCols1, int numRows2, int numCols2) {
   return (numCols1 == numRows2);
}

/*****************************************************************************************************************************/

/*
   abstract:
      load zero matrix
   restrictions:
      none
*/
void MatrLoadZero(double& matr[], int numRows, int numCols) {
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         matr[MatrIndiciesToOffset(row, col, numRows, numCols)] = ZERO;
      }
   }
}

/*
   abstract:
      load identity matrix
   restrictions:
      matrix should be square
*/
void MatrLoadIdentity(double& matr[], int numRows, int numCols) {
   // if matrix is not square - don't modify it
   if (!MatrIsSquare(numRows, numCols)) {
      return;
   }
   // reinitialize all elements of matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         // default element value is zero
         double elemValue = ZERO;
         // if element is on main diagonal - the value is one
         if (MatrIsElemOnMainDiagonal(row, col)) {
            elemValue = ONE;
         }
         matr[MatrIndiciesToOffset(row, col, numRows, numCols)] = elemValue;
      }
   }
}

/*
   abstract:
      load matrix with random numbers in range [rangeLowLimit;rangeHighLimit]
   restrictions:
      none
*/
void MatrLoadInRangeRandom(double& matr[], int numRows, int numCols, double rangeLowLimit, double rangeHighLimit) {
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         double elemValue = MathInRangeRandom(rangeLowLimit, rangeHighLimit);
         matr[MatrIndiciesToOffset(row, col, numRows, numCols)] = elemValue;
      }
   }
}

/*****************************************************************************************************************************/

/*
   abstract:
      check whether matrix is zero using given comparison tolerance
   restrictions:
      tolerance should be greater or equal zero
*/
bool MatrIsZero(double& matr[], int numRows, int numCols, double tolerance) {
   // analyze all elements of matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         // retrive element value and compare it with zero
         double elemValue = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
         if (!MathDoublesEqual(elemValue, ZERO, tolerance)) {
            return (false);
         }
      }
   }
   // all elements are zero
   return (true);
}

/*
   abstract:
      check whether matrix is diagonal using given comparison tolerance
      matrix is assumed to be diagonal if all elements outside of main diagonal are zero
   restrictions:
      matrix should be square
      tolerance should be greater or equal zero
*/
bool MatrIsDiagonal(double& matr[], int numRows, int numCols, double tolerance) {
   // if matrix is not square - it can't be diagonal
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // analyze all elements outside of main diagonal
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         if (!MatrIsElemOnMainDiagonal(row, col)) {
            // return on any non-zero element outside of main diagonal
            double elemValue = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
            if (!MathDoublesEqual(elemValue, ZERO, tolerance)) {
               return (false);
            }
         }
      }
   }
   // matrix is diagonal (or zero - see assumption in abstract)
   return (true);
}

/*
   abstract:
      check whether matrix is identity using given comparison tolerance
   restrictions:
      matrix should be square
      tolerance should be greater or equal zero
*/
bool MatrIsIdentity(double& matr[], int numRows, int numCols, double tolerance) {
   // if matrix is not square - it can't be identity
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // analyze all elements of matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         // retrive real element value
         double elemValue = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
         // compute value which we should see in identity matrix
         double elemValueAssumed = ZERO;
         if (MatrIsElemOnMainDiagonal(row, col)) {
            elemValueAssumed = ONE;
         }
         // if real and assumed values are different - matrix is not identity
         if (!MathDoublesEqual(elemValue, elemValueAssumed, tolerance)) {
            return (false);
         }
      }
   }
   // matrix is identity
   return (true);
}

/*
   abstract:
      check whether matrix is symmetric (a[i][j]=a[j][i]) using given comparison tolerance
   restrictions:
      matrix should be square
      tolerance should be greater or equal zero
*/
bool MatrIsSymmetric(double& matr[], int numRows, int numCols, double tolerance) {
   // if matrix is not square - it can't be symmetric
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // retrive and compare values of elements which can be equal in symmetric matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = row + 1; col < numCols; col++) {
         double elemValue1 = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
         double elemValue2 = matr[MatrIndiciesToOffset(col, row, numRows, numCols)];
         if (!MathDoublesEqual(elemValue1, elemValue2, tolerance)) {
            // if any pair of values is not equal - matrix is not symmetric
            return (false);
         }
      }
   }
   // matrix is symmetric
   return (true);
}

/*
   abstract:
      check whether matrix is antisymmetric (a[i][j]=-a[j][i]; a[k][k]=0) using given comparison tolerance
   restrictions:
      matrix should be square
      tolerance should be greater or equal zero
*/
bool MatrIsAntisymmetric(double& matr[], int numRows, int numCols, double tolerance) {
   // if matrix is not square - it can't be cosymmetric
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // analyze one triangle of indicies
   for (int row = 0; row < numRows; row++) {
      for (int col = row; col < numCols; col++) {
         // if indicies doesn't belong main diagonal
         if (!MatrIsElemOnMainDiagonal(row, col)) {
            // check a pair of elements
            double elemValue1 = MINUS_ONE * matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
            double elemValue2 = matr[MatrIndiciesToOffset(col, row, numRows, numCols)];
            if (!MathDoublesEqual(elemValue1, elemValue2, tolerance)) {
               return (false);
            }
         }
         else {
            // indicies belong main diagonal - element should be zero
            double elemValue = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
            if (!MathDoublesEqual(elemValue, ZERO, tolerance)) {
               return (false);
            }
         }
      }
   }
   // matrix is cosymmetric
   return (true);
}

/*
   abstract:
      check whether matricies are equal
   restrictions:
      matricies should have same size
      tolerance should be greater or equal zero
*/
bool MatrEqual(double& matr1[], double& matr2[], int numRows, int numCols, double tolerance) {
   // analyze all elements
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         // retrive and compare elements with same indicies
         double elemMatr1 = matr1[MatrIndiciesToOffset(row, col, numRows, numCols)];
         double elemMatr2 = matr2[MatrIndiciesToOffset(row, col, numRows, numCols)];
         // if elements differ - matricies are not equal
         if (!MathDoublesEqual(elemMatr1, elemMatr2, tolerance)) {
            return (false);
         }
      }
   }
   // matrices are equal
   return (true);
}

/*****************************************************************************************************************************/

/*
   abstract:
      adds scalar value to each element of matrix
   restrictions:
      none
*/
void MatrAddScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]) {
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         tmpResult[offset] = matr[offset] + scalar;
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      subtracts scalar value from each element of matrix
   restrictions:
      none
*/
void MatrSubScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]) {
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         tmpResult[offset] = matr[offset] - scalar;
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      multiplies each element of matrix by scalar value
   restrictions:
      none
*/
void MatrMulByScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]) {
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         tmpResult[offset] = matr[offset] * scalar;
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      divides each element of matrix by scalar value
   restrictions:
      scalar should be non-zero
*/
void MatrDivByScalar(double& matr[], int numRows, int numCols, double scalar, double& result[]) {
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   double oneOverScalar = ONE / scalar;
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         tmpResult[offset] = matr[offset] * oneOverScalar;
      }
   }
   MatrCopy(tmpResult, result);
}

/*****************************************************************************************************************************/

/*
   abstract:
      adds second matrix to first
   restrictions:
      matricies should have same size
*/
void MatrAddMatr(double& matr1[], double& matr2[], int numRows, int numCols, double& result[]) {
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         tmpResult[offset] = matr1[offset] + matr2[offset];
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      subtracts second matrix from first
   restrictions:
      matricies should have same size
*/
void MatrSubMatr(double& matr1[], double& matr2[], int numRows, int numCols, double& result[]) {
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(row,col,numRows,numCols);
         tmpResult[offset]= matr1[offset] - matr2[offset];
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      multiplies first matrix by second
   restrictions:
      matricies should be compatible for multiplication
*/
void MatrMulMatr(
      double& matr1[], int numRows1, int numCols1,
      double& matr2[], int numRows2, int numCols2,
      double& result[], int& numRowsRes, int& numColsRes) {
   // calculate size of result matrix
   numRowsRes = numRows1;
   numColsRes = numCols2;
   // calculate number of multiplications/additions in each inner dot product
   int numElemsOfDotProduct = numCols1;
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRowsRes, numColsRes);
   // compute all elements of result matrix
   for (int row = 0; row < numRowsRes; row++) {
      for (int col = 0; col < numColsRes; col++) {
         // offset of element in result matrix
         int productElemOffset = MatrIndiciesToOffset(row, col, numRowsRes, numColsRes);
         // compute inner dot product (row from first matrix multiplied by column of second)
         tmpResult[productElemOffset] = ZERO;
         for (int elem = 0; elem < numElemsOfDotProduct; elem++) {
            int srcElem1Offset = MatrIndiciesToOffset(row, elem, numRows1, numCols1);
            int srcElem2Offset = MatrIndiciesToOffset(elem, col, numRows2, numCols2);
            tmpResult[productElemOffset] += matr1[srcElem1Offset] * matr2[srcElem2Offset];
         }
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      computes trace (sum of diagonal elements) of given matrix
      source matrix is not modified
   restrictions:
      matrix should be square
*/
double MatrTrace(double& matr[], int numRows, int numCols) {
   if (!MatrIsSquare(numRows, numCols)) {
      return (EMPTY_VALUE);
   }
   double trace = ZERO;
   for (int elem = 0; elem < numRows; elem++) {
      trace += matr[MatrIndiciesToOffset(elem, elem, numRows, numCols)];
   }
   return (trace);
}

/*
   abstract:
      performs matrix transpose
   restrictions:
      none
*/
void MatrTranspose(double& matr[], int numRows, int numCols, double& result[], int& numRowsRes, int& numColsRes) {
   // compute size of transposed matrix
   numRowsRes = numCols;
   numColsRes = numRows;
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRowsRes, numColsRes);
   // write elements back from temporary matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         tmpResult[MatrIndiciesToOffset(col, row, numRowsRes, numColsRes)] = matr[MatrIndiciesToOffset(row, col, numRows, numCols)];
      }
   }
   MatrCopy(tmpResult, result);
}

/*****************************************************************************************************************************/

/*
   abstract:
      performs Gaussian triangulation of matrix
      returns size of subset of linear-independent rows (matrix rank)
   restrictions:
      tolerance should be greater or equal zero
*/
int MatrGaussianElimination(double& matr[], int numRows, int numCols, double& result[], double tolerance) {
   // create temporary matrix for result
   double tmpResult[];
   MatrCopy(matr, tmpResult);
   // current processed element position
   int currentRow = 0;
   int currentCol = 0;
   // process whole matrix
   while ((currentRow < numRows) && (currentCol < numCols)) {
      // perform pivoting (search for the element of column with maximal absolute value)
      double bestAbs = ZERO;
      int bestPivotingRow = currentRow;
      for (int searchedRow = currentRow; searchedRow < numRows; searchedRow++) {
         int offset = MatrIndiciesToOffset(searchedRow, currentCol, numRows, numCols);
         double currentAbs = MathAbs(tmpResult[offset]);
         if (currentAbs > bestAbs) {
            bestAbs = currentAbs;
            bestPivotingRow = searchedRow;
         }
      }
      // ensure that pivoting element is non-zero
      if (bestAbs <= tolerance) {
         // pivoting element is zero - fill the column with exact zeros
         for (int zeroFilledRow = currentRow; zeroFilledRow < numRows; zeroFilledRow++) {
            tmpResult[MatrIndiciesToOffset(zeroFilledRow, currentCol, numRows, numCols)] = ZERO;
         }
         // continue with next column
         currentCol++;
         continue;
      }
      // swap rows if needed
      if (currentRow != bestPivotingRow) {
         for (int swappedCol = currentCol; swappedCol < numCols; swappedCol++) {
            int offset1 = MatrIndiciesToOffset(currentRow, swappedCol, numRows, numCols);
            int offset2 = MatrIndiciesToOffset(bestPivotingRow, swappedCol, numRows, numCols);
            double tmp = tmpResult[offset1];
            tmpResult[offset1] = tmpResult[offset2];
            tmpResult[offset2] = MINUS_ONE * tmp;
         }
      }
      // perform elimination
      for (int eliminatedRow = currentRow + 1; eliminatedRow < numRows; eliminatedRow++) {
         // precompute coeffitient for elimination
         int numeratorOffset = MatrIndiciesToOffset(eliminatedRow, currentCol, numRows, numCols);
         int denominatorOffset = MatrIndiciesToOffset(currentRow, currentCol, numRows, numCols);
         double coeff = -tmpResult[numeratorOffset] / tmpResult[denominatorOffset];
         // add current row multiplied by precomputed coeffitient to eliminated row
         tmpResult[MatrIndiciesToOffset(eliminatedRow, currentCol, numRows, numCols)] = ZERO;
         for (int colOfEliminatedRow = currentCol + 1; colOfEliminatedRow < numCols; colOfEliminatedRow++) {
            int updatedElementOffset = MatrIndiciesToOffset(eliminatedRow, colOfEliminatedRow, numRows, numCols);
            int addedElementOffset = MatrIndiciesToOffset(currentRow, colOfEliminatedRow, numRows, numCols);
            tmpResult[updatedElementOffset] += coeff * tmpResult[addedElementOffset];
         }
      }
      // continue with next row/column
      currentRow++;
      currentCol++;
   }
   MatrCopy(tmpResult, result);
   // return number of non-zero rows (matrix rank)
   return (currentRow);
}

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
      ) {
   // store matricies sizes in more convenient way
   int numEquations = numRows;
   int numVars = numCols;
   int numRHS = numColsRHS;
   int numVarsAndRHS = numVars + numRHS;
   // allocate temporary matrix for the result
   numRowsRoots = numVars;
   numColsRoots = numRHS;
   double tmpRoots[];
   MatrSetSize(tmpRoots, numVars, numRHS);
   // create extended matrix system
   double extMatr[];
   MatrCopy(matr, extMatr);
   MatrResize(extMatr, numRows, numCols, numEquations, numVarsAndRHS);
   for (int appendedCol = 0; appendedCol < numRHS; appendedCol++) {
      for (int row = 0; row < numEquations; row++) {
         int offsetFrom = MatrIndiciesToOffset(row, appendedCol, numRowsRHS, numColsRHS);
         int offsetTo = MatrIndiciesToOffset(row, numVars + appendedCol, numEquations, numVarsAndRHS);
         extMatr[offsetTo] = rhs[offsetFrom];
      }
   }
   // perform Gaussian elimination with pivoting (applied only to matrix
   // system part, not to the whole extended system matrix).
   // current processed element position
   int currentRow = 0;
   int currentCol = 0;
   while ((currentRow < numEquations) && (currentCol < numVars)) {
      // perform pivoting (search for the element of column with maximal absolute value)
      double bestAbs = ZERO;
      int bestPivotingRow = currentRow;
      for (int searchedRow = currentRow; searchedRow < numEquations; searchedRow++) {
         int offset = MatrIndiciesToOffset(searchedRow, currentCol, numEquations, numVarsAndRHS);
         double currentAbs = MathAbs(extMatr[offset]);
         if (currentAbs > bestAbs) {
            bestAbs = currentAbs;
            bestPivotingRow = searchedRow;
         }
      }
      // ensure that pivoting element is non-zero
      if (bestAbs <= tolerance) {
         // pivoting element is zero - fill the column with exact zeros
         for (int zeroFilledRow = currentRow; zeroFilledRow < numEquations; zeroFilledRow++) {
            extMatr[MatrIndiciesToOffset(zeroFilledRow, currentCol, numEquations, numVarsAndRHS)] = ZERO;
         }
         // continue with next column
         currentCol++;
         continue;
      }
      // swap rows if needed
      if (currentRow != bestPivotingRow) {
         for (int swappedCol = currentCol; swappedCol < numVarsAndRHS; swappedCol++) {
            int offset1 = MatrIndiciesToOffset(currentRow, swappedCol, numEquations, numVarsAndRHS);
            int offset2 = MatrIndiciesToOffset(bestPivotingRow, swappedCol, numEquations, numVarsAndRHS);
            double tmp = extMatr[offset1];
            extMatr[offset1] = extMatr[offset2];
            extMatr[offset2] = MINUS_ONE * tmp;
         }
      }
      // perform elimination
      for (int eliminatedRow = currentRow + 1; eliminatedRow < numEquations; eliminatedRow++) {
         // precompute coeffitient for elimination
         int numeratorOffset = MatrIndiciesToOffset(eliminatedRow, currentCol, numEquations, numVarsAndRHS);
         int denominatorOffset = MatrIndiciesToOffset(currentRow, currentCol, numEquations, numVarsAndRHS);
         double coeff = -extMatr[numeratorOffset] / extMatr[denominatorOffset];
         // add current row multiplied by precomputed coeffitient to eliminated row
         extMatr[MatrIndiciesToOffset(eliminatedRow, currentCol, numEquations, numVarsAndRHS)] = ZERO;
         for (int colOfEliminatedRow = currentCol + 1; colOfEliminatedRow < numVarsAndRHS; colOfEliminatedRow++) {
            int updatedElementOffset = MatrIndiciesToOffset(eliminatedRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            int addedElementOffset = MatrIndiciesToOffset(currentRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            extMatr[updatedElementOffset] += coeff * extMatr[addedElementOffset];
         }
      }
      // continue with next row/column
      currentRow++;
      currentCol++;
   }
   // remember matrix rank for additional analysis
   int numNonZeroRows = currentRow;
   // ensure that all the rows on bottom are zero (especially in the part of right-hand sides)
   for (int rowOnBottom = numNonZeroRows; rowOnBottom < numEquations; rowOnBottom++) {
      if (!MatrRowIsZero(extMatr, numEquations, numVarsAndRHS, rowOnBottom, tolerance)) {
         return (false);
      }
   }
   // perform Gauss-Jordan reduction
   currentRow = numNonZeroRows - 1;
   currentCol = numVars - 1;
   while ((currentRow > 0) && (currentCol >= 0)) {
      // if current variable is unused in system of simultaneous equations - skip it
      if (MatrColIsZero(matr, numEquations, numVarsAndRHS, currentCol, tolerance)) {
         currentCol--;
         continue;
      }
      // search for first non-zero coeffitient for current variable
      while (MathDoublesEqual(extMatr[MatrIndiciesToOffset(currentRow, currentCol, numEquations, numVarsAndRHS)], ZERO, tolerance)) {
         currentRow--;
      }
      // eliminate this this variable from the rest of equations
      for (eliminatedRow = currentRow - 1; eliminatedRow >= 0; eliminatedRow--) {
         // precompute coeffitient for elimination
         numeratorOffset = MatrIndiciesToOffset(eliminatedRow, currentCol, numEquations, numVarsAndRHS);
         denominatorOffset = MatrIndiciesToOffset(currentRow, currentCol, numEquations, numVarsAndRHS);
         coeff = -extMatr[numeratorOffset] / extMatr[denominatorOffset];
         // add currentRow row multiplied by coeff to eliminatedRow row
         extMatr[MatrIndiciesToOffset(eliminatedRow, currentCol, numEquations, numVarsAndRHS)] = ZERO;
         for (colOfEliminatedRow = 0; colOfEliminatedRow < currentCol; colOfEliminatedRow++) {
            updatedElementOffset = MatrIndiciesToOffset(eliminatedRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            addedElementOffset = MatrIndiciesToOffset(currentRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            extMatr[updatedElementOffset] += coeff * extMatr[addedElementOffset];
         }
         for (colOfEliminatedRow = currentCol + 1; colOfEliminatedRow < numVarsAndRHS; colOfEliminatedRow++) {
            updatedElementOffset = MatrIndiciesToOffset(eliminatedRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            addedElementOffset = MatrIndiciesToOffset(currentRow, colOfEliminatedRow, numEquations, numVarsAndRHS);
            extMatr[updatedElementOffset] += coeff * extMatr[addedElementOffset];
         }
      }
      // move to the next variable for elimination
      currentCol--;
      // prevent continious elimination from the same equation - move to the next equation
      currentRow--;
   }
   // if rank is less then number of equations - there is no unique solution
   if (numNonZeroRows < numEquations) {
      return (false);
   }
   // calculate roots for each right-hand side
   // analyze all non-zero rows
   for (currentRow = 0; currentRow < numNonZeroRows; currentRow++) {
      // search for first non-zero element in coeffitients part of extended matrix
      int computedVariable = -1;
      for (int checkedCol = 0; checkedCol < numVars; checkedCol++) {
         int checkedElementOffset = MatrIndiciesToOffset(currentRow, checkedCol, numEquations, numVarsAndRHS);
         if (MathAbs(extMatr[checkedElementOffset]) > tolerance) {
            if (computedVariable != -1) {
               // There should be only one non-zero element in each column. If there is more then
               // one non-zero element - the solution may exist, but it seems, that algorithm is
               // incorrect. The assertion should always be fine by the way (no crashes during testing).
               Print("LibMatrix: MatrGJBatchSolve() crash.");
               Print("LibMatrix: Assertion (computedVariable==-1) failure.");
               Print("LibMatrix: Contact library developer.");
            }
            // if we found first non-zero element - remember it
            computedVariable = checkedCol;
         }
      }
      // check whether selected variable values can be computed
      if (computedVariable == -1) {
         // equation is destroyed during elimination, but right-hand side is non-zero - no solution
         if (!MatrRowIsZero(extMatr, numEquations, numVarsAndRHS, currentRow, tolerance)) {
            return (false);
         }
      }
      else {
         // calculate selected variable value for each right-hand side
         int offsetCoeffOfReducedMatr = MatrIndiciesToOffset(currentRow, computedVariable, numEquations, numVarsAndRHS);
         double coeffAtVar = ONE / extMatr[offsetCoeffOfReducedMatr];
         for (int indexRHS = 0; indexRHS < numRHS; indexRHS++) {
            int offsetVar = MatrIndiciesToOffset(computedVariable, indexRHS, numVars, numRHS);
            int offsetCoeffRHS = MatrIndiciesToOffset(currentRow, numVars + indexRHS, numEquations, numVarsAndRHS);
            tmpRoots[offsetVar] = extMatr[offsetCoeffRHS] * coeffAtVar;
         }
      }
   }
   // batch solution found
   MatrCopy(tmpRoots, roots);
   return (true);
}

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
      double& result[], int& numRowsRes, int& numColsRes) {
   // create matrix for minor elements
   numRowsRes = numRows - 1;
   numColsRes = numCols - 1;
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRowsRes, numColsRes);
   // extract minor skipping 1 row and 1 column
   int rowOfMinor = 0;
   for (int row = 0; row < numRows; row++) {
      if (row != rowToExclude) {
         int colOfMinor = 0;
         for (int col = 0; col < numCols; col++) {
            if (col != colToExclude) {
               int offsetMatr = MatrIndiciesToOffset(row, col, numRows, numCols);
               int offsetMinor = MatrIndiciesToOffset(rowOfMinor, colOfMinor, numRowsRes, numColsRes);
               tmpResult[offsetMinor] = matr[offsetMatr];
               colOfMinor++;
            }
         }
         rowOfMinor++;
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      computes algebraic complement of given element (a determinant of its minor with sign)
      Gaussian elimination with given pivoting tolerance is used to compute determinant
   restrictions:
      tolerance should be greater or equal zero
      matrix should be square
      both numRows and numCols should be greater or equal 2
*/
double MatrAlgebraicComplement(double& matr[], int numRows, int numCols, int elemRow, int elemCol, double tolerance) {
   // if matrix is not square - it's impossible to compute determinant
   if (!MatrIsSquare(numRows, numCols)) {
      return (EMPTY_VALUE);
   }
   // extract minor from matrix
   double minor[];
   int minorRows, minorCols;
   MatrMinor(matr, numRows, numCols, elemRow, elemCol, minor, minorRows, minorCols);
   // compute minor's determinant
   double det = MatrDet(minor, minorRows, minorRows, tolerance);
   // return determinant with sign
   return (det * MathPow(MINUS_ONE, elemRow + elemCol));
}

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
bool MatrInvertUsingMinors(double& matr[], int numRows, int numCols, double& result[], double tolerance) {
   // if matrix is not square - it couldn't be inverted
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // if determinant of source matrix is zero - it couldn't be inverted
   double det = MatrDet(matr, numRows, numCols, tolerance);
   if (MathDoublesEqual(det, ZERO, tolerance)) {
      return (false);
   }
   // compute connected matrix and divide it by determinant of source matrix - the result is inverted matrix
   MatrComputeConnectedMatr(matr, numRows, numCols, result, tolerance);
   MatrDivByScalar(result, numRows, numCols, det, result);
   return (true);
}

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
bool MatrInvertUsingGJ(double& matr[], int numRows, int numCols, double& result[], double tolerance) {
   // if matrix is not square - it couldn't be inverted
   if (!MatrIsSquare(numRows, numCols)) {
      return (false);
   }
   // create temporary matrix for result
   double tmpResult[];
   // create right hand side matrix
   double identity[];
   MatrSetSize(identity, numRows, numCols);
   MatrLoadIdentity(identity, numRows, numCols);
   // apply Gauss-Jordan reduction to find inversion
   bool gjResult = MatrGJBatchSolve(matr, numRows, numCols, identity, numRows, numCols, tmpResult, numRows, numCols, tolerance);
   MatrCopy(tmpResult, result);
   return (gjResult);
}

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
double MatrDet(double& matr[], int numRows, int numCols, double tolerance) {
   // if matrix is not square - determinant can't be computed
   if (!MatrIsSquare(numRows, numCols)) {
      return (EMPTY_VALUE);
   }
   // perform triangulation on source matrix
   double matrTriangulated[];
   MatrGaussianElimination(matr, numRows, numCols, matrTriangulated, tolerance);
   // compute determinant using triangulated matrix
   double det = ONE;
   for (int elem = 0; elem < numRows; elem++) {
      det *= matrTriangulated[MatrIndiciesToOffset(elem, elem, numRows, numCols)];
   }
   return (det);
}

/*
   abstract:
      computes determinant of triangulated matrix
      source matrix is not modified
   restrictions:
      matrix should be square
*/
double MatrDetTriang(double& matr[], int numRows, int numCols) {
   // if matrix is not square - determinant can't be computed
   if (!MatrIsSquare(numRows, numCols)) {
      return (EMPTY_VALUE);
   }
   // compute determinant
   double det = ONE;
   for (int elem = 0; elem < numRows; elem++) {
      det *= matr[MatrIndiciesToOffset(elem, elem, numRows, numCols)];
   }
   return (det);
}

/*****************************************************************************************************************************/

/*
   abstract:
      computes matrix rank using Gaussian elimination with given pivoting tolerance
      source matrix is not modified
   restrictions:
      tolerance should be greater or equal zero
*/
double MatrRank(double& matr[], int numRows, int numCols, double tolerance) {
   // temporary matrix
   double matrTriangulated[];
   // apply Gaussian elimination method and return number of non-zero rows
   return (MatrGaussianElimination(matr, numRows, numCols, matrTriangulated, tolerance));
}

/*
   abstract:
      computes rank of triangulated matrix
      source matrix is not modified
      tolerance is used when searching for zero row
   restrictions:
      tolerance should be greater or equal zero
*/
double MatrRankTriang(double& matr[], int numRows, int numCols, double tolerance) {
   // matrix is triangulated - we just need to find first zero row
   int row = 0;
   while (row < numRows) {
      if (MatrRowIsZero(matr, numRows, numCols, row, tolerance)) {
         break;
      }
      row++;
   }
   return (row);
}

/*****************************************************************************************************************************/

/*
   abstract:
      computes connected matrix (can be used for matrix inversion with complexity O(n^2)*Odet=O(N^5))
   restrictions:
      tolerance should be greater or equal zero
      matrix should be square
      both numRows and numCols should be greater or equal 2
*/
void MatrComputeConnectedMatr(double& matr[], int numRows, int numCols, double& result[], double tolerance) {
   // if matrix is not square - it's impossible to compute determinant
   if (!MatrIsSquare(numRows, numCols)) {
      return;
   }
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   // compute algebraic complement of each element and transpose matrix (at one time)
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         int offset = MatrIndiciesToOffset(col, row, numRows, numCols);
         tmpResult[offset] = MatrAlgebraicComplement(matr, numRows, numCols, row, col, tolerance);
      }
   }
   MatrCopy(tmpResult, result);
}

/*
   abstract:
      performs per-element linear interpolation of two matricies
   restrictions:
      matricies should have same size
*/
void MatrLerpMatr(double& matr1[], double& matr2[], int numRows, int numCols, double balance, double& result[]) {
   // create temporary matrix for result
   double tmpResult[];
   MatrSetSize(tmpResult, numRows, numCols);
   // for each element of result matrix
   for (int row = 0; row < numRows; row++) {
      for (int col = 0; col < numCols; col++) {
         // precompute element offset
         int offset = MatrIndiciesToOffset(row, col, numRows, numCols);
         // perform interpolation
         tmpResult[offset] = MathLerp(matr1[offset], matr2[offset], balance);
      }
   }
   MatrCopy(tmpResult, result);
}

/*****************************************************************************************************************************/

/*
   abstract:
      prints matrix to log
      source matrix is not modified
   restrictions:
      none
*/
void MatrPrint(double& matr[], int numRows, int numCols) {
   for (int row = numRows - 1; row >= 0; row--) {
      string line="";
      for (int col = 0; col < numCols; col++) {
         line = StringConcatenate(line, " ", DoubleToStr(matr[MatrIndiciesToOffset(row, col, numRows, numCols)], 16));
      }
      Print(line);
   }
}

/*
   abstract:
      saves matrix to file
   restrictions:
      file should be opened with flags FILE_BIN|FILE_WRITE
*/
void FileWriteMatr(double& matr[], int numRows, int numCols, int handle) {
   FileWriteInteger(handle, numRows, LONG_VALUE);
   FileWriteInteger(handle, numCols, LONG_VALUE);
   FileWriteArray(handle, matr, 0, numRows * numCols);
}

/*
   abstract:
      saves matrix to file
   restrictions:
      file should be opened with flags FILE_BIN|FILE_READ
*/
void FileReadMatr(double& matr[], int& numRows, int& numCols, int handle) {
   numRows = FileReadInteger(handle, LONG_VALUE);
   numCols = FileReadInteger(handle, LONG_VALUE);
   MatrSetSize(matr, numRows, numCols);
   FileReadArray(handle, matr, 0, numRows * numCols);
}

/*****************************************************************************************************************************/


