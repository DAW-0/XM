
/*
   LibMatrixEx1.mq4
   Copyright © 2009 Logunov Eugene (lea)
   rootrat@mail.ru
*/

#property copyright "Copyright © 2009 Logunov Eugene (lea)"
#property link "rootrat@mail.ru"

// display the parameter window at the start
#property show_inputs

// include the library
#include <LibMatrix.mqh>

// polynomial refresh rate
extern int delay = 100;
// degree of polynomial
extern int degree = 5;
// initial distance between the control lines
extern int linesMargin = 15;
// width of chart lines
extern int linesWidth = 3;
// color of the vertical lines over the interpolation interval
extern color colVLineInt = Red;
// color of the vertical line at the end of the extrapolation interval
extern color colVLineExtr = Green;
// color of polynomial lines within the interpolation interval
extern color colInt = DimGray;
// color of polynomial lines within the extrapolation interval
extern color colExt = DarkGray;

// internal function (calculation of the X and Y moments)
double SXY(double& x[], double& y[], int numPoints, int dx, int dy) {
   double s = 0;
   for (int i = 0; i < numPoints; i++) {
      s += MathPow(x[i], dx) * MathPow(y[i], dy);
   }
   return (s / numPoints);
}

// creating polynomial regression
bool Regression(double& x[], double& y[], int numPoints, int polyDegree, double& poly[]) {
   // create system matrix
   double A[];
   int numRowsA = polyDegree + 1;
   int numColsA = polyDegree + 1;
   MatrSetSize(A, numRowsA, numColsA);
   // fill the matrix
   for (int row = 0; row < numRowsA; row++) {
      for (int col = 0; col < numColsA; col++) {
         int offset = MatrIndiciesToOffset(row, col, numRowsA, numColsA);
         A[offset] = SXY(x, y, numPoints, row + col, 0);
      }
   }
   // create a right hand side vector
   double B[];
   int numRowsB = polyDegree + 1;
   int numColsB = 1;
   MatrSetSize(B, numRowsB, numColsB);
   // fill the right hand side vector
   for (row = 0; row < numRowsB; row++) {
      offset = MatrIndiciesToOffset(row, 0, numRowsB, numColsB);
      B[offset] = SXY(x, y, numPoints, row, 1);
   }
   // solve a system of linear algebraic equations
   int numRowsX, numColsX;
   bool status =
      MatrGJBatchSolve(
         A, numRowsA, numColsA,
         B, numRowsB, numColsB,
         poly, numRowsX, numColsX,
         DEFAULT_TOLERANCE
      );
   if (!status) {
      Print("Error solving the system");
   }
   return (status);
}

// calculation of the polynomial value at a given point
double EvalPoly(double x,double& poly[],int degree) {
   double value = poly[degree];
   for (int i=degree - 1; i >= 0; i--) {
      value = (value * x) + poly[i];
   }
   return (value);
}

// flag of successful completion of initialization
bool initOk = true;
// number of lines for plotting the chart
int numLines = 0;

// script initialization
int init() {
   // trying to create vertical lines for interpolation and extrapolation interval selection
   if (!ObjectCreate("RegTime_End", OBJ_VLINE, 0, WindowTimeOnDropped(), 0)) {
      Print("Error creating the RegTime_End vertical line");
      initOk = false;
      return (0);
   }
   if (!ObjectCreate("RegTime_Start", OBJ_VLINE, 0, WindowTimeOnDropped() - linesMargin * (Period() * 60), 0)) {
      ObjectDelete("RegTime_End");
      Print("Error creating the RegTime_Start vertical line");
      initOk = false;
      return (0);
   }
   if (!ObjectCreate("RegTime_ExtEnd", OBJ_VLINE, 0, WindowTimeOnDropped() + linesMargin * (Period() * 60), 0)) {
      ObjectDelete("RegTime_End");
      ObjectDelete("RegTime_Start");
      Print("Error creating the RegTime_ExtEnd vertical line");
      initOk = false;
      return (0);
   }
   // set colors of the created lines
   ObjectSet("RegTime_End", OBJPROP_COLOR, colVLineInt);
   ObjectSet("RegTime_Start", OBJPROP_COLOR, colVLineInt);
   ObjectSet("RegTime_ExtEnd", OBJPROP_COLOR, colVLineExtr);
   // refresh the window
   WindowRedraw();
   return (0);
}

// script deinitialization
int deinit() {
   // exit if deinitialization failed
   if (!initOk) {
      return (0);
   }
   // delete the old chart lines
   for (int i = 0; i < numLines; i++) {
      ObjectDelete("RegPlot_" + i);
   }
   // delete the vertical lines
   ObjectDelete("RegTime_End");
   ObjectDelete("RegTime_Start");
   ObjectDelete("RegTime_ExtEnd");
   // refresh the window
   WindowRedraw();
   return (0);
}

// main script function (loop function)
int start() {
   // exit if deinitialization failed
   if (!initOk) {
      return (0);
   }
   // previous control line positions
   datetime timeEndPrev = -1;
   datetime timeStartPrev = -1;
   datetime timeExtEndPrev = -1;
   while (!IsStopped()) {
      // idle time for reducing the system load
      Sleep(delay);
      // check the presence of control lines in the chart (if they are not present, terminate)
      if (ObjectFind("RegTime_End") == -1) {
         Print("The RegTime_End line not found");
         break;
      }
      if (ObjectFind("RegTime_Start") == -1) {
         Print("The RegTime_Start line not found");
         break;
      }
      if (ObjectFind("RegTime_ExtEnd") == -1) {
         Print("The RegTime_ExtEnd line not found");
         break;
      }
      // read the position of control lines
      datetime timeEnd = ObjectGet("RegTime_End", OBJPROP_TIME1);
      datetime timeStart = ObjectGet("RegTime_Start", OBJPROP_TIME1);
      datetime timeExtEnd = ObjectGet("RegTime_ExtEnd", OBJPROP_TIME1);
      // check the need in updating the plot
      if ((timeEnd == timeEndPrev) && (timeStart == timeStartPrev) && (timeExtEnd == timeExtEndPrev)) {
         continue;
      }
      timeEndPrev = timeEnd;
      timeStartPrev = timeStart;
      timeExtEndPrev = timeExtEnd;
      // delete the old lines
      for (int i = 0; i < numLines; i++) {
         ObjectDelete("RegPlot_" + i);
      }
      numLines = 0;
      // determine bar shifts
      int shiftEnd = iBarShift(Symbol(), Period(), timeEnd);
      int shiftStart = iBarShift(Symbol(), Period(), timeStart);
      int shiftExtEnd = iBarShift(Symbol(), Period(), timeExtEnd);
      // calculate the number of points for interpolation and extrapolation
      int numBarsInt = shiftStart - shiftEnd + 1;
      int numBarsExt = shiftEnd - shiftExtEnd + 1;
      // prepare the data
      double x[];
      double y[];
      ArrayResize(x, numBarsInt);
      ArrayResize(y, numBarsInt);
      for (i = 0; i < numBarsInt; i++) {
         x[i] = i;
         y[i] = Close[shiftEnd + i];
      }
      // create a polynomial
      double poly[];
      if (Regression(x, y, numBarsInt, degree, poly)) {
         // plot the polynomial over the interpolation interval
         for (i = 0; i < numBarsInt - 1; i++) {
            ObjectCreate(
               "RegPlot_" + numLines,
               OBJ_TREND,
               0,
               iTime(Symbol(), Period(), shiftEnd + i),
               EvalPoly(i, poly, degree),
               iTime(Symbol(), Period(), shiftEnd + i + 1),
               EvalPoly(i + 1, poly, degree)
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_RAY,
               false
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_COLOR,
               colInt
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_WIDTH,
               linesWidth
               );
            numLines++;
         }
         // plot the polynomial over the extrapolation interval
         for (i = -1; i > -numBarsExt; i--) {
            ObjectCreate(
               "RegPlot_" + numLines,
               OBJ_TREND,
               0,
               iTime(Symbol(), Period(), shiftEnd + i),
               EvalPoly(i, poly, degree),
               iTime(Symbol(), Period(), shiftEnd + i + 1),
               EvalPoly(i + 1, poly, degree)
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_RAY,
               false
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_COLOR,
               colExt
               );
            ObjectSet(
               "RegPlot_" + numLines,
               OBJPROP_WIDTH,
               linesWidth
               );
            numLines++;
         }
      }
      // refresh the window
      WindowRedraw();
   }
   return (0);
}

