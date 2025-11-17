#include <math.h>

// Multiply C = A * B, where A is a m x p matrix, and B is a p x n matrix.
// All matrices A, B, C must be pre-allocated (say, using malloc or similar).
// The memory storage for C must *not* overlap in memory with either A or B. 
// That is, you **cannot** do C = A * C, or C = C * B. However, A and B can overlap, and so C = A * A is fine, as long as the memory buffer for A is not overlaping in memory with that of C.
// Very important: All matrices are stored in **column-major** format.
// Example. Suppose 
//      [ 1 8 2 ]
//  A = [ 3 5 7 ]
//      [ 0 2 4 ]
//  Then, the storage in memory is
//   1, 3, 0, 8, 5, 2, 2, 7, 4. 
void MultiplyMatrices(int m, int p, int n, const double * A, const double * B, double * C)
{
  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
    {
      double entry = 0.0;
      for(int k=0; k<p; k++)
        entry += A[k * m + i] * B[j * p + k];
      C[m * j + i] = entry;
    }
  }
}

double vecLen(const double v[3]) {
  return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void normalize(double v[3]) {
  double L = vecLen(v);
  if (L > 1e-8f) {
    v[0] /= L;
    v[1] /= L;
    v[2] /= L;
    }
}

void scale(double v[3], const double s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

void deep_copy(double dest[3], const double origin[3]) {
  dest[0] = origin[0];
  dest[1] = origin[1];
  dest[2] = origin[2];
}

double dot(const double a[3], const double b[3]) {
  return((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]));
}

void cross(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];
}

void subtract(const double a[3], const double b[3], double out[3]) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

void add(const double a[3], const double b[3], double out[3]) {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
}

void add(const double a[3], const double scalar, double out[3]) {
  out[0] = a[0] + scalar;
  out[1] = a[1] + scalar;
  out[2] = a[2] + scalar;
}

void divide(const double a[3], const double scalar, double out[3]) {
  out[0] = a[0] / scalar;
  out[1] = a[1] / scalar;
  out[2] = a[2] / scalar;
}

double dist(const double a[3], const double b[3]) {
    double res[3];
    subtract(b, a, res);
    return std::sqrt(std::pow(res[0], 2) + std::pow(res[1], 2) + std::pow(res[2], 2));
}

double magnitude(const double a[3]){
  return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}