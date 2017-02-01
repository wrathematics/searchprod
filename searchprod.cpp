/*  Copyright (c) 2017, Schmidt
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <algorithm> // only for tests

using std::vector;


#define restrict __restrict__ // for g++
typedef const char *const restrict chr_r;
typedef const int *const restrict int_r;
typedef const double *const restrict dbl_r;


extern "C" void dgemm_(chr_r transa, chr_r transb, int_r m, int_r n, int_r k,
  dbl_r alpha, dbl_r A, int_r LDA, dbl_r B, int_r LDB, dbl_r beta,
  double *const restrict C, int_r LDC);


class Params
{
  public:
    Params(const int stride, const double thresh) : stride(stride), thresh(thresh) {};
    int stride;
    double thresh;
};

class Matrix
{
  public:
    Matrix(const int nr, const int nc) : nr(nr), nc(nc) {x.resize(nr*nc);};
    int nr;
    int nc;
    const int size() const{return x.size();}
    double operator [](int i) const {return x[i];}
    const double *data() const{return &x[0];}
    double *data(){return &x[0];}
    void print() const;
    
    void fill()
    {
      for (int j=0; j<nc; j++)
      {
        for (int i=0; i<nr; i++)
          x[i + nr*j] = ((double) rand() / (double) RAND_MAX);
      }
    }
    
  private:
    vector<double> x;
};

void Matrix::print() const
{
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<nc; j++)
      printf("%.3f ", x[i + nr*j]);
    
    putchar('\n');
  }
}



// A * B^T
static inline void ABT(const int Apos, const int Bpos, const int len, const Matrix &A, const Matrix &B, Matrix &C)
{
  char transa = 'N';
  char transb = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  
  
  // DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  dgemm_(&transa, &transb, &len, &len, &A.nc, &alpha, A.data() + Apos, &A.nr, B.data() + Bpos, &B.nr, &beta, C.data(), &C.nr);
}

static inline void subsearch(const int Apos, const int Bpos, const Params &p, const Matrix &A, const Matrix &B, Matrix &C, vector<double> &x)
{
  ABT(Apos, Bpos, p.stride, A, B, C);
  
  
  for (int i=0; i<C.size(); i++)
  {
    if (std::abs(C[i]) > p.thresh)
      x.push_back(C[i]);
  }
}

static void search(const Params &p, const Matrix &A, const Matrix &B, vector<double> &x)
{
  Matrix C(p.stride, p.stride);
  int B_search = B.nr;
  int Bpos = 0;
  
  
  while (B_search > 0)
  {
    int A_search = A.nr;
    int Apos = 0;
    
    while (A_search > 0)
    {
      subsearch(Apos, Bpos, p, A, B, C, x);
      
      Apos += p.stride;
      A_search -= p.stride;
    }
    
    Bpos += p.stride;
    B_search -= p.stride;
  }
}



int main()
{
  Params p(5, 1.2);
  const int m = 10;
  const int n = 3;
  Matrix A(m, n);
  Matrix B(m, n);
  
  A.fill();
  B.fill();
  
  vector<double> x;
  search(p, A, B, x);
  std::sort(x.begin(), x.end());
  for (int i=0; i<x.size(); i++)
    printf("%.3f ", x[i]);
  printf("\n\n");
  
  // compute all of them
  p.stride = m;
  vector<double> truth;
  search(p, A, B, truth);
  
  std::sort(truth.begin(), truth.end());
  for (int i=0; i<truth.size(); i++)
    printf("%.3f ", truth[i]);
  
  putchar('\n');
  
  
  return 0;
}
