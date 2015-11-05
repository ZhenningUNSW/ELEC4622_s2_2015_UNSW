/*****************************************************************************/
// File: dft.cpp
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#define _USE_MATH_DEFINES  // Makes `M_PI' available as a constant.

#include <stdlib.h>
#include <cmath>
#include "dft.h"
#include <complex>

/*****************************************************************************/
/*                            my_direct_dft::init                            */
/*****************************************************************************/

void my_direct_dft::init(int N, bool is_forward)
{
  cleanup(); // Delete any pre-existing buffers.
  this->N = N;
  real_buf = new float[N];
  imag_buf = new float[N];
  real_trig = new double[N];
  imag_trig = new double[N];
  for (int n=0; n < N; n++)
    {
      real_trig[n] = cos(n * 2.0 * M_PI / N);
      imag_trig[n] = sin(n * 2.0 * M_PI / N);
      if (is_forward)
        imag_trig[n] = -imag_trig[n];
    }
}

/*****************************************************************************/
/*                   my_direct_dft::perform_transform                        */
/*****************************************************************************/

void my_direct_dft::perform_transform(float *real, float *imag, int stride)
{
  // First copy the input values to the real and imaginary temporary buffers.
  int n, k;
  float *rp, *ip;
  for (rp=real, ip=imag, n=0; n < N; n++, rp+=stride, ip+=stride)
    { real_buf[n] = *rp;  imag_buf[n] = *ip; }

  // Now compute each output coefficient in turn
  for (rp=real, ip=imag, k=0; k < N; k++, rp+=stride, ip+=stride)
    {
      int index = 0; // This holds n*k mod N; it indexes the trig tables
      double real_sum=0.0, imag_sum=0.0;
      for (n=0; n < N; n++, index+=k)
        {
          if (index >= N)
            index -= N;
          real_sum += real_buf[n]*real_trig[index]
                    - imag_buf[n]*imag_trig[index];
          imag_sum += real_buf[n]*imag_trig[index]
                    + imag_buf[n]*real_trig[index];
        }
      *rp = (float) real_sum;
      *ip = (float) imag_sum;
    }
}

void my_fft::init(int N, bool is_forward)
{
	cleanup(); // Delete any pre-existing buffers.
/*	this->N = N;
	buf = new _complex[N];
	trig = new _complex[N];
	for (int n = 0; n < N; n++)	
	{
		trig[n].x = cos(n * 2.0 * M_PI / N);
		trig[n].y = sin(n * 2.0 * M_PI / N);
		if (is_forward)
			trig[n].y = -trig[n].y;
	}
*/
}

_complex* my_fft::perform_transform(_complex *x, int N)
{
	_complex * X = new _complex[N];
	_complex  *D, *E;
	int k;

	if (N == 1) {
		X[0] = x[0];
		return X;
	}

	_complex *e = new _complex[N / 2];
	_complex *d = new _complex[N / 2];
	for (k = 0; k < N / 2; k++) {
		e[k] = x[2 * k];
		d[k] = x[2 * k + 1];
	}

	E = this->perform_transform(e, N / 2);
	D = this->perform_transform(d, N / 2);

	for (k = 0; k < N / 2; k++) {
		_complex Wn;
		_complex d;
		Wn.x = cos(-2 * M_PI * k / N);
		Wn.y = sin(-2 * M_PI * k / N);
		d.x = D[k].x * Wn.x - Wn.y * D[k].y;
		d.y = D[k].y * Wn.x + D[k].x * Wn.y;
		/* Multiply entries of D by the twiddle factors e^(-2*pi*i/N * k) */
		X[k].x = E[k].x + d.x;
		X[k].y = E[k].y + d.y;
		X[k + N / 2].x = E[k].x - d.x;
		X[k + N / 2].y = E[k].y - d.y;
	}

	delete[] D;
	delete[] E;
	return X;
}