/*****************************************************************************/
// File: vector_filter.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/
#include <smmintrin.h>//SSE4.1(include tmmintrin.h)
#include <emmintrin.h> // Include SSE2 processor intrinsic functions
#include <stdlib.h>
#include "aligned_image_comps.h"

/*****************************************************************************/
/*                     my_aligned_image_comp::vector_filter                  */
/*****************************************************************************/

void my_aligned_image_comp::vector_filter(my_aligned_image_comp *in)
{
#define FILTER_EXTENT 4
#define FILTER_TAPS (2*FILTER_EXTENT+1)

  // Create the vertical filter PSF as a local array on the stack.
  __m128 filter_buf[FILTER_TAPS];
  __m128 *mirror_psf = filter_buf+FILTER_EXTENT;

  float hk = 1.0F / FILTER_TAPS;
          // `mirror_psf' points to the central tap in the filter
  for (int t=-FILTER_EXTENT; t <= FILTER_EXTENT; t++)
    mirror_psf[t] = _mm_set1_ps(hk);

  // Check for consistent dimensions
  assert(in->border >= FILTER_EXTENT);
  assert((this->height <= in->height) && (this->width <= in->width));
  assert(((stride & 3) == 0) && ((in->stride & 3) == 0));
  int vec_stride_in = in->stride / 4;
  int vec_stride_out = this->stride / 4;
  int vec_width_out = (this->width+3)/4; // Big enough to cover the width

  // Do the filtering
  __m128 *line_out = (__m128 *) buf;
  __m128 *line_in = (__m128 *)(in->buf);
//  const int imm8 = 0xFFFFFFFF;
  for (int r=0; r < height; r++,
       line_out+=vec_stride_out, line_in+=vec_stride_in)
    for (int c=0; c < vec_width_out; c++)
      {
        __m128 *ip = (line_in+c) - vec_stride_in*FILTER_EXTENT;
        __m128 sum = _mm_setzero_ps();
		for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++, ip += vec_stride_in)
			sum = _mm_add_ps(sum,_mm_mul_ps(mirror_psf[y],*ip));
			//sum = _mm_dp_ps(mirror_psf[y], *ip, imm8);
        line_out[c] = sum;
      }
}
