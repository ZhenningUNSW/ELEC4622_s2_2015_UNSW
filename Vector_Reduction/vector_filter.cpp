/*****************************************************************************/
// File: vector_filter.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <emmintrin.h> // Include SSE2 processor intrinsic functions
#include <smmintrin.h>
#include <stdlib.h>
#include "aligned_image_comps.h"

/*****************************************************************************/
/*                     my_aligned_image_comp::vector_filter                  */
/*****************************************************************************/

void my_aligned_image_comp::vector_filter(my_aligned_image_comp *in, int filter_length)
{
#define FILTER_EXTENT 14
#define FILTER_TAPS (2*FILTER_EXTENT+1)
	int shifted_filter_cof[FILTER_TAPS] = { 0, -3, 0, 35, 41, -69, -171, 0, 355, 305,
											-425, -978, 0, 2576, 4890, 2576, 0, -978,
											-425, 305, 355, 0, -171, -69, 41, 35, 0, -3 };
	int filter_cof[FILTER_TAPS] = { -1, -5, 15, 49, 0, -139, -128, 185, 419, 0, -817, -734, 1170,
									3922, 5243, 3922, 1170, -734, -817, 0, 419, 185, -128, -139, 0,
									49, 15, -5, -1 };
	// Create the vertical filter PSF as a local array on the stack.
	__m128i filter_buf[FILTER_TAPS];
	__m128i *mirror_psf = filter_buf + FILTER_EXTENT;
	__m128i filter_buf_2[FILTER_TAPS];
	__m128i *mirror_psf_2 = filter_buf_2 + FILTER_EXTENT;

	// `mirror_psf' points to the central tap in the filter
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf[t] = _mm_setzero_si128();
		mirror_psf_2[t] = _mm_setzero_si128();
	}
	for (int t = -filter_length; t <= filter_length; ++t) {
		mirror_psf[t] = _mm_set1_epi32(filter_cof[t]);
		mirror_psf_2[t] = _mm_set1_epi32(shifted_filter_cof[t]);
	}
	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((this->height <= in->height) && (this->width <= in->width));
	assert(((stride & 3) == 0) && ((in->stride & 3) == 0));
	int vec_stride_in = in->stride / 4;
	int vec_stride_out = this->stride / 4;
	int vec_width_out = (this->width + 3) / 4; // Big enough to cover the width

	// Do the filtering
	__m128i *line_out = (__m128i *) buf;
	__m128i *line_in = (__m128i *)(in->buf);
	for (int r = 0; r < height; r++,
		line_out += vec_stride_out, line_in += vec_stride_in)
		for (int c = 0; c < vec_width_out; c++){
			__m128i *ip = (line_in + c) - vec_stride_in*filter_length;
			__m128i sum = _mm_setzero_si128();
			for (int y = -filter_length; y <= filter_length; y++, ip += vec_stride_in)
				sum = _mm_add_epi32(sum, _mm_mullo_epi32(mirror_psf[y], *ip));
			line_out[c] = _mm_srai_epi32(sum, 15);
		}
}
