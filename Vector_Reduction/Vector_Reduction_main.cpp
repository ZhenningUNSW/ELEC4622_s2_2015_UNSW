/*****************************************************************************/
// File: vertical_filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "aligned_image_comps.h"
#include "math.h"
#include "time.h"

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*              my_aligned_image_comp::perform_boundary_extension            */
/*****************************************************************************/

void my_aligned_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	__int32 *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	__int32 *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	__int32 *left_edge = buf - border*stride;
	__int32 *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}

/*****************************************************************************/
/*                        my_aligned_image_comp::filter                      */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int filter_length)
{
#define FILTER_EXTENT 14
#define FILTER_TAPS (2*FILTER_EXTENT+1)
	__int32 shifted_filter_cof[FILTER_TAPS] = { 0, -9,	0,	88,	103, -171, -427,
		0, 886, 763, -1063, -2445,	0,	6438,
		12221,	12221,	6438,	0, -2445, -1063	,
		763, 886, 0, -427 ,-171,103,88, 0, -9 };
	__int32 filter_cof[FILTER_TAPS] = { -2, -12, 37, 124, 0, -347, -321,	461,
		1047,	0, -2042, -1835,	2923,	9800,	13102,
		9800,	2923, -1835, -2042,	0,	1047,	461, -321,
		-347,	0,	124,	37, -12, -2 };
	// Create the vertical filter PSF as a local array on the stack.
	__int32 filter_buf[FILTER_TAPS];
	__int32 filter_buf_2[FILTER_TAPS];

	__int32 *mirror_psf = filter_buf + FILTER_EXTENT;
	__int32 *mirror_psf_2 = filter_buf_2 + FILTER_EXTENT;
	__int32 *p1 = filter_cof + FILTER_EXTENT;
	__int32 *p2 = shifted_filter_cof + FILTER_EXTENT;
	// `mirror_psf' points to the central tap in the filter
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf[t] = 0;
		mirror_psf_2[t] = 0;
	}

	for (int t = -filter_length; t <= filter_length; t++) {
		mirror_psf[t] = p1[t];
		mirror_psf_2[t] = p2[t];
	}
	
	if (filter_length == 0) {
		mirror_psf[filter_length] = 1 << 15;
		mirror_psf_2[filter_length] = 1 << 15;
	}
/*	else {
		int sum1 = 0, sum2 = 0;
	    for (int t = -filter_length; t <= filter_length; t++) {
			sum1 += mirror_psf[t];
			sum2 += mirror_psf_2[t];
		}
		sum1 = sum1 >> 15;
		sum2 = sum2 >> 15;
		for (int t = -filter_length; t <= filter_length; t++) {
			if (sum1 != 0)
				mirror_psf[t] /= sum1;
			if (sum2 != 0)
				mirror_psf_2[t] /= sum2;
		}
	}
*/

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((this->height <= in->height) && (this->width <= in->width));

	my_aligned_image_comp *intermedia = new my_aligned_image_comp;

	intermedia->init(this->height, in->width, 15);

	// Perform the convolution
	for (int r = 0, e = 0; r < intermedia->height; r += 2, e += 5) {
		for (int c = 0; c < intermedia->width; c++)
		{
			__int32 *ip = in->buf + ((e + 2)*in->stride) + c;
			__int32 *op = intermedia->buf + (r + 1)*intermedia->stride + c;
			__int32 sum = 0;
			for (int y = -filter_length; y <= filter_length; ++y)
				sum += ip[y*in->stride] * mirror_psf_2[y];
			*op = sum >> 15;
		}
	}

	for (int r = 0, e = 0; r < intermedia->height; r += 2, e += 5) {
		for (int c = 0; c < intermedia->width; c++) {
			__int32 *ip = in->buf + (e*in->stride) + c;
			__int32 *op = intermedia->buf + r*intermedia->stride + c;
			__int32 sum = 0;
			for (int y = -filter_length; y <= filter_length; y++)
				sum += ip[y*in->stride] * mirror_psf[y];
			*op = sum >> 15;
		}
	}

	intermedia->perform_boundary_extension();

	for (int r = 0; r < this->height; ++r) {
		for (int c = 0, e = 0; c < this->width; c += 2, e += 5) {
			__int32 *ip = intermedia->buf + r*intermedia->stride + e;
			__int32 *op = buf + r*stride + c;
			__int32 sum = 0;
			for (int y = -filter_length; y <= filter_length; y++)
				sum += ip[y] * mirror_psf[y];
			*op = sum >> 15;
		}
	}

	for (int r = 0; r < this->height; ++r) {
		for (int c = 0, e = 0; c < this->width; c += 2, e += 5) {
			__int32 *ip = intermedia->buf + r*intermedia->stride + e + 2;
			__int32 *op = buf + r*stride + c + 1;
			__int32 sum = 0;
			for (int y = -filter_length; y <= filter_length; y++)
				sum += ip[y] * mirror_psf_2[y];
			*op = sum >> 15;
		}
	}
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file>\n", argv[0]);
		return -1;
	}

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int filter_length = 0;
		filter_length = atoi(argv[3]);

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;
		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 15); // Leave a border of 4

		int r; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				__int32 *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps)
					dst[c] = (__int32)*src; // The cast to type "float" is not
						  // strictly required here, since bytes can always be
						  // converted to floats without any loss of information.
			}
		}
		bmp_in__close(&in);
		int out_width = (int)ceilf(width * 0.4f);
		int out_height = (int)ceilf(height * 0.4f);
		// Allocate storage for the filtered output
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			output_comps[n].init(out_height, out_width, 0); // Don't need a border for output

		  // Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++)
			input_comps[n].perform_boundary_extension();
		printf("Strat filtering!\n");
		for (int j = 0; j < 1000; ++j) {
#if 1
			for (n = 0; n < num_comps; n++)
				output_comps[n].filter(input_comps + n, filter_length);
#else
			for (n = 0; n < num_comps; n++)
				output_comps[n].vector_filter(input_comps + n, filter_length);
#endif
		}
		printf("Filtering finished!\n");
		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], out_width, out_height, num_comps)) != 0)
			throw err_code;
		for (r = out_height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = line + n; // Points to first sample of component n
				__int32 *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < out_width; c++, dst += num_comps) {
					if (src[c] > 255)
						src[c] = 255;
					else if (src[c] < 0)
						src[c] = 0;
					*dst = (io_byte)src[c];
				}// The cast to type "io_byte" is
						// required here, since floats cannot generally be
						// converted to bytes without loss of information.  The
						// compiler will warn you of this if you remove the cast.
						// There is in fact not the best way to do the
						// conversion.  You should fix it up in the lab.
			}
			bmp_out__put_line(&out, line);
	}
		bmp_out__close(&out);
		delete[] line;
		delete[] input_comps;
		delete[] output_comps;

		clock_t end_time = clock();
		float elaps = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
		printf_s("The runtime is %f seconds!\n\r", elaps);
		system("pause");

}
	catch (int exc) {
		if (exc == IO_ERR_NO_FILE)
			fprintf(stderr, "Cannot open supplied input or output file.\n");
		else if (exc == IO_ERR_FILE_HEADER)
			fprintf(stderr, "Error encountered while parsing BMP file header.\n");
		else if (exc == IO_ERR_UNSUPPORTED)
			fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
				"simple example supports only 8-bit and 24-bit data.\n");
		else if (exc == IO_ERR_FILE_TRUNC)
			fprintf(stderr, "Input or output file truncated unexpectedly.\n");
		else if (exc == IO_ERR_FILE_NOT_OPEN)
			fprintf(stderr, "Trying to access a file which is not open!(?)\n");
		return -1;
	}
	return 0;
}
