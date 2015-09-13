/*****************************************************************************/
// File: vertical_filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/
#include <time.h>
#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>


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
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border*stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}

/*****************************************************************************/
/*                        filter coefficient generation                      */
/*****************************************************************************/

void my_aligned_image_comp::filter_generation(float *mirror_psf, int filter_length, float shift)
{
#define FILTER_EXTENT 14
	assert((2 * filter_length + 1) > FILTER_EXTENT);
	assert(abs(shift) > 1);

	if (shift == 0)
		mirror_psf[0] = 1;
	else if (shift < 0) {
		mirror_psf[0] = (1 + shift);
		mirror_psf[-1] = -shift;
	}
	else {
		mirror_psf[0] = (1 - shift);
		mirror_psf[1] = shift;
	}	
}

/*****************************************************************************/
/*                        my_aligned_image_comp::filter                      */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int filter_length, int mode)
{
#define FILTER_EXTENT 1
#define FILTER_TAPS 2*FILTER_EXTENT+1


	// Create the vertical filter PSF as a local array on the stack.
	float *filter_buf_1 = new float[FILTER_TAPS];
	float *filter_buf_2 = new float[FILTER_TAPS];
	float *filter_buf_3 = new float[FILTER_TAPS];
	float *filter_buf_4 = new float[FILTER_TAPS];
	float *filter_buf_5 = new float[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf_1[i] = 0.0F;
		filter_buf_2[i] = 0.0F;
		filter_buf_3[i] = 0.0F;
		filter_buf_4[i] = 0.0F;
		filter_buf_5[i] = 0.0F;
	}

	float *mirror_kernal_1 = filter_buf_1 + FILTER_EXTENT;
	float *mirror_kernal_2 = filter_buf_2 + FILTER_EXTENT;
	float *mirror_kernal_3 = filter_buf_3 + FILTER_EXTENT;
	float *mirror_kernal_4 = filter_buf_4 + FILTER_EXTENT;
	float *mirror_kernal_5 = filter_buf_5 + FILTER_EXTENT;

	filter_generation(mirror_kernal_1, filter_length, 0);
	filter_generation(mirror_kernal_2, filter_length, 0.4F);
	filter_generation(mirror_kernal_3, filter_length, -0.2F);
	filter_generation(mirror_kernal_4, filter_length, 0.2F);
	filter_generation(mirror_kernal_5, filter_length, -0.4F);

/*	for (int t = -filter_length; t < filter_length; ++t)
		printf("t = %d, 1:%f;2:%f;3:%f;4:%f;5:%f\n", t, mirror_kernal_1[t],
			mirror_kernal_2[t], mirror_kernal_3[t], mirror_kernal_4[t],
			mirror_kernal_5 [t]
			);
*/
	float *filter_kernal[5] =
			{mirror_kernal_1, mirror_kernal_2, mirror_kernal_3, mirror_kernal_4, mirror_kernal_5};


	// Check for consistent dimensions
	assert(in->border > FILTER_EXTENT);
//	assert((this->height <= in->height) && (this->width <= in->width));

#define verticle_mode 1
#define horizon_mode 2
#define numKernals 5

	// Perform the convolution
	if (mode == verticle_mode) {
		for (int r = 0, e = 0; r < height; r += 5, e += 2) {
			for (int c = 0; c < in->width; c++){
				for (int w = 0, offset = 0; w < numKernals; ++w) {
					if (r + w >  height )
						break;
					if ((w == 2)||(w == 3))
						offset = 1;
					else if (w == 4)
						offset = 2;
					float *ip = in->buf + ((e + offset)*in->stride) + c;
					float *op = this->buf + (r + w)*stride + c;
					float sum = 0.0F;
					for (int y = -filter_length; y <= filter_length; ++y) {
						sum += ip[y] * filter_kernal[w][y];
					}
					*op = sum;
				}
			}
		}
	}

	if (mode == horizon_mode) {

		for (int r = 0; r < in->height; ++r) {
			for (int c = 0, e = 0; c < width; c += 5, e += 2) {
				for (int w = 0, offset = 0; w < numKernals; ++w) {
					if (w + c > width)
						break;
					if ((w == 2) || (w == 3))
						offset = 1;
					else if (w == 4)
						offset = 2;
					float *ip = in->buf + r*in->stride + e + offset;
					float *op = buf + r*stride + c + w;
					float sum = 0.0F;
					for (int y = -filter_length; y <= filter_length; y++)
						sum += ip[y] * filter_kernal[w][y];
					*op = sum;
				}
			}
		}
	}
	delete[] filter_buf_1;
	delete[] filter_buf_2;
	delete[] filter_buf_3;
	delete[] filter_buf_4;
	delete[] filter_buf_5;
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
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <filter taps>\n", argv[0]);
		return -1;
	}

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int filter_tap = atoi(argv[3]);
		filter_tap = 1;

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;
		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 1); // Leave a border of 4

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
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps)
					dst[c] = (float)*src; // The cast to type "float" is not
						  // strictly required here, since bytes can always be
						  // converted to floats without any loss of information.
			}
		}
		delete[] line;
		bmp_in__close(&in);

		// Allocate storage for the filtered output
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermediea_comps =
			new my_aligned_image_comp[num_comps];

		int out_width, out_height;
		out_width = (int)ceilf(width * 2.5F);
		out_height = (int)ceilf(height * 2.5F);

		for (n = 0; n < num_comps; n++) {
			intermediea_comps[n].init(out_height, width, 1);
			output_comps[n].init(out_height, out_width, 0); // Don't need a border for output
		}

		  // Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
			intermediea_comps[n].perform_boundary_extension();
		}

		printf("Start filtering!\r\n");

#define verticle_mode 1
#define horizon_mode 2
		for (int j = 0; j < 1; ++j) {
#if 1
			for (n = 0; n < num_comps; n++)
				intermediea_comps[n].filter(input_comps + n, filter_tap, verticle_mode);

			delete[] input_comps;
			for (n = 0; n < num_comps; n++) {	
				intermediea_comps[n].perform_boundary_extension();
				output_comps[n].filter(intermediea_comps + n, filter_tap, horizon_mode);
			}
#else
			for (n = 0; n < num_comps; n++)
				output_comps[n].vector_filter(input_comps + n);
#endif
		}

		delete[] intermediea_comps;
		printf("Filtering end!\r\n");

		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], out_width, out_height, num_comps)) != 0)
			throw err_code;

		io_byte *out_line = new io_byte[out_width*num_comps];

		for (r = out_height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = out_line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < out_width; c++, dst += num_comps) {
					if (src[c] > 255) {
						src[c] = 255;
					}
					else if (src[c] < 0) {
						src[c] = 0;
					}
						
					*dst = (io_byte)src[c];
				}
// The cast to type "io_byte" is
// required here, since floats cannot generally be
// converted to bytes without loss of information. The
// compiler will warn you of this if you remove the cast.
// There is in fact not the best way to do the
// conversion. You should fix it up in the lab.
			}
			bmp_out__put_line(&out, out_line);
		}	
		
		
		bmp_out__close(&out);
		delete[] out_line;		
		delete[] output_comps;
		clock_t end_time = clock();
		float elaps = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
		printf_s("The runtime is %f seconds!\n\r", elaps);
		




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
