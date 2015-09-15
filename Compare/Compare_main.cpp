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

/*****************************************************************************/
/*                      void compute_difference                              */
/*****************************************************************************/

void compute_difference(my_aligned_image_comp *in1, my_aligned_image_comp *in2, my_aligned_image_comp *out)
{
	int width, height;
	width = in1->width;
	height = in1->height;
	float sum_diff = 0.0f;
	float mse = 0.0f;
	for (int i = 0; i < width; ++i)
		for (int j = 0; j < height; ++j) {
			float *source_1 = in1->buf + j * in1->width + i;
			float *source_2 = in2->buf + j * in2->width + i;
			float *dest = out->buf + j * width + i;
			float temp = *source_1 - *source_2;
			*dest = 0.5f * temp;
			sum_diff += temp;
			temp *= temp;
			mse += temp;
		}

	float rip_width = 1.0f / width;
	float rip_height = 1.0f / height;
	sum_diff *= rip_width;
	sum_diff *= rip_height;
	mse *= rip_width;
	mse *= rip_height;
	
	float PSNR = 0.0f;
	PSNR = 10.0f * log10f(255.0f * 255.0f / mse);
	printf("The Mean Error is %f \n The MSE is %f \n The PSNR is %f dB\n\r", sum_diff, mse, PSNR);

}
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
/*                        my_aligned_image_comp::filter                      */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int filter_length, int mode)
{
#define FILTER_EXTENT 14
#define FILTER_TAPS 2*FILTER_EXTENT+1
#define PI 3.141592653589793F

	// Create the vertical filter PSF as a local array on the stack.
	float filter_buf[FILTER_TAPS];
	float filter_buf_2[FILTER_TAPS];

	for (int i = 0; i < FILTER_TAPS; ++i) {
		filter_buf[i] = 0.0F;
		filter_buf_2[i] = 0.0F;
	}

	float *mirror_kernal_1 = filter_buf + FILTER_EXTENT;
	float *mirror_kernal_2 = filter_buf_2 + FILTER_EXTENT;

	// `mirror_kernal_1' points to the central tap in the filter
	for (int t = -filter_length; t <= filter_length; ++t) {
		mirror_kernal_2[t] = sinf(0.4F * PI * (t - .5F)) / (0.4F * PI * (t - .5F)) \
			* 0.5F * (1 + cosf( PI * (t - .5F) / (filter_length + 0.5F)));

		if (filter_length == 0)
			mirror_kernal_2[t] = 1;
		if (t == 0)
			mirror_kernal_1[t] = 1;
		else if (t > 0)
			mirror_kernal_1[t] = mirror_kernal_1[-t];
		else
			mirror_kernal_1[t] = sinf(0.4F * PI * t) / (0.4F * PI * t) \
				* 0.5F * (1 + cosf(PI * t / (filter_length + 0.5F)));
	}

	float gain_1 = 0, gain_2 = 0;
	for (int t = -filter_length; t <= filter_length; t++) {
		gain_1 += mirror_kernal_1[t];
		gain_2 += mirror_kernal_2[t];
	}

	gain_1 = 1 / gain_1;
	gain_2 = 1 / gain_2;

	for (int t = -filter_length; t <= filter_length; t++) {
		mirror_kernal_1[t] = mirror_kernal_1[t] * gain_1;
		mirror_kernal_2[t] = mirror_kernal_2[t] * gain_2;
	}

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
//	assert((this->height <= in->height) && (this->width <= in->width));

/*	float *intermed_handle = new float[in->stride];
	float *intermed_buf = intermed_handle + in->border;
*/
	// Perform the convolution
	if (mode == 1) {
		for (int r = 0, e = 0; r < height; r += 2, e += 5) {
			if (e > in->height)
				break;
			for (int c = 0; c < width; c++)
			{
				float *ip = in->buf + ((e + 2)*in->stride) + c;
				float *op = buf + (r + 1)*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; ++y)
					sum += ip[y] * mirror_kernal_2[y];
					//sum += ip[y] * mirror_kernal_1[y];
				*op = sum;
			}
		}

		for (int r = 0, e = 0; r < height; r += 2, e += 5) {
			if (e >  in->height)
				break;
			for (int c = 0; c < width; c++){
				float *ip = in->buf + (e*in->stride) + c;
				float *op = buf + r*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y] * mirror_kernal_1[y];
				*op = sum;
			}
		}
	}

	if (mode == 2) {
		for (int r = 0; r < height; ++r) {
			for (int c = 0, e = 0; c < width; c += 2, e += 5) {
				if (e > in->stride)
					break;
				float *ip = in->buf + r*in->stride + e ;
				float *op = buf + r*stride + c;
				float sum = 0.0F;
				for (int y = -filter_length; y <= filter_length; y++)
					sum += ip[y] * mirror_kernal_1[y];
				*op = sum;
			}

			for (int c = 0, e = 0; c < width; c += 2, e += 5) {
				if (e > in->stride)
					break;
				float *ip = in->buf + r*in->stride + e + 2;
				float *op = buf + r*stride + c + 1;
				float sum = 0.0F;
				for (int y = -filter_length ; y <= filter_length; y++)
					sum += ip[y] * mirror_kernal_2[y];
					//sum += ip[y] * mirror_kernal_1[y];
				*op = sum;
			}
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
		fprintf(stderr, "Usage: %s <in1 bmp file> <in2 bmp file> <out bmp file>\n", argv[0]);
		return -1;
	}

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		bmp_in in2;

		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;
		if ((err_code = bmp_in__open(&in2, argv[2])) != 0)
			throw err_code;

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;
		my_aligned_image_comp *input1 =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp *input2 =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++) {
			input1[n].init(height, width, 0); // Leave a border of 4
			input2[n].init(in2.rows, in2.cols, 0);
		}
		int r; // Declare row index
		io_byte *line = new io_byte[width * num_comps];
		io_byte *line2 = new io_byte[in2.cols * num_comps];
		
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			if ((err_code = bmp_in__get_line(&in2, line2)) != 0)
				throw err_code;

			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				io_byte *src2 = line2 + n;

				float *dst = input1[n].buf + r * input1[n].stride;
				float *dst2 = input2[n].buf + r * input1[n].stride;
				
				for (int c = 0; c < width; c++, src += num_comps) {
					dst[c] = (float)*src;
					dst2[c] = (float)*src2;
				}
			}
		}

		bmp_in__close(&in);
		bmp_in__close(&in2);

		// Allocate storage for the filtered output
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];

		for (n = 0; n < num_comps; ++n)
			output_comps[n].init(height, width, 0);

		for (n = 0; n < num_comps; n++) {	
			compute_difference(input1+n, input2+n, output_comps + n);
		}
		
		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[3], width, height, num_comps)) != 0)
			throw err_code;
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < width; c++, dst += num_comps) {
					if (src[c] > 255) {
						src[c] = 255;
					}
					else if (src[c] < 0) {
						src[c] = 0;
					}
//					printf("%6.3f ",src[c]);
					*dst = (io_byte)src[c];
				}
// The cast to type "io_byte" is
// required here, since floats cannot generally be
// converted to bytes without loss of information. The
// compiler will warn you of this if you remove the cast.
// There is in fact not the best way to do the
// conversion. You should fix it up in the lab.
			}
			bmp_out__put_line(&out, line);
		}
		bmp_out__close(&out);
		delete[] line;
		delete[] line2;
		delete[] input1;
		delete[] input2;
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
