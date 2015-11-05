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
#include "dft.h"
#include <complex>

#define X_ECCESS -6
#define Y_ECCESS -7
#define NOT_A_ODD_NUMBER -8

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
		for (c = 0; c < width; c++) {		
			first_line[-r*stride + c] = first_line[r*stride + c];
		}
	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[-r*stride + c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border*stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[c];
			right_edge[c] = right_edge[-c];
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
		mirror_kernal_2[t] = 0.4f*sinf( 0.4f *PI * ( t - .5F)) / (0.4f *PI *  (t - .5F)) \
			* 0.5F * (1 + cosf( PI * (t - .5F) / (filter_length + 0.5F)));

		if (filter_length == 0)
			mirror_kernal_2[t] = 1;
		if (t == 0)
			mirror_kernal_1[t] = 0.4f * 1;
		else if (t > 0)
			mirror_kernal_1[t] = mirror_kernal_1[-t];
		else
			mirror_kernal_1[t] = 0.4f*sinf(0.4f *PI * t) / (0.4f * PI * t) \
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
					sum += ip[y*in->stride] * mirror_kernal_2[y];
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
					sum += ip[y*in->stride] * mirror_kernal_1[y];
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
/*

	//Create a single line buffer and copy the value form the intermediate result
	for (int j = 0; j < this->height; ++j) {
		for (int i = 0; i < this->width; ++i)
			intermed_buf[i] = this->buf[j * this->stride + i];

		//Perform boundry extension with mirror padding
		for (int i = 0; i < in->border; ++i) {
			intermed_buf[-i - 1] = intermed_buf[i];
			intermed_buf[i + this->width] = intermed_buf[-i + this->width - 1];
		}

		//Perform convolution
		for (int c = 0; c < width; c++)
		{
			//float *ip = in->buf + r*in->stride + c;
			float *op = this->buf + j*stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
				sum += intermed_buf[c + y] * mirror_kernal_1[y];
			*op = sum;
		}
	}
*/
}

void my_aligned_image_comp::hanning_window()
{
	int half_N = width / 2;
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			buf[i * stride + j] = buf[i * stride + j] * \
				0.25F * (1 + cosf(PI * (j - half_N) / half_N)) *\
				(1 + cosf(PI * (i - half_N) / half_N));
		}
	}
	
}

void my_aligned_image_comp::mean_substract()
{
	float mean = 0, sum = 0;

	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j){
			sum += buf[i * stride + j];	
		}

	mean = sum / width / height;
	int half_N = width / 2;
	for (int i = 0; i < height; ++i)
		for (int j = 0; j < width; ++j) {
			buf[i * stride + j] -= mean;
		}
}

void my_aligned_image_comp::perform_dft(float alpha)
{
	float *dft_real = new float[height*width];
	float *dft_imag = new float[height*width];

	// Process the image, plane by plane
		// First copy all samples to the `dft_real' buffer
	for (int r = 0; r < height; r++)
		for (int c = 0; c < width; c++)
		{
			dft_real[r*width + c] = this->buf[r*stride + c];
			dft_imag[r*width + c] = 0.0F;
		}
		// Next, perform the 2D DFT

		// Put your code here
		my_direct_dft *dft = new my_direct_dft[1];
		dft->init(height, 1);
		for (int c = 0; c < width; c++)
			dft->perform_transform((dft_real + c), (dft_imag + c), width);

		for (int r = 0; r < height; r++)
			dft->perform_transform((dft_real + r * width), (dft_imag + r * width), 1);

		// Write DFT magnitudes to the output image, possibly using a log
		// function to see the values more clearly.

		// Put your code here\

		// DFT magniture divided by N^2

		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++){
				buf[r*stride + c] = dft_imag[r*width + c] * dft_imag[r*width + c] 
								  + dft_real[r*width + c] * dft_real[r*width + c];
				buf[r*stride + c] /= (width * width);
			}

		if ( alpha > 0) {
			for (int r = 0; r < height; r++)
				for (int c = 0; c < width; c++)
					buf[r*stride + c] *= alpha;
		}

/*		FILE *fp;
		fopen_s(&fp, "Logdata.m", "w");
		fprintf_s(fp, "[");
		for (int r = 0; r < height; ++r) {
			for (int c = 0; c < width; ++c) {
				fprintf_s(fp, "(%.3f)+(%.3fi) ",
					dft_real[r * width + c], dft_imag[r * width + c]);
			}
			fprintf_s(fp, ";\n");
		}
		fprintf_s(fp, "]");
		*/
		
		// Re-arrange the quadrants.
		// Reuse the buffer dft_real
		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++) 
				dft_real[r*width + c] = buf[r*stride + c];

		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
				buf[r*stride + c] = 0;

/*		FILE *fp;
		fopen_s(&fp, "Logdata.m", "w");
		fprintf_s(fp, "abspic = [");
		for (int r = 0; r < height; ++r) {
			for (int c = 0; c < width; ++c) {
				fprintf_s(fp, "%f,", dft_real[r*width + c]);
			}
			fprintf_s(fp, ";\n");
		}
		fprintf_s(fp, "]");
*/

		int half_N = width / 2;
		// Now start to re-arrange
		for (int r = 0; r < half_N; r++)
			for (int c = 0; c < half_N; c++)
				buf[r*stride + c] = dft_real[(r+half_N)*width + c + half_N];

		for (int r = half_N; r < width; r++)
			for (int c = half_N; c < width; c++)
				buf[r*stride + c] = dft_real[(r - half_N) *width + (c-half_N)];

		for (int r = 0; r < half_N; r++)
			for (int c = half_N; c < width; c++)
				buf[r*stride + c] = dft_real[(r+half_N) *width + (c-half_N)];

		for (int r = 0; r < half_N; r++)
			for (int c = half_N; c < width; c++)
				buf[c*stride + r] = dft_real[(c - half_N) *width + (r + half_N)];

}

void my_aligned_image_comp::perform_fft(float alpha)
{
	int N = 0;
	float f = (float)width;
	unsigned int const t = 1U << ((*(unsigned int *)&f >> 23) - 0x7f);
	N = t << (t < width);

	_complex *line = new _complex[N];
	_complex *image_buf = new _complex[N * N];

	for (int c = 0; c < N; ++c) {
		line[c].x = 0;
		line[c].y = 0;
	}

	for (int r = 0; r < N; ++r)
		for (int c = 0; c < N; ++c) {
			image_buf[r * N + c].x = buf[r * stride + c];
			image_buf[r * N + c].y = 0;

		}

/*	for (int r = 0; r < N; ++r) {
		image_buf[(N - 1) * N + r].x = image_buf[(N - 2) * N + r].x;
		image_buf[(N - 1) * N + r].y = image_buf[(N - 2) * N + r].y;
	}

	for (int r = 0; r < N; ++r) {
		image_buf[r * N + N - 1].x = image_buf[r * N + N - 2].x;
		image_buf[r * N + N - 1].y = image_buf[r * N + N - 2].y;
	}
*/
	for (int r = 0; r < N; ++r) {
		for (int c = 0; c < N; ++c) {
			line[c].x = image_buf[r * N + c].x;
			line[c].y = 0;
		}

		my_fft *fft = new my_fft;
//		fft->init(N, 1);
		line = fft->perform_transform(line, N);

		for (int c = 0; c < N; ++c) {
			image_buf[r * N + c].x = line[c].x;
			image_buf[r * N + c].y = line[c].y;
		}
	}
	
	for (int c = 0; c < N; ++c) {
		 for (int r = 0; r < N; ++r) {
			line[r].x = image_buf[r * N + c].x;
			line[r].y = image_buf[r * N + c].y;
		}

		my_fft *fft = new my_fft;
//		fft->init(N, 1);
		line = fft->perform_transform(line, N);

		for (int r = 0; r < N; ++r) {
			image_buf[r * N + c].x = line[r].x;
			image_buf[r * N + c].y = line[r].y;

		}
	}
/*
	FILE *fp;
	fopen_s(&fp, "Logdata.m", "w");
	fprintf_s(fp, "myfft=[");
	for (int c = 0; c < N; ++c) {
		for (int r = 0; r < N; ++r) {
			fprintf_s(fp, "(%f)+(%fi) ", image_buf[c * N + r].x, image_buf[c * N + r].y);
		}
		fprintf_s(fp, ";\n");
	}
	fprintf_s(fp, "];");
*/
	delete[] line;
/*
	_complex *sorted_image = new _complex[N * N];
	int half_N = N / 2;
	// Now start to re-arrange
	for (int r = 0; r < half_N; r++)
		for (int c = 0; c < half_N; c++)
			sorted_image[r * N + c] = image_buf[(r + half_N) * N + c + half_N];

	for (int r = half_N; r < N; r++)
		for (int c = half_N; c < N; c++)
			sorted_image[r * N + c] = image_buf[(r - half_N) * N + (c - half_N)];

	for (int r = 0; r < half_N; r++)
		for (int c = half_N; c < N; c++)
			sorted_image[r * N + c] = image_buf[(r + half_N) * N + (c - half_N)];

	for (int r = 0; r < half_N; r++)
		for (int c = half_N; c < N; c++)
			sorted_image[c * N + r] = image_buf[(c - half_N) * N + (r + half_N)];

	delete[] image_buf;
	*/
	for (int r = 0; r < width; r++)
		for (int c = 0; c < width; c++)
			buf[r * stride + c] = (image_buf[r * N + c].x * image_buf[r * N + c].x\
								+ image_buf[r * N + c].y * image_buf[r * N + c].y)\
								* alpha;

}

void my_aligned_image_comp::calculate_estimation(float *table)
{
	int N_16 = width / 16;
	int N_8 = width / 8;
	int N_4 = width / 4;
	int N_2 = width / 2;
	
	float *image_buf = new float[width * height];

	for (int r = 0; r < height; ++r)
		for (int c = 0; c < width; ++c) {
			image_buf[r * width + c] = buf[r * stride + c];
		}
	int width_new = width + 1;
	float *sorted_image = new float[width_new * width_new];
	float *center = sorted_image + N_2 * width_new + N_2;

	// Now start to re-arrange
	for (int r = 0; r <= N_2; r++)
		for (int c = 0; c <= N_2; c++)
			center[r * width_new + c] = image_buf[c + r * width];

	for (int r = 0; r <= N_2; r++)
		for (int c = -N_2; c < 0; c++)
			center[r * width_new + c] = image_buf[r * width + c + width];

	for (int r = -N_2; r < 0; r++)
		for (int c = 0; c <= N_2; c++)
			center[r * width_new + c] = image_buf[(r + width) * width + c];

	for (int r = -N_2; r < 0; r++)
		for (int c = -N_2; c < 0; c++)
			center[r * width_new + c] = image_buf[(r + width) * width + (c + width)];

	FILE *fp;
	fopen_s(&fp, "Logdata.m", "w");
	fprintf_s(fp, "myfft=[");
	for (int r = -N_2; r <= N_2; ++r) {
		for (int c = -N_2; c <= N_2; ++c) {
			fprintf_s(fp, "(%.7f)  ", center[r * width_new + c]);
		}
		fprintf_s(fp, ";\n");
	}
	fprintf_s(fp, "];");

	float s_r = 0, s_g = 0, s_b = 0, s_dc = 0;
	for (int r = -N_16; r <= N_16; ++r) {
		for (int c = -N_16; c <= N_16; ++c) {
			s_dc += center[r * width_new + c];
		}
	}

	for (int r = -N_8; r <= N_8; ++r) {
		for (int c = -N_8; c <= N_8; ++c) {
			s_b += center[r * width_new + c];
		}
	}

	for (int r = -N_4; r <= N_4; ++r) {
		for (int c = -N_4; c <= N_4; ++c) {
			s_r += center[r * width_new + c];
		}
	}

	for (int r = -N_2; r <= N_2; ++r) {
		for (int c = -N_2; c <= N_2; ++c) {
			s_g += center[r * width_new + c];
		}
	}

	float r = 0, g = 0, b = 0, temp;
	b = s_b - s_dc;
	temp = (2 * N_8 + 1) * (2 * N_8 + 1) - (2 * N_16 + 1) * (2 * N_16 + 1);
	temp = temp * temp;
	table[0] = b / temp;
	printf("The value of B = %f \r\n", table[0]);

	g = s_g - s_r;
	temp = (2 * N_2 + 1) * (2 * N_2 + 1) - (2 * N_4 + 1) * (2 * N_4 + 1);
	temp = temp * temp;
	table[1] = g / temp;
	printf("The value of G = %f \r\n", table[1]);

	r = s_r - s_b;
	temp = (2 * N_4 + 1) * (2 * N_4 + 1) - (2 * N_8 + 1) * (2 * N_8 + 1);
	temp = temp * temp;
	table[2] = r / temp;
	printf("The value of R = %f \r\n", table[2]);
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
	if (argc != 7)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <x> <y> <N> <alpha>\n", argv[0]);
		return -1;
	}
	
	printf("Program started.\r\n");

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int location_x = atoi(argv[3]);
		int location_y = atoi(argv[4]);
		int out_size = atoi(argv[5]);
		float alpha = atof(argv[6]);

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;

		if (out_size % 2 == 1)
		{
			err_code = NOT_A_ODD_NUMBER;
			throw err_code;
		}

		int half_outsize = out_size / 2;

		if (width < location_x) {
			err_code = X_ECCESS;
			throw err_code;
		}

		if (height < location_y) {
			err_code = Y_ECCESS;
			throw err_code;
		}

		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, half_outsize+1); // Leave a border of output size

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

		for (n = 0; n < num_comps; n++) {
			output_comps[n].init(out_size, out_size, 0);// Don't need a border for output
		}

		  // Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
		}

		for (int k = 0; k < num_comps; ++k)
			for (int i = -half_outsize; i < half_outsize; ++i)
				for (int j = -half_outsize; j < half_outsize; ++j) {
					output_comps[k].buf[(i+ half_outsize) * output_comps[k].stride + (j+ half_outsize)]
						= input_comps[k].buf[i * input_comps[k].stride + j 
						+ location_x + location_y * input_comps[k].stride];
				}

		for (n = 0; n < num_comps; n++) {
			output_comps[n].mean_substract();
			output_comps[n].hanning_window();
			//output_comps[n].perform_dft(alpha);
			output_comps[n].perform_fft(alpha);
			float *table = new float[3];
			output_comps[n].calculate_estimation(table);
		}


		 
		// Write the image back out again
		bmp_out out;
		io_byte *out_line = new io_byte[out_size*num_comps];
		if ((err_code = bmp_out__open(&out, argv[2], out_size, out_size, num_comps)) != 0)
			throw err_code;
		for (r = out_size - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = out_line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < out_size; c++, dst += num_comps) {
					
					if (src[c] > 255) {
						src[c] = 255;
					}
					else if (src[c] < 0) {
						src[c] = 0;
					}
						
					*dst = (io_byte)src[c];
				}
			}
			bmp_out__put_line(&out, out_line);
		}
		bmp_out__close(&out);
		
		delete[]out_line;
		
		delete[] input_comps;
		delete[] output_comps;

		clock_t end_time = clock();
		float elaps = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
		printf_s("The runtime is %f seconds!\n\r", elaps);
		system("Pause");

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
		else if (exc == X_ECCESS)
			fprintf(stderr, "The given location eccess x dimension!\n");
		else if (exc == Y_ECCESS)
			fprintf(stderr, "The given location eccess y dimension!\n");
		else if (exc == NOT_A_ODD_NUMBER)
			fprintf(stderr, "The given output dimension is not a odd number!\n");
		return -1;
	}
	return 0;
}
