/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include "motion.h"
#include "time.h"
#include "cmath"

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */
#define nominal_block_size 8 
/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	int *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	int *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	int *left_edge = buf - border*stride;
	int *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}

/*****************************************************************************/
/* STATIC                         find_motion                                */
/*****************************************************************************/

static mvector
find_motion(my_image_comp *ref, my_image_comp *tgt,
	int start_row, int start_col, int block_width, int block_height, int pixel)
	/* This function finds the motion vector which best describes the motion
	   between the `ref' and `tgt' frames, over a specified block in the
	   `tgt' frame.  Specifically, the block in the `tgt' frame commences
	   at the coordinates given by `start_row' and `start_col' and extends
	   over `block_width' columns and `block_height' rows.  The function finds
	   the translational offset (the returned vector) which describes the
	   best matching block of the same size in the `ref' frame, where
	   the "best match" is interpreted as the one which minimizes the sum of
	   absolute differences (SAD) metric. */
{
	mvector vec, best_vec;
	float resolution = 1 / (float)pixel;
	int sad, best_sad = INT_MAX;
	for (vec.y = -8; vec.y <= 8; vec.y+= resolution)
		for (vec.x = -8; vec.x <= 8; vec.x+= resolution)
		{
			int ref_row = ((float)pixel)*((float)start_row - vec.y);
			int ref_col = ((float)pixel)*((float)start_col - vec.x);
			if ((ref_row < 0) || (ref_col < 0) ||
				((ref_row + pixel * block_height) > ref->height) ||
				((ref_col + pixel * block_width) > ref->width))
				continue; // Translated block not containe within reference frame
			int r, c;
			int *rp = ref->buf + ref_row*ref->stride + ref_col;
			int *tp = tgt->buf + start_row*tgt->stride + start_col;
			for (sad = 0, r = block_height; r > 0; r--,
				rp += (pixel*ref->stride), tp += tgt->stride)
				for (c = 0; c < block_width; c++)
				{
					int diff = tp[c] - rp[pixel*c];
					sad += (diff < 0) ? (-diff) : diff;
				}
			if (sad < best_sad)
			{
				best_sad = sad;
				best_vec = vec;
			}
		}

	return best_vec;
}

/*****************************************************************************/
/* STATIC                         motion_comp                                */
/*****************************************************************************/

static void
motion_comp(my_image_comp *ref, my_image_comp *tgt, mvector vec,
	int start_row, int start_col, int block_width, int block_height, int pixel)
	/* This function transfers data from the `ref' frame to a block within the
	   `tgt' frame, thereby realizing motion compensation.  The motion in
	   question has already been found by `find_motion' and is captured by
	   the `vec' argument.  The block in the `tgt' frame commences
	   at the coordinates given by `start_row' and `start_col' and extends
	   for `block_width' columns and `block_height' rows. */
{
	int r, c;
	int ref_row = ((float)pixel)*((float)start_row - vec.y);
	int ref_col = ((float)pixel)*((float)start_col - vec.x);
	int *rp = ref->buf + ref_row*ref->stride + ref_col;
	int *tp = tgt->buf + start_row*tgt->stride + start_col;
	for (r = block_height; r > 0; r--,
		rp += (pixel*ref->stride), tp += tgt->stride)
		for (c = 0; c < block_width; c++)
			tp[c] = rp[pixel*c];
}

static void draw_line(my_image_comp *ref, my_image_comp *tgt, mvector vec,
	int start_row, int start_col, int block_width, int block_height) {
	int stride = tgt->stride;
	int r = block_height >> 1, c = block_width >> 1;
	int *center = tgt->buf + (r + start_row) * tgt->stride + start_col + c;
	int x = vec.x, y = vec.y;

	if (x == 0 && y > 0) {
		for (int index = 0; index <= y; ++index) {
			center[index * stride] = 0;
		}
	}
	else if (x == 0 && y < 0) {
		for (int index = y; index <= 0; ++index) {
			center[index * stride] = 0;
		}
	}

	if (y == 0 && x > 0) {
		for (int index = 0; index <= x; ++index) {
			center[index] = 0;
		}
	}
	else if (y == 0 && x < 0) {
		for (int index = x; index <= 0; ++index) {
			center[index] = 0;
		}
	}

	if (y > 0 && x > 0) {
		if (vec.y >= x) {
			int block_each_x = (y + x - 1) / x;
			center[y * stride + x] = 0;
			for (c = 0, r = 0; c < x; ++c) {
				for (int index = block_each_x; index > 0; --index, ++r) {
					if (r > y)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
		else {
			int block_each_y = (x + y - 1) / y;
			center[y * stride + x] = 0;
			for (r = 0, c = 0; r < y; ++r) {
				for (int index = block_each_y; index > 0; --index, ++c) {
					if (c > x)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
	}

	if (x < 0 && y < 0) {
		if (-y >= -x) {
			int block_each_x = (-y - x - 1) / (-x);
			center[y * stride + x] = 0;
			for (c = 0, r = 0; c > x; --c) {
				for (int index = block_each_x; index > 0; --index, --r) {
					if (r < y)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
		else {
			int block_each_y = (-x - y - 1) / -y;
			center[y * stride + x] = 0;
			for (r = 0, c = 0; r > y; --r) {
				for (int index = block_each_y; index > 0; --index, --c) {
					if (c < x)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
	}

	if (y < 0 && x > 0) {
		if (-y >= x) {
			int block_each_x = (-y + x - 1) / x;
			center[y * stride + x] = 0;
			for (c = 0, r = 0; c < x; ++c) {
				for (int index = block_each_x; index > 0; --index, --r) {
					if (r < y)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
		else {
			int block_each_y = (x - y - 1) / (-y);
			center[y * stride + x] = 0;
			for (r = 0, c = 0; r > y; --r) {
				for (int index = block_each_y; index > 0; --index, ++c) {
					if (c > x)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
	}

	if (y > 0 && x < 0) {
		if (y >= -x) {
			int block_each_x = (y - x - 1) / (-x);
			center[y * stride + x] = 0;
			for (c = 0, r = 0; c > x; --c) {
				for (int index = block_each_x; index > 0; --index, ++r) {
					if (r > y)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
		else {
			int block_each_y = (-x + y - 1) / y;
			center[y * stride + x] = 0;
			for (r = 0, c = 0; r < y; ++r) {
				for (int index = block_each_y; index > 0; --index, --c) {
					if (c < x)
						break;
					center[r * stride + c] = 0;
				}
			}
		}
	}

}

static void interpolation(my_image_comp *in, my_image_comp *out) {

	for (int r = 0; r < out->height; ++r)
		for (int c = 0; c < out->width; ++c) {
			out->buf[r  * out->stride + c] = 0;
		}

	out->perform_boundary_extension();

	for (int r = 0; r < in->height; ++r)
		for (int c = 0; c < in->width; ++c) {
			int k = r << 1;
			int h = c << 1;
			out->buf[k * out->stride + h] = in->buf[r * in->stride + c];
		}

	
	for (int r = 0; r < out->height; r+=2)
		for (int c = 0; c < out->width; c+=2) {
			out->buf[r * out->stride + c + 1] = out->buf[r * out->stride + c]\
				 + out->buf[r * out->stride + c + 2];
			out->buf[r * out->stride + c + 1] = out->buf[r * out->stride + 1 + c] >> 1;
		}

	for (int r = 0; r < out->height+2; r+=2)
		for (int c = 0; c < out->width+2; ++c) {
			out->buf[(r+1) * out->stride + c] = out->buf[r * out->stride + c]\
				+ out->buf[(r+2) * out->stride + c];
			out->buf[(r+1) * out->stride + c] = out->buf[(r+1) * out->stride + c] >> 1;
		}

	
/*
	FILE *fp;
	fopen_s(&fp, "Logdata.m", "w");
	fprintf_s(fp, "mypic = [");
	for (int r = 0; r < out->height; ++r) {
		for (int c = 0; c < out->width; ++c) {
			fprintf_s(fp, "%d, ",
				out->buf[ r * in->height + c]);
		}
		fprintf_s(fp, ";\n");
	}
	fprintf_s(fp, "];");
	*/
	
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
		fprintf(stderr,
			"Usage: %s <bmp frame 1> <bmp frame 2> <bmp MC out>\n",
			argv[0]);
		return -1;
	}

	clock_t start_time = clock();

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in[2];
		if ((err_code = bmp_in__open(&in[0], argv[1])) != 0)
			throw err_code;
		if ((err_code = bmp_in__open(&in[1], argv[2])) != 0)
			throw err_code;

		int width = in[0].cols, height = in[0].rows;
		if ((width != in[1].cols) || (height != in[1].rows))
		{
			fprintf(stderr, "The two input frames have different dimensions.\n");
			return -1;
		}
		my_image_comp mono[2];
		mono[0].init(height, width, 4); // Leave a border of 4 (in case needed)
		mono[1].init(height, width, 4); // Leave a border of 4 (in case needed)

		mono[0].perform_boundary_extension();

		int n, r, c;
		int num_comps = in[0].num_components;
		io_byte *line = new io_byte[width*num_comps];
		for (n = 0; n < 2; n++)
		{
			for (r = height - 1; r >= 0; r--)
			{ // "r" holds the true row index we are reading, since the image
			  // is stored upside down in the BMP file.
				if ((err_code = bmp_in__get_line(&(in[n]), line)) != 0)
					throw err_code;
				io_byte *src = line; // Points to first sample of component n
				int *dst = mono[n].buf + r * mono[n].stride;
				for (c = 0; c < width; c++, src += num_comps)
					dst[c] = *src;
			}
			bmp_in__close(&(in[n]));
		}

		// Allocate storage for the motion compensated output
		my_image_comp output;
		output.init(height, width, 0); // Don't need a border for output

		my_image_comp ref_frame;
		ref_frame.init(height * 2, width * 2, 4);
		
		interpolation(&mono[0], &ref_frame);


		my_image_comp output_comps[3];

		for (int i = 0; i < 3; ++i) {
			output_comps[i].init(height, width, 0);
		}

		for (int r = 0; r < height; ++r)
			for (int c = 0; c < width; ++c) {
				output_comps[0].buf[r * output_comps[0].stride + c] = (mono[1].buf[r * mono[1].stride + c] >> 1) + 128;
				output_comps[1].buf[r * output_comps[1].stride + c] = (mono[1].buf[r * mono[1].stride + c] >> 1) + 128;
				output_comps[2].buf[r * output_comps[2].stride + c] = (mono[1].buf[r * mono[1].stride + c] >> 1) + 128;
			}

		// Now perform simple motion estimation and compensation
		int nominal_block_width = nominal_block_size;
		int nominal_block_height = nominal_block_size;
		int block_width, block_height;
		FILE *fp;
		fopen_s(&fp, "Logdata.m", "w");
		fprintf_s(fp, "myvec = [");
		for (r = 0; r < height; r += block_height)
		{
			block_height = nominal_block_height;
			if ((r + block_height) > height)
				block_height = height - r;
			for (c = 0; c < width; c += block_width)
			{
				block_width = nominal_block_width;
				if ((c + block_width) > width)
					block_width = width - c;
				mvector vec = find_motion(&(ref_frame), &(mono[1]),
					r, c, block_width, block_height, 2);
				
				fprintf_s(fp, "%f,%f;\n",vec.x, vec.y);
				
				motion_comp(&(ref_frame), &output, vec,
					r, c, block_width, block_height, 2);
				draw_line(&(mono[1]), &output_comps[1], vec,
					r, c, block_width, block_height);
			
			}
		}
		fprintf_s(fp, "];");
		/*
		for (int r = 0; r < height; ++r)
			for (int c = 0; c < width; ++c) {
				output_comps[0].buf[r * output_comps[0].stride + c] = (output.buf[r * output.stride + c]);
				output_comps[1].buf[r * output_comps[1].stride + c] = (output.buf[r * output.stride + c]);
				output_comps[2].buf[r * output_comps[2].stride + c] = (output.buf[r * output.stride + c]);
			}
		*/
		float sum = 0;
		for (int r = 0; r < height; ++r)
			for (int c = 0; c < width; ++c) {
				float diff = 0;
				diff = mono[1].buf[r * mono[1].stride + c] - output.buf[r * output.stride + c];
				diff *= diff;
				sum += diff;
			}

		sum = sum / (width * height);
		printf("The MSE is %f \r\n", sum);

		// Write the motion compensated image out
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[3], width, height, 3)) != 0)
			throw err_code;

		io_byte *out_line = new io_byte[width * 3];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < 3; n++)
			{
				io_byte *dst = out_line + n; // Points to first sample of component n
				int *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < width; c++, dst += 3) {

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
		delete[] line;

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
		return -1;
	}
	return 0;
}
