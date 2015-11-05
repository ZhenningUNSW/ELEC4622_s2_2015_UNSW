/*****************************************************************************/
// File: io_bmp.h
// Author: David Taubman
// Last Revised: 30 July, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#ifndef IO_BMP_H
#define IO_BMP_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// Simple data types
typedef unsigned long io_uint32;
typedef long io_int32;
typedef unsigned char io_byte;

// Error codes
#define IO_ERR_NO_FILE          ((int) -1) /* If file not found */
#define IO_ERR_FILE_HEADER      ((int) -2) /* If header has an error */
#define IO_ERR_FILE_TRUNC       ((int) -3) /* If file ends unexpectely */
#define IO_ERR_UNSUPPORTED      ((int) -4) /* Exception code if file uses an
                                              unsupported format. */
#define IO_ERR_FILE_NOT_OPEN    ((int) -5) /* If trying to read/write a file
                                              which is not open, or has come
                                              to the end. */

// Structures defined here:
struct bmp_header;
struct bmp_in_state;
struct bmp_out_state;

/*****************************************************************************/
/*                              bmp_header                                   */
/*****************************************************************************/

struct bmp_header {
    io_uint32 size; // Size of this structure: must be 40
    io_int32 width; // Image width
    io_int32 height; // Image height; -ve means top to bottom.
    io_uint32 planes_bits; // Planes in 16 LSB's (must be 1); bits in 16 MSB's
    io_uint32 compression; // Only accept 0 here (uncompressed RGB data)
    io_uint32 image_size; // Can be 0
    io_int32 xpels_per_metre; // We ignore these
    io_int32 ypels_per_metre; // We ignore these
    io_uint32 num_colours_used; // Entries in colour table (0 means use default)
    io_uint32 num_colours_important; // 0 means all colours are important.
  };
  /* Notes:
        This header structure must be preceded by a 14 byte field, whose
     first 2 bytes contain the string, "BM", and whose next 4 bytes contain
     the length of the entire file.  The next 4 bytes must be 0. The final
     4 bytes provides an offset from the start of the file to the first byte
     of image sample data.
        If the bit_count is 1, 4 or 8, the structure must be followed by
     a colour lookup table, with 4 bytes per entry, the first 3 of which
     identify the blue, green and red intensities, respectively. */

/*****************************************************************************/
/*                                 bmp_in                                    */
/*****************************************************************************/

struct bmp_in {
    int num_components, rows, cols;
    int num_unread_rows;
    int line_bytes; // Number of bytes in each line, not including padding
    int alignment_bytes; // Bytes at end of each line to make a multiple of 4.
    FILE *in;
  };

extern int bmp_in__open(bmp_in *state, const char *fname);
  /* Opens the image file with the indicated name, initializing the supplied
     `state' structure to hold working state information for subsequent use
     with `bmp_in__close()' and `bmp_in__get_line()'.
        If an error occurs, the function returns one of the error codes
     `IO_ERR_NO_FILE', `IO_ERR_FILE_HEADER', `IO_ERR_FILE_TRUNC' or
     `IO_ERR_UNSUPPORTED' which are defined at the top of this header file.
        Otherwise, the function returns 0 for success. */

extern void bmp_in__close(bmp_in *state);
  /* You should use this function to close any image opened with
     `bmp_in__close'. */

extern int bmp_in__get_line(bmp_in *state, io_byte *line);
  /* Reads the next line of image data from the file opened using the most
     recent successful call to `bmp_in__open' (with the same `state'
     structure), storing the recovered samples in the supplied `line' buffer.
     Components are interleaved in BGR order within the `line' buffer.
        If successful, the function returns 0.  If the file terminates
     unexpectedly, the `IO_ERR_FILE_TRUNC' error code is returned.  If the
     file is not currently open, or the end has been reached, the
     `IO_ERR_FILE_NOT_OPEN' error code is returned. */

/*****************************************************************************/
/*                                bmp_out                                    */
/*****************************************************************************/

struct bmp_out {
    int num_components, rows, cols;
    int num_unwritten_rows;
    int line_bytes; // Number of bytes in each line, not including padding
    int alignment_bytes; // Number of 0's at end of each line.
    FILE *out;
  };

extern int bmp_out__open(bmp_out *state, const char *fname,
                         int width, int height, int num_components);
  /* Opens an image file with the indicated name for writing, initializing
     the supplied `state' structure to hold working state information for
     subsequent use with `bmp_out__close()' and `bmp_out__put_line()'.
        The `num_components' value should be 1 for a monochrome image and 3
     for a colour image.
        The function returns 0 if successful, `IO_ERR_NO_FILE' if the file
     cannot be opened, or else `IO_ERR_SUPPORTED' if an illegal combination
     of parameters is supplied. */

extern void bmp_out__close(bmp_out *state);
  /* You should use this function to close any image opened with
     `bmp_in__close'. */

extern int bmp_out__put_line(bmp_out *state, io_byte *line);
  /* Writes the next line of image data to the file opened using the most
     recent successful call to `bmp_out__open' (with the same `state'
     structure), writing the samples supplied via the `line' buffer.
     Components should be interleaved in BGR order within the `line' buffer.
        If successful, the function returns 0.  If the file cannot be written
     (e.g., the disk may be full), the `IO_ERR_FILE_TRUNC' error code is
     returned.  If the file is not currently open, or the end has been
     reached, the `IO_ERR_FILE_NOT_OPEN' error code is returned. */

#endif // IO_BMP_H
