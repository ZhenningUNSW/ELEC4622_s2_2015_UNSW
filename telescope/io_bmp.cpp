/*****************************************************************************/
// File: io_bmp.cpp
// Author: David Taubman
// Last Revised: 30 July, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <string.h>  // Import `memset' function
#include "io_bmp.h"

/* ========================================================================= */
/*                             Internal Functions                            */
/* ========================================================================= */

/*****************************************************************************/
/* STATIC                       to_little_endian                             */
/*****************************************************************************/

static void
  to_little_endian(io_int32 *words, int num_words)
{
  io_int32 test = 1;
  io_byte *first_byte = (io_byte *) &test;
  if (*first_byte)
    return; // Machine uses little-endian architecture already.
  io_int32 tmp;
  for (; num_words--; words++)
    {
      tmp = *words;
      *words = ((tmp >> 24) & 0x000000FF) +
               ((tmp >> 8)  & 0x0000FF00) +
               ((tmp << 8)  & 0x00FF0000) +
               ((tmp << 24) & 0xFF000000);
    }
}

/*****************************************************************************/
/* INLINE                      from_little_endian                            */
/*****************************************************************************/

static inline void
  from_little_endian(io_int32 *words, int num_words)
{
  to_little_endian(words,num_words);
}


/* ========================================================================= */
/*                                   bmp_in                                  */
/* ========================================================================= */

/*****************************************************************************/
/*                                bmp_in__open                               */
/*****************************************************************************/

int bmp_in__open(bmp_in *state, const char *fname)
{
  memset(state,0,sizeof(bmp_in)); // Start by reseting everything
  if ((state->in = fopen(fname,"rb")) == NULL)
    return(IO_ERR_NO_FILE);

  io_byte magic[14];
  bmp_header header;
  fread(magic,1,14,state->in);
  if ((magic[0] != 'B') || (magic[1] != 'M'))
    return(IO_ERR_FILE_HEADER);
  if (fread(&header,1,40,state->in) != 40)
    return(IO_ERR_FILE_TRUNC);

  from_little_endian((io_int32 *) &header,10);
  state->cols = header.width;
  state->rows = header.height;
  int bit_count = (header.planes_bits>>16);
  if (bit_count == 24)
    state->num_components = 3;
  else if (bit_count == 8)
    state->num_components = 1;
  else
    return(IO_ERR_UNSUPPORTED);
  int palette_entries_used = header.num_colours_used;
  if (state->num_components != 1)
    palette_entries_used = 0;
  else if (header.num_colours_used == 0)
    palette_entries_used = (1<<bit_count);
  int header_size = 54 + 4*palette_entries_used;

  int offset = magic[13];
  offset <<= 8; offset += magic[12];
  offset <<= 8; offset += magic[11];
  offset <<= 8; offset += magic[10];
  if (offset < header_size)
    return(IO_ERR_FILE_HEADER);
  if (palette_entries_used)
    fseek(state->in,palette_entries_used*4,SEEK_CUR); // Skip over palette
  if (offset > header_size)
    fseek(state->in,offset-header_size,SEEK_CUR);
  state->num_unread_rows = state->rows;
  state->line_bytes = state->num_components * state->cols;
  state->alignment_bytes =
    (4-state->line_bytes) & 3; // Pad to a multiple of 4 bytes
  return 0;
}

/*****************************************************************************/
/*                                bmp_in__close                              */
/*****************************************************************************/

void bmp_in__close(bmp_in *state)
{
  if (state->in != NULL)
    fclose(state->in);
  memset(state,0,sizeof(bmp_in));
}

/*****************************************************************************/
/*                              bmp_in__get_line                             */
/*****************************************************************************/

int bmp_in__get_line(bmp_in *state, io_byte *line)
{
  if ((state->in == NULL) || (state->num_unread_rows <= 0))
    return(IO_ERR_FILE_NOT_OPEN);
  state->num_unread_rows--;
  if (fread(line,1,(size_t) state->line_bytes,state->in) !=
      (size_t) state->line_bytes)
    return(IO_ERR_FILE_TRUNC);
  if (state->alignment_bytes > 0)
    {
      io_byte buf[3];
      fread(buf,1,(size_t) state->alignment_bytes,state->in);
    }
  return 0;
}


/* ========================================================================= */
/*                                  bmp_out                                  */
/* ========================================================================= */

/*****************************************************************************/
/*                               bmp_out__open                             */
/*****************************************************************************/

int bmp_out__open(bmp_out *state, const char *fname,
                  int width, int height, int num_components)
{
  memset(state,0,sizeof(bmp_out)); // Start by reseting everything
  state->num_components = num_components;
  state->rows = state->num_unwritten_rows = height;
  state->cols = width;
  io_byte magic[14];
  bmp_header header;
  int header_bytes = 14+sizeof(header);
  assert(header_bytes == 54);
  if (num_components == 1)
    header_bytes += 1024; // Need colour LUT.
  else if (num_components != 3)
    return(IO_ERR_UNSUPPORTED);
  state->line_bytes = num_components * width;
  state->alignment_bytes = (4-state->line_bytes) & 3;
  int file_bytes = header_bytes +
    (state->line_bytes+state->alignment_bytes)*state->rows;
  magic[0] = 'B'; magic[1] = 'M';
  magic[2] = (io_byte) file_bytes;
  magic[3] = (io_byte)(file_bytes>>8);
  magic[4] = (io_byte)(file_bytes>>16);
  magic[5] = (io_byte)(file_bytes>>24);
  magic[6] = magic[7] = magic[8] = magic[9] = 0;
  magic[10] = (io_byte) header_bytes;
  magic[11] = (io_byte)(header_bytes>>8);
  magic[12] = (io_byte)(header_bytes>>16);
  magic[13] = (io_byte)(header_bytes>>24);
  header.size = 40;
  header.width = width;
  header.height = height;
  header.planes_bits = 1; // Set `planes'=1 (mandatory)
  header.planes_bits |= ((num_components==1)?8:24) << 16; // Set bits per pel.
  header.compression = 0;
  header.image_size = 0;
  header.xpels_per_metre = header.ypels_per_metre = 0;
  header.num_colours_used = header.num_colours_important = 0;
  to_little_endian((io_int32 *) &header,10);
  if ((state->out = fopen(fname,"wb")) == NULL)
    return(IO_ERR_NO_FILE);

  fwrite(magic,1,14,state->out);
  fwrite(&header,1,40,state->out);
  if (num_components == 1)
    for (int n=0; n < 256; n++)
      { fputc(n,state->out); fputc(n,state->out);
        fputc(n,state->out); fputc(0,state->out); }
  return 0;
}

/*****************************************************************************/
/*                               bmp_out__close                              */
/*****************************************************************************/

void bmp_out__close(bmp_out *state)
{
  if (state->out != NULL)
    fclose(state->out);
  memset(state,0,sizeof(bmp_out));
}

/*****************************************************************************/
/*                              bmp_out__put_line                            */
/*****************************************************************************/

int bmp_out__put_line(bmp_out *state, io_byte *line)
{
  if ((state->out == NULL) || (state->num_unwritten_rows <= 0))
    return(IO_ERR_FILE_NOT_OPEN);
  state->num_unwritten_rows--;
  if (fwrite(line,1,(size_t) state->line_bytes,state->out) !=
      (size_t) state->line_bytes)
    throw(IO_ERR_FILE_TRUNC);
  if (state->alignment_bytes > 0)
    {
      io_byte buf[3] = {0,0,0};
        fwrite(buf,1,(size_t) state->alignment_bytes,state->out);
    }
  return 0;
}
