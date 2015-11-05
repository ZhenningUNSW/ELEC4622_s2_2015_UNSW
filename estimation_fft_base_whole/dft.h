/*****************************************************************************/
// File: dft.h
// Author: David Taubman
// Last Revised: 28 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

/*****************************************************************************/
/* CLASS                         my_direct_dft                               */
/*****************************************************************************/

class my_direct_dft {
  public: // Publically available member functions of the class.
    my_direct_dft()
      { N=0; real_buf=imag_buf=NULL; real_trig=imag_trig=NULL; }
      /* After construction, be sure to call `init', before using the
         `forward_transform' or `inverse_transform' functions. */
    ~my_direct_dft() { cleanup(); }
    void init(int N, bool is_forward);
      /* You must call this function before using the `perform_transform'
         function, to set up internal buffers and trigonometric tables.  You
         should try to call this function as infrequently as possible -- i.e.,
         try to use the transform functions multiple times between calls to
         this function, if that is appropriate for the application.  This is
         because setting up the trigonometric tables can be time consuming.
         The `N' argument is the length of the DFT.  The `is_forward' argument
         is true if you want `perform_transform' to perform a forward DFT;
         otherwise an inverse transform will be performed.  You may call the
         function multiple times to change the value of N.
       */
    void perform_transform(float *real, float *imag, int stride);
      /* This function performs a 1D DFT (or inverse DFT), taking its source
         data from the supplied `real' and `imag' arrays and writing the result
         back to these same arrays.  The number of elements to be processed
         by the DFT/IDFT is given by the `N' value supplied in the most
         recent call to `init'.  The transform is a forward DFT if the
         `forward' argument to `init' was true; otherwise an inverse DFT is
         performed.  Successive samples in the `real' and `imag'
         arrays are separated by `stride'.  This allows for a variety of
         different buffer organizations. */
  private: // Internal helper functions
    void cleanup()
      {
        if (real_buf != NULL) delete[] real_buf;
        if (imag_buf != NULL) delete[] imag_buf;
        if (real_trig != NULL) delete[] real_trig;
        if (imag_trig != NULL) delete[] imag_trig;
        real_buf = imag_buf = NULL;  real_trig = imag_trig = NULL;
        N = 0;
      }
  private: // Data members
    int N;
    float *real_buf; // Working buffer with `N' entries
    float *imag_buf; // Working buffer with `N' entries
    double *real_trig; // Holds cos(2*pi*n/N) for n=0,...,N-1
    double *imag_trig; // For inverse transforms, this array holds sin(2*pi*n/N)
                       // For forward transforms, it holds sin(-2*pi*n/N).
  };


class my_fft {
public: // Publically available member functions of the class.
	my_fft() {
		N = 0; buf = NULL; trig = NULL;
	}
	/* After construction, be sure to call `init', before using the
	`forward_transform' or `inverse_transform' functions. */
	~my_fft() { cleanup(); }
	void init(int N, bool is_forward);
	/* You must call this function before using the `perform_transform'
	function, to set up internal buffers and trigonometric tables.  You
	should try to call this function as infrequently as possible -- i.e.,
	try to use the transform functions multiple times between calls to
	this function, if that is appropriate for the application.  This is
	because setting up the trigonometric tables can be time consuming.
	The `N' argument is the length of the DFT.  The `is_forward' argument
	is true if you want `perform_transform' to perform a forward DFT;
	otherwise an inverse transform will be performed.  You may call the
	function multiple times to change the value of N.
	*/
	_complex* perform_transform(_complex *buf, int N);
	/* This function performs a 1D DFT (or inverse DFT), taking its source
	data from the supplied `real' and `imag' arrays and writing the result
	back to these same arrays.  The number of elements to be processed
	by the DFT/IDFT is given by the `N' value supplied in the most
	recent call to `init'.  The transform is a forward DFT if the
	`forward' argument to `init' was true; otherwise an inverse DFT is
	performed.  Successive samples in the `real' and `imag'
	arrays are separated by `stride'.  This allows for a variety of
	different buffer organizations. */
private: // Internal helper functions
	void cleanup()
	{
		if (buf != NULL) delete[] buf;
		if (trig != NULL) delete[] trig;
		buf = NULL;
		trig = NULL;
		N = 0;
	}
private: // Data members
	int N;
	_complex *buf;
	_complex *trig;
};