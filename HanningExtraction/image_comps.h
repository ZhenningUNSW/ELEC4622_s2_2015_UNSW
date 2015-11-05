/*****************************************************************************/
// File: image_comps.h
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

/*****************************************************************************/
/* STRUCT                        my_image_comp                               */
/*****************************************************************************/

struct my_image_comp {
    // Data members: (these occupy space in the structure's block of memory)
    int width;
    int height;
    int stride;
    int border; // Extra rows/cols to leave around the boundary
    float *handle; // Points to start of allocated memory buffer
    float *buf; // Points to the first real image sample
    // Function members: (these do not occupy any space in memory)
    my_image_comp()
      { width = height = stride = border = 0;  handle = buf = NULL; }
    ~my_image_comp()
      { if (handle != NULL) delete[] handle; }
    void init(int height, int width, int border)
      {
        this->width = width;  this->height = height;  this->border = border;
        stride = width + 2*border;
        if (handle != NULL)
          delete[] handle; // Delete mem allocated by any previous `init' call
        handle = new float[stride*(height+2*border)];
        buf = handle + (border*stride) + border;
      }
    void perform_boundary_extension();
       // This function is implemented in "filtering_main.cpp".
  };
  /* Notes:
       * This "struct" is actually a C++ class.  For some of you, this may
         be your first introduction to classes, but don't worry.
       * C++ classes can be declared using either "struct" or "class"; for
         the moment, you don't need to worry about the subtle differences.
         We are using "struct" here to remind you that there is very little
         difference between a C structure and a C++ class.  To the computer,
         both are blocks of memory containing variables.
       * The main extra thing that comes with classes is the ability to
         associate the data with functions (or rather to associate functions
         with the data).  As we are using the functions here, they are just
         the same as any other C function, except that the syntax for invoking
         them involves the instance of the structure that they are to work on.
       * Every function declared as part of a class has an implicit argument
         (in the machine code it is the first argument of the real function)
         named "this".  This first argument is a pointer to the instance of
         the structure on which the function was invoked.
       * For example, if "comp" is a variable of type "my_image_comp", the
         statement "comp.init(256,256,4);" essentially translates into
         something like "my_image_comp__init(&comp,256,256,4);".  There is
         no magic here.  It's just what we have been doing in the "io_bmp"
         code, except that the syntax is a little less clumsy.
       * There are only two other things you need to know:
         a) Inside the implementation of a class function, you can refer to
            a member variable via the "this" pointer or just by name -- we
            do both in the present sample.
         b) Every class has two special functions (even if they are not
            explicitly implemented), the constructor and destructor.  The
            constructor function has the same name as the class and no
            return type (not even void).  It is invoked when any instance
            of the class is instantiated, either by declaring it as a
            variable or allocating it with "new".  The destructor has the
            name of the class, with a prepended tilde, and again no return
            type.  It is invoked when an instance of the class goes out of
            scope (e.g., the function in which it was declared as a variable
            returns) or an instance allocated with "new" is deleted with
            "delete".  If an array of class instances is allocated with
            new, each member's constructor is invoked, and each member's
            destructor is later invoked when "delete" is used to destroy
            the array.  This is why the special form "delete[]" needs to
            be used for arrays -- otherwise, the run-time system will not
            know to access the hidden information left behind by "new",
            indicating how many objects there are to destroy.
  */
