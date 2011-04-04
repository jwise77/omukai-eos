// $Id: test_Parse.cpp 1970 2011-02-01 20:22:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-06
/// @brief    Test program for reading in parameters then displaying them

#include <stdio.h>
#include <string.h>
#include "test.hpp"
#include "parameters.hpp"

const int debug = true;

int main()
{
  cello_parameters_read(stdin);
  cello_parameters_print();
}
