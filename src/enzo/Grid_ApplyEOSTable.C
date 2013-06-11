/***********************************************************************
/
/  USE AN EQUATION OF STATE TO OVERRIDE TEMPERATURE AT HIGH DENSITY
/
/  written by: Daisuke Yamasawa & John Wise
/  date:       June, 2013
/  modified1:  
/
/  NOTES:   The equation of state is a function of density and 
/           metallicity.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

int grid::ApplyEOSTable(void)
{

  

  return SUCCESS;

}
