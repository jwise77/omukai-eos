/***********************************************************************
/
/  READ A EQUATION OF STATE TABLE (OMUKAI COOLING)
/
/  written by: Daisuke Yamasawa & John Wise
/  date:       June, 2013
/  modified1:  
/
/  PURPOSE:
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

int InitializeEquilibriumCoolDataEOS(const char *table_name)
{

  char *delims = (char*) "\t";
  char *value;
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int i, j, count, index;
  float dummy;

  fptr = fopen(table_name, "r");
  if(fptr == NULL){
    ENZO_VFAIL("Error opening %s", table_name);
  }

  index = 0;
  while((fgets(line, MAX_LINE_LENGTH, fptr)) != NULL){
    sscanf(line, "# Density Bins = %d", &CoolData.EOS_NumberOfDensityBins);
    if (sscanf(line, "# Metallicity Bins = %d", &CoolData.EOS_NumberOfMetallicityBins) > 0) {
      CoolData.EOS_Table = new float*[CoolData.EOS_NumberOfMetallicityBins];
      for (i = 0; i < CoolData.EOS_NumberOfMetallicityBins; i++)
	CoolData.EOS_Table[i] = new float[CoolData.EOS_NumberOfDensityBins];
    }
    sscanf(line, "# Density Range = %"FSYM" %"FSYM, &CoolData.EOS_DensityRange[0], &CoolData.EOS_DensityRange[1]);
    sscanf(line, "# Metallicity Range = %"FSYM" %"FSYM, &CoolData.EOS_MetallicityRange[0], &CoolData.EOS_MetallicityRange[1]);
	  
    if(line[0] != '#'){
      value = strtok(line, delims);
      dummy = atof(value);
      value = strtok(NULL, delims);
      count = 0;
      for (j = 0; j < CoolData.EOS_NumberOfMetallicityBins; j++) {
	CoolData.EOS_Table[j][index] = atof(value);
	value = strtok(NULL, delims);
      }
      index++;
    }
  }
  fclose(fptr);

  return SUCCESS;

}
