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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::ApplyEOSTable(const float* TotalMetals)
{

  const float Zsolar = 0.0204;

  /* If TotalMetals == NULL, then consider zero metallicity case */

  bool MetalsExist = (TotalMetals != NULL);

  /* TotalMetals is metal density in code units, use GetUnits to
     convert to cgs. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
    
  /* Find fields: density, total energy, velocity1-3. */

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum, B1Num, B2Num, B3Num);

  /* Find Multi-species fields. */

  if (MultiSpecies)
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			  HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  /* Loop over grid.  Check which cells have rho >
     MinimumDensityForEOSTable, and then replace the cell temperature
     with the EOS temperature */

  int i, dim, size;
  for (dim = 0; size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  int dindex, zindex, dindex1, zindex1;
  float dfactor, zfactor, drho, dz, drho_inv, dz_inv, rho_thres, z0, rho0, Z_H, 
    log_rho, dom, new_temperature, mu;
  rho_thres = MinimumDensityForEOSTable * mh / DensityUnits;
  drho = (CoolData.EOS_DensityRange[1] - CoolData.EOS_DensityRange[0]) / 
    (CoolData.EOS_NumberOfDensityBins - 1);
  dz = (CoolData.EOS_MetallicityRange[1] - CoolData.EOS_MetallicityRange[0]) / 
    (CoolData.EOS_NumberOfMetallicityBins - 1);
  drho_inv = 1.0 / drho;
  dz_inv = 1.0 / dz;
  dom = DensityUnits / mh;

  for (i = 0; i < size; i++) {
    
    // If density is greater than threshold value, modify thermal energy.
    if (BaryonField[DensNum][i] > rho_thres) {

      // Calculate interpolation values for density
      log_rho = log10(BaryonField[DensNum][i] * dom);
      log_rho = min(max(log_rho, CoolData.EOS_DensityRange[0]), CoolData.EOS_DensityRange[1]);
      dindex = int((log_rho - CoolData.EOS_DensityRange[0]) * drho_inv);
      rho0 = CoolData.EOS_DensityRange[0] + drho * dindex;
      dfactor = (log_rho - rho0) * drho_inv;

      // Calculate interpolation values for metallicity
      if (MetalsExist) {
	Z_H = log10(TotalMetals[i] / BaryonField[DensNum][i] / Zsolar);
	Z_H = min(max(Z_H, CoolData.EOS_MetallicityRange[0]), CoolData.EOS_MetallicityRange[1]);
	zindex = int((Z_H - CoolData.EOS_MetallicityRange[0]) * dz_inv);
	z0 = CoolData.EOS_MetallicityRange[0] + dz * zindex;
	zfactor = (Z_H - z0) * dz_inv;
      } else {
	zindex = 0;  // lowest (hopefully nearly zero!) metallicity
	zfactor = 0.0;
      }

      dindex1 = min(dindex + 1, CoolData.EOS_NumberOfDensityBins-1);
      zindex1 = min(zindex + 1, CoolData.EOS_NumberOfMetallicityBins-1);
      new_temperature = 
	CoolData.EOS_Table[zindex ][dindex ] * (1 - zfactor) * (1 - dfactor) +
	CoolData.EOS_Table[zindex1][dindex ] * (zfactor)     * (1 - dfactor) +
	CoolData.EOS_Table[zindex ][dindex1] * (1 - zfactor) * (dfactor) +
	CoolData.EOS_Table[zindex1][dindex1] * (zfactor)     * (dfactor);

      if (MultiSpecies == 0) {
	mu = Mu;
      } else {
	if (MultiSpecies > 0)
	  mu = BaryonField[DeNum][i] + BaryonField[HINum][i] + BaryonField[HIINum] + 
	    4.0 * (BaryonField[HeINum][i] + BaryonField[HeIINum][i] + BaryonField[HeIIINum][i]);
	if (MultiSpecies > 1)
	  mu += BaryonField[HMNum][i] + 2.0 * (BaryonField[H2INum][i] + BaryonField[H2IINum][i]);
	if (MultiSpecies > 2)
	  mu += BaryonField[DINum][i] + BaryonField[DIINum][i] + 2.0 * BaryonField[HDINum][i];
	mu /= BaryonField[DensNum][i];
      } // ENDELSE
      
      
    } // ENDIF overdense

  } // ENDFOR cells

  return SUCCESS;

}
