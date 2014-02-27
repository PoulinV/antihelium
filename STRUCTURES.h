#ifndef STRUCTURES_H
#define STRUCTURES_H

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "COMMON.h"
/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

// DECLARATION DES STRUCTURES


struct Structure_Propagation
{
	double PUISSANCE_COEFF_DIFF;
	double DIFFUSION_0_GV;
	double E_DIFFUS;
	double VENT_GALACTIQUE;
	double V_ALFEN;
	double PHI_FISK;
};


struct Structure_Nuclei
{
	double q_i[NDIM+1];
	double BESSEL_COEF_i[NDIM+1];
	double BESSEL_COEF_Enuc_i[DIM_TAB_PROTON+1][NDIM+1];
};


struct Structure_Cross_Section
{
	double DSPBAR_SUR_DEPBAR_H_ON_H  [DIM_TAB_PBAR+1][DIM_TAB_PROTON+1];
	double DSPBAR_SUR_DEPBAR_H_ON_HE [DIM_TAB_PBAR+1][DIM_TAB_PROTON+1];
	double DSPBAR_SUR_DEPBAR_HE_ON_H [DIM_TAB_PBAR+1][DIM_TAB_PROTON+1];
	double DSPBAR_SUR_DEPBAR_HE_ON_HE[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1];
};


struct Structure_Primary_Source_Term
{

	double X_PBAR[N_x_pbar_scan+1];

	double DNPBAR_ON_DTPBAR_CHANNEL[number_channels+1][N_x_pbar_scan+1];

	double DNPBAR_ON_DEPBAR   [DIM_TAB_PBAR+1];
	double PRIMARY_SOURCE_TERM[DIM_TAB_PBAR+1];
};


struct Structure_Pbar
{
	double BESSEL_PBARi[NDIM+1];
	double TABLE_Abar_i[DIM_TAB_PBAR+1][NDIM+1];
	double BESSEL_PBAR_PRI_Epbar_i[DIM_TAB_PBAR+1][NDIM+1];
	double BESSEL_PBAR_SEC_Epbar_i[DIM_TAB_PBAR+1][NDIM+1];
	double BESSEL_PBAR_TER_Epbar_i[DIM_TAB_PBAR+1][NDIM+1];
	double BESSEL_PBAR_TOT_Epbar_i[DIM_TAB_PBAR+1][NDIM+1];	
};


/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/
#endif	