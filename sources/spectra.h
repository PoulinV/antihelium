#ifndef SPECTRA_H
#define SPECTRA_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "../COMMON.h"
#include "STRUCTURES.h"

#include "BESSEL_PRELIM.h"
#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"
#include "SOLAR_MOD.h"
#include "PROTON.h"
#include "HELIUM.h"
#include "PRIMARY_PBAR.h"
#include "ANTI_PROTON.h"


/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

void print_total_pbar_spectra_MIN_MED_MAX(struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium, struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation, struct Structure_Primary_Source_Term* pt_Primary_Source_Term, double alpha_i[NDIM+1]);

void print_PBAR_IS_SPECTRUM(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1]);
void print_PBAR_TOA_SPECTRUM(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double T_TOA[DIM_TAB_PBAR+1]);


/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif