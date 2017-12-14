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

void print_PBAR_IS_SPECTRUM(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double A_nuclei);
void print_PBAR_TOA_SPECTRUM(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_TOA[DIM_TAB_PBAR+1], double A_nuclei);

void print_PROTON_IS_SPECTRUM(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1]);
void print_PROTON_TOA_SPECTRUM(double PROTON_TOA_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double T_PROTON_TOA[DIM_TAB_PROTON_SPECTRUM+1]);

void print_PROTON_SPECTRUM_exp(void);

void print_PBAR_OVER_P_IS_SPECTRUM(double PBAR_OVER_P_IS_SPECTRUM[DIM_TAB_PBAR+1]);
void print_PBAR_OVER_P_TOA_SPECTRUM(double PBAR_OVER_P_TOA_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_OVER_P_TOA[DIM_TAB_PBAR+1]);

void print_PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY(double PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2]);
void print_PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY(double PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2], double T_PBAR_OVER_P_TOA[DIM_TAB_PBAR+1]);

void print_HELIUM_SPECTRUM_exp(void);

void print_secondary_source_term(double PRIMARY_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Cross_Section* pt_Cross_Section, double TARGET_DENSITY,  double A_nuclei);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif
