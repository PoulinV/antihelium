#ifndef ANTI_PROTON_H
#define ANTI_PROTON_H

#define NRANSI

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "../COMMON.h"
#include "STRUCTURES.h"

#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"
#include "SOLAR_MOD.h"

#include "TRIDAG.h"
#include "nrutil.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

void calculation_BESSEL_PBAR_SECONDARY_Epbar_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium, struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TERTIARY_Epbar_i(double alpha_i[NDIM+1], struct Structure_Pbar* pt_Pbar);
void calculation_BESSEL_PBAR_SUM_123_Epbar_i(struct Structure_Pbar* pt_Pbar);

void calculation_BESSEL_PBAR_TOT_direct_inversion_A(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TOT_direct_inversion_B(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TOT_direct_inversion_TD_NR(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(double time_max,long n_time,struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PBAR_TOT_diffusion_soluce_B(double time_max,long n_time,struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation);

void inversion_tridiagonal(double a[DIM_TAB_PBAR+1],double b[DIM_TAB_PBAR+1],
double c[DIM_TAB_PBAR+1],double r[DIM_TAB_PBAR+1],double u[DIM_TAB_PBAR+1]);

extern void gaussj(float **a, int n, float **b, int m);

void PBAR_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PBAR+1]);
void PBAR_BESSEL_TABLES_123_initialization(struct Structure_Pbar* pt_Pbar);

void PBAR_IS_SPECTRUM_calculation(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1]);
void PBAR_TOA_SPECTRUM_calculation(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double PBAR_TOA_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_TOA[DIM_TAB_PBAR+1], struct Structure_Propagation* pt_Propagation);

void tertiary_component_effect_calculation(struct Structure_Pbar* pt_Pbar, double alpha_i[NDIM+1]);
void ELDR_effect_calculation(struct Structure_Propagation* pt_Propagation, struct Structure_Pbar* pt_Pbar, double alpha_i[NDIM+1]);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#undef NRANSI

#endif