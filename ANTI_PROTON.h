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

#include "COMMON.h"
#include "STRUCTURES.h"

#include "DAVID_CROSS_SECTIONS_torsten.h"
#include "DIFFUSION_PROPAGATION_julien_clump_1012.h"
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

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#undef NRANSI

#endif