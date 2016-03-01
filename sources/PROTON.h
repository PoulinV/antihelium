#ifndef PROTON_H
#define PROTON_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "../COMMON.h"
#include "STRUCTURES.h"

#include "besselj0_next.h"
#include "besselj1_next.h"
#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"
#include "SOLAR_MOD.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

double flux_proton_EXP(double E_proton);

double Q_proton_tot(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation);
void calcul_method_A_BESSEL_Pi(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation);
void calcul_method_B_BESSEL_Pi(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_PROTON_Ep_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation);

double flux_proton_TH(double r,double z,double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation);
double pbar_emissivity_per_H_solar(double E_pbar);

void PROTON_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1]);

void PROTON_IS_SPECTRUM_calculation(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1]);
void PROTON_TOA_SPECTRUM_calculation(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double PROTON_TOA_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double T_PROTON_TOA[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Propagation* pt_Propagation);

double fit_proton_flux_AMS02_manuela(double R);
double fit_proton_flux_AMS02_manuela_demodulated(double E_proton);
double fit_proton_flux_AMS02_manuela_demodulated_2(double EnIS);
void FFA_IS_to_TOA_modified(double A,double Z,double PHI,double EnIS,double flux_IS,double *EnTOA,double *flux_TOA);

double fit_proton_flux_AMS02_yoann_demodulated_2(double EnIS);
double fit_proton_flux_AMS02_yoann(double R);



/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif