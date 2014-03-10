#ifndef PROTON_0809_H
#define PROTON_0809_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "COMMON.h"
#include "STRUCTURES.h"

#include "besselj0_next.h"
#include "besselj1_next.h"
#include "DIFFUSION_PROPAGATION_julien_clump_1012.h"
#include "DAVID_CROSS_SECTIONS_torsten.h"

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

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif