#ifndef HELIUM_0809_H
#define HELIUM_0809_H

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

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES CONSTANTES

#define A_HE 4.0
/* Le nombre de nuclons dans un noyau d'hlium est not A_HE. */

#define Z_HE 2.0
/* La charge lectrique d'un noyau d'hlium est not Z_HE. */

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

double flux_helium_EXP(double E_nucleon);

void calcul_method_B_BESSEL_HEi(double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_HELIUM_Ep_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);

double flux_helium_TH(double r,double z,double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);

extern double besselj0(double x);
extern double besselj1(double x);

extern double sigma_total_pH(double E_proton);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif