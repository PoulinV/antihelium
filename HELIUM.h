#ifndef HELIUM_H
#define HELIUM_H

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
#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES CONSTANTES

#define A_HE 4.0
/* Le nombre de nucleons dans un noyau d'helium est note A_HE. */

#define Z_HE 2.0
/* La charge electrique d'un noyau d'helium est note Z_HE. */

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

double flux_helium_EXP(double E_nucleon);

void calcul_method_B_BESSEL_HEi(double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);
void calculation_BESSEL_HELIUM_Ep_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);

double flux_helium_TH(double r,double z,double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif