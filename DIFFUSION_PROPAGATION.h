#ifndef DIFFUSION_PROPAGATION_H
#define DIFFUSION_PROPAGATION_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "besselj0_next.h"
#include "besselj1_next.h"

#include "COMMON.h"
#include "STRUCTURES.h"



/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES CONSTANTES

#define K_BOLTZMANN (8.617385e-14)
/* La constante de BOLTZMANN est exprime en [GeV Kelvin^{-1}]. */

#define H_BAR (6.58211915e-25)
/* La constante de PLANCK est exprime en [GeV s]. */

#define MASSE_ELECTRON 0.511e-03
/* La masse de l'lectron est exprime en [GeV]. */

#define RADIUS_ELECTRON 2.818e-13
/* Le rayon classique de l'lectron est exprim en [cm]. */

#define V_ION_H  19.e-9
#define V_ION_HE 44.e-9
/* Les potentiels d'ionisation de l'hydrogne et de l'hlium sont
respectivement nots V_ION_H et V_ION_HE. Ils sont exprims en [GeV]. */

#define DENSITE_FREE_ELECTRON 0.033
/* La densit d'lectrons libres dans le milieu interstellaire -- ISM -- est
exprime en [cm^{-3}]. */

#define T_ELECTRONIC 3.e5
/* La temprature du plasma lectronique est exprime en [Kelvin] -- voir
Manheim & al. */

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

double K_space_diffusion(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation);
double D_energy_diffusion(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation);
double b_energy_losses(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation);
double b_energy_losses_for_wei_calculations(double energy,double mass,double Z_em,double *ionisation,double *coulomb,double *adiabatic,double *drift, struct Structure_Propagation* pt_Propagation);

double GENERIC_FLUX(double r,double z,
double energy,double mass,double Z_em,double alpha_i[NDIM+1],double BESSEL_COEFFICIENTi[NDIM+1], struct Structure_Propagation* pt_Propagation);

double GENERIC_FLUX_04(double r,double z,
double energy,double mass,double Z_em,double alpha_i[NDIM+1],double BESSEL_COEFFICIENTi[NDIM+1], struct Structure_Propagation* pt_Propagation);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif