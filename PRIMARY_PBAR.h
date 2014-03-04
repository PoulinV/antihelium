#ifndef PRIMARY_PBAR_H
#define PRIMARY_PBAR_H

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

#define RHO_CHI_SOLAR 0.3
/* La densite de masse des neutralinos dans le voisinage solaire est exprimee en [GeV cm^{-3}]. */
#define RHO_CHI_0 1.0
/* La valeur de reference pour la densite de masse des neutralinos est exprimee en [GeV cm^{-3}]. */

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

void   calculation_BESSEL_PBAR_PRIMARY_Epbar_i(long n_vert, long n_rad, double alpha_i[NDIM+1],
                                               struct Structure_Pbar* pt_Pbar,
                                               struct Structure_Propagation* pt_Propagation,
                                               struct Structure_Primary_Source_Term* pt_Primary_Source_Term);

double rapport_rho_chi_sur_rho_0              (double rr,double z);
double rapport_rho_chi_sur_rho_0_Einasto      (double rr,double z);

void   DNPBAR_ON_DTPBAR_gaelle_read_file      (struct Structure_Primary_Source_Term* pt_Primary_Source_Term);
double dNpbar_on_dEpbar_primary_calculation   (double mass_chi, int channel, struct Structure_Primary_Source_Term* pt_Primary_Source_Term);
void   primary_source_calculation             (double mass_chi, struct Structure_Primary_Source_Term* pt_Primary_Source_Term);

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif