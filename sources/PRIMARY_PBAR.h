#ifndef PRIMARY_PBAR_H
#define PRIMARY_PBAR_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include <sys/types.h> 
#include <sys/stat.h> 

#include "../COMMON.h"
#include "STRUCTURES.h"

#include "besselj0_next.h"
#include "besselj1_next.h"
#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"
#include "ANTI_PROTON.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES CONSTANTES

//#define RHO_CHI_SOLAR 0.3
// La densite de masse des neutralinos dans le voisinage solaire est exprimee en [GeV cm^{-3}]. 
//#define RHO_CHI_0 1.0
// La valeur de reference pour la densite de masse des neutralinos est exprimee en [GeV cm^{-3}]. 

//#define RC_SMBH 0.1
// Renormalization radius expressed in [kpc]

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

void   calculation_BESSEL_PBAR_PRIMARY_Epbar_i(long n_vert, long n_rad, double alpha_i[NDIM+1],
                                               struct Structure_Pbar* pt_Pbar,
                                               struct Structure_Propagation* pt_Propagation,
                                               struct Structure_Primary_Source_Term* pt_Primary_Source_Term);

double rapport_rho_chi_sur_rho_0_old          (double rr,double z);
double rapport_rho_chi_sur_rho_0              (double rr,double z);

void   DNPBAR_ON_DTPBAR_gaelle_read_file      (double mass_chi, struct Structure_Primary_Source_Term* pt_Primary_Source_Term);
double dNpbar_on_dEpbar_primary_calculation   (double mass_chi, int channel, struct Structure_Primary_Source_Term* pt_Primary_Source_Term);
void   primary_source_calculation             (double mass_chi, struct Structure_Primary_Source_Term* pt_Primary_Source_Term);

void DM_preliminary(struct Structure_Primary_Source_Term* pt_Primary_Source_Term);

void primary_spectra_BCGS_2014(struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation, struct Structure_Primary_Source_Term* pt_Primary_Source_Term, double alpha_i[NDIM+1]);

void DM_source_term_calculation(struct Structure_Primary_Source_Term* pt_Primary_Source_Term);


/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif