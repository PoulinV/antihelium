#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

/*
* CE MODULE RENVOIE EN MILLIBARN LA SECTION EFFICACE DES PROTONS, DES ANTIPROTONS ET
* DES ANTIDEUTERONS SUR L'HYDROGENE, SECTION EFFICACE PARAMETREE PAR LES DONNEES DU CERN
* (http://pdg.lbl.gov)
*
* THE PARAMETERIZATION IS
* sigma(p) = A + B*p^n + C*ln^2(p) + D*ln(p)
* WHERE SIGMA IS IN MILLIBARN AND p IS IN GEV/C. THE BEST-FIT COEFFICIENTS A, B, C AND D,
* AND THE EXPONENT n ARE TABULATED BELOW. ALSO GIVEN IS THE RANGE OF MOMENTUM OVER WHICH
* THE FIT WAS DONE.
*/

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "../COMMON.h"
#include "STRUCTURES.h"

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES CONSTANTES

/*
* VALABLE POUR L'INTERVALLE: 3.0-2100 GEV/C
*/
#define p_seuil_pp_tot 2.83
#define A_pp_tot 48.0 	    /* erreur=+/- 0.1   */
#define B_pp_tot 0.0
#define n_pp_tot 0.0
#define C_pp_tot 0.522			/* erreur=+/- 0.005 */
#define D_pp_tot -4.51			/* erreur=+/- 0.05  */

/*
* VALABLE POUR L'INTERVALLE: 5.0-1.73x10^6 GEV/C
*/
#define p_seuil_pbarp_tot 4.733
#define A_pbarp_tot 38.4		/* erreur=+/- 4.4   */
#define B_pbarp_tot 77.6		/* erreur=+/- 2.8   */
#define n_pbarp_tot -0.64		/* erreur=+/- 0.07  */
#define C_pbarp_tot 0.26		/* erreur=+/- 0.05  */
#define D_pbarp_tot -1.2		/* erreur=+/- 0.9   */

/*
* VALABLE POUR L'INTERVALLE: 5.0-1.73x10^6 GEV/C
*/
#define p_seuil_pbarp_el 5.029
#define A_pbarp_el 10.2			/* erreur=+/- 0.7   */
#define B_pbarp_el 52.7			/* erreur=+/- 1.8   */
#define n_pbarp_el -1.16		/* erreur=+/- 0.05  */
#define C_pbarp_el 0.125		/* erreur=+/- 0.014 */
#define D_pbarp_el -1.28		/* erreur=+/- 0.2   */

/*
* VALABLE POUR L'INTERVALLE: 9.5-280 GEV/C
*/
#define p_seuil_max_pbard_tot 9.5
#define A_pbard_tot 112			/* erreur=+/- 13    */
#define B_pbard_tot 125			/* erreur=+/- 8     */
#define n_pbard_tot -1.08		/* erreur=+/- 0.15  */
#define C_pbard_tot 1.14 		/* erreur=+/- 0.49  */
#define D_pbard_tot -12.4		/* erreur=+/- 4.9   */


/*
* POUR LES SECTIONS EFFICACES LORSQUE L'IMPULSION EST INFERIEURE A CELLE DONNEE
* CI-DESSUS, J'AJUSTE LES DONNEES DU CERN AVEC UNE FONCTION POLYNOMIALE EN ln(p).
* LES DETAILS DES 'OUTPUTS' DE CES AJUSTEMENTS PEUVENT ETRES TROUVES DANS LE
* FICHIER fit_sigma_p_et_pbar.doc DANS LE MEME REPERTOIRE.
*
* sigma = P1 + P2*ln(p) + P3*ln^2(p) + P4*ln^3(p) + P5*ln^4(p) +
*         P6*ln^5(p) + P7*ln^6(p)
*/

/*
* VALABLE POUR L'INTERVALLE: 0.14-3.0 GEV/C
*/ 
#define  P1_pp_tot 29.038		  /* erreur=+/- 0.62113E-01 */
#define  P2_pp_tot 37.126		  /* erreur=+/- 0.16886	    */
#define  P3_pp_tot 28.260		  /* erreur=+/- 0.26031	    */
#define  P4_pp_tot -52.708		/* erreur=+/- 0.30967 	  */
#define  P5_pp_tot -28.658		/* erreur=+/- 0.30778	    */
#define  P6_pp_tot 16.918		  /* erreur=+/- 0.16909	    */
#define  P7_pp_tot 14.320		  /* erreur=+/- 0.14092 	  */

/*
* VALABLE POUR L'INTERVALLE: 0.18-5.0 GEV/C
*/ 
#define  P1_pbarp_tot 117.44		/* erreur=+/- 0.15108 	*/
#define  P2_pbarp_tot -49.027		/* erreur=+/- 0.24424 	*/
#define  P3_pbarp_tot 27.982		/* erreur=+/- 0.41893 	*/
#define  P4_pbarp_tot -21.626		/* erreur=+/- 0.51670  	*/
#define  P5_pbarp_tot -0.18134	/* erreur=+/- 0.26325 	*/
#define  P6_pbarp_tot 4.3295		/* erreur=+/- 0.23804 	*/
#define  P7_pbarp_tot 0.0

/*
* VALABLE POUR L'INTERVALLE: 0.15-5.0 GEV/C
*/
#define  P1_pbarp_el 45.091		  /* erreur=+/- 0.23099 	*/
#define  P2_pbarp_el -14.050		/* erreur=+/- 0.43102 	*/
#define  P3_pbarp_el 0.56917		/* erreur=+/- 0.54940 	*/
#define  P4_pbarp_el -12.832		/* erreur=+/- 0.55406 	*/
#define  P5_pbarp_el 4.9824		  /* erreur=+/- 0.41239 	*/
#define  P6_pbarp_el 3.7379		  /* erreur=+/- 0.22140 	*/
#define  P7_pbarp_el -1.5916		/* erreur=+/- 0.11591 	*/

/*
* Pour l'interaction antiproton sur cible de deuterium ou antideuterium
* sur proton, l'approximation devient
*
* sigma = P1 + P2*p + P3*p^2 + P4*p^3 + P5*p^4 + P6*p^5 + P7*p^6 +
*         P8*p^7 + P9*p^8
*
* VALABLE POUR L'INTERVALLE: 0.3-3.3 GEV/C
*/
#define  P1_pbard_min_tot 1208.7
#define  P2_pbard_min_tot -4552.7
#define  P3_pbard_min_tot 9040.5
#define  P4_pbard_min_tot -9949.2
#define  P5_pbard_min_tot 6513.0
#define  P6_pbard_min_tot -2585.9
#define  P7_pbard_min_tot 606.73
#define  P8_pbard_min_tot -76.699
#define  P9_pbard_min_tot 3.9769

/*
* VALABLE POUR L'INTERVALLE: 2.-9.5 GEV/C
*/
#define p_seuil_min_pbard_tot 2.
#define  P1_pbard_max_tot  5.30429
#define  P2_pbard_max_tot -0.31054


#define mb_cm2 (1.e-27)         /* Constante de conversion mb en cm2 */
#define uma_Gev (931.48e-3)     /* Constante de conversion uma en Gev */
#define m_p 1.00794
#define MASSE_DEUT 1.876        /* La masse du deuteron est exprimee en [GeV] */

/*
* CODE DE LAURENT DEROME
*/
#define PMASS 0.93827231
//in units of [GeV]

/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/

//	DECLARATION DES FONCTIONS

double sigma_total_pH(double E_proton);
double sigma_inelastic_pH_TAN_and_NG(double E_proton);

double sigma_total_pbarH(double E_pbar);
double sigma_elastic_pbarH(double E_pbar);
double sigma_inelastic_pbarH_TAN_and_NG(double E_pbar);
double sigma_inelastic_NOANN_pbarH_TAN_and_NG(double E_pbar);

double sigma_total_dbarH(double E_dbar);

double E_d3S_on_d3P_PBAR_ECM(double E_CMF,double pLstar,double pTstar);
double E_d3S_on_d3P_PBAR_LAB(double E_proton,double pL,double pT,
double *E_CMF,double *pLstar,double *pTstar);
double dSpbar_sur_dEpbar_DIRECTE(double E_proton,double E_pbar,long n_step_theta);
double dSpbar_sur_dEpbar_SIMPSON(double E_proton,double E_pbar,long n_step_theta);

void DSPBAR_SUR_DEPBAR_H_ON_H_write_file(struct Structure_Cross_Section* pt_Cross_Section);
void CLEANING_ALL_THE_DSPBAR_SUR_DEPBAR(struct Structure_Cross_Section* pt_Cross_Section);
void DSPBAR_SUR_DEPBAR_H_ON_H_read_file(struct Structure_Cross_Section* pt_Cross_Section);
void DSPBAR_SUR_DEPBAR_H_ON_HE_read_file(struct Structure_Cross_Section* pt_Cross_Section);
void DSPBAR_SUR_DEPBAR_HE_ON_H_read_file(struct Structure_Cross_Section* pt_Cross_Section);
void DSPBAR_SUR_DEPBAR_HE_ON_HE_read_file(struct Structure_Cross_Section* pt_Cross_Section);

double GetInvarMul_pAapX(double At,double p1,double y,double mt);
double GetReactionCrossSection(double At,double Pproj);

double E_d3S_on_d3P_PBAR_H_ON_HE_LAB(double E_proton,double pL,double pT);
double dSpbar_sur_dEpbar_SIMPSON_H_ON_HE(double E_proton,double E_pbar,long n_step);
void   DSPBAR_SUR_DEPBAR_H_ON_HE_write_file(struct Structure_Cross_Section* pt_Cross_Section);

double E_d3S_on_d3P_PBAR_HE_ON_H_LAB(double E_nucleon,double pL,double pT);
double dSpbar_sur_dEpbar_SIMPSON_HE_ON_H(double E_nucleon,double E_pbar,long n_step);
void   DSPBAR_SUR_DEPBAR_HE_ON_H_write_file(struct Structure_Cross_Section* pt_Cross_Section);

void DSPBAR_SUR_DEPBAR_HE_ON_HE_write_file(struct Structure_Cross_Section* pt_Cross_Section);

double invariant_multiplicity_pH_apX(double p1,double y,double mt);
double E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(double E_proton,double pL,double pT);
double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_duperray(double E_proton,double E_pbar,long n_step);

double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_tan_ng_mass_T(double E_proton,double E_pbar,long n_step);

double sigma_total_pH_MDGS(double E_proton);
double sigma_elastic_pH_MDGS(double E_proton);
double sigma_inelastic_pH_MDGS(double E_proton);
double invariant_multiplicity_pH_apX_MDGS_F12(double p1,double y,double mt);
double E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(double E_proton,double pL,double pT);
double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_MDGS_F12(double E_proton,double E_pbar,long n_step);


/**********************************************************************************************************************************************************************************************************/
/**********************************************************************************************************************************************************************************************************/
#endif