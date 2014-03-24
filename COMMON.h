#ifndef COMMON_H
#define COMMON_H

/********************************************************************************************/
/********************************************************************************************/
#define NDIM 100
/* Nombre de fonctions J0(alpha_i*rho) utilisees dans le developpement en serie de Bessel. */

/********************************************************************************************/
#define MASSE_PROTON 0.938272013
#undef  MASSE_PROTON
#define MASSE_PROTON 0.938
/* La masse du proton est exprimee en [GeV]. */

#define E_PROTON_MIN (7.*MASSE_PROTON)
/* E_PROTON_MIN designe l'energie TOTALE MINIMALE des protons. Nous prenons pour l'instant
le seuil de production des antiprotons lors d'une collision proton + proton -----> pbar + X.
ATTENTION : quand nous considererons l'helium, ce seuil est susceptible de s'abaisser.
L'energie est exprimee en [GeV]. */

#define E_PROTON_MAX (1.e6)
/* E_PROTON_MAX designe l'energie TOTALE MAXIMALE des protons. Elle est exprimee en [GeV]. */

#define DIM_TAB_PROTON 500
/* DIM_TAB_PROTON est le nombre d'intervalles en energie TOTALE de proton.
ATTENTION : pour pouvoir realiser une integration a la SIMPSON, il est imperatif
d'avoir ici un nombre PAIR. */
#if ((DIM_TAB_PROTON % 2) == 1)
  #error DIM_TAB_PROTON DOIT ETRE PAIR = LE MODIFIER EN CONSEQUENCE !
#endif

/********************************************************************************************/
#define T_PBAR_MIN (0.1)
/* T_PBAR_MIN est la valeur MINIMALE de l'energie CINETIQUE des antiprotons
que l'on considere dans le probleme. Elle est exprimee en [GeV]. */

#define T_PBAR_MAX (10000.)
/* T_PBAR_MAX est la valeur MAXIMALE de l'energie CINETIQUE des antiprotons
que l'on considere dans le probleme. Elle est exprimee en [GeV]. */

#define DIM_TAB_PBAR 250
/* DIM_TAB_PBAR est le nombre d'intervalles en energie CINETIQUE des antiprotons.
ATTENTION : pour pouvoir realiser une integration a la SIMPSON, il est imperatif
d'avoir ici un nombre PAIR. */
#if ((DIM_TAB_PBAR % 2) == 1)
  #error DIM_TAB_PBAR DOIT ETRE PAIR = LE MODIFIER EN CONSEQUENCE !
#endif

/********************************************************************************************/
/*
* Liste des fichiers et de leurs adresses ou sont stockees les valeurs de
* la section efficace differentielle de production des antiprotons
* \frac{d \sigma_{\pbar}}{d E_{\pbar}} au cours de la reaction generique
*
* {P ou ALPHA} + {H ou HE}_{milieu interstellaire au repos} -----> PBAR + X
*
* L'energie CINETIQUE de l'antiproton produit au cours de cette reaction est
* notee T_pbar et varie de T_PBAR_MIN a T_PBAR_MAX en prenant (DIM_TAB_PBAR + 1)
* valeurs differentes.
*
* L'energie TOTALE E_proton des protons incidents OU ALORS l'energie TOTALE
* PAR NUCLEON E_nucleon des particules ALPHA incidentes -- noyaux d'helium --
* varie de E_PROTON_MIN a E_PROTON_MAX en prenant (DIM_TAB_PROTON + 1) valeurs
* differentes.
*
* - FILE_NAME_H_ON_H   contient le tableau DSPBAR_SUR_DEPBAR_H_ON_H  [i_pbar][i_proton]
* - FILE_NAME_H_ON_HE  contient le tableau DSPBAR_SUR_DEPBAR_H_ON_HE [i_pbar][i_proton]
* - FILE_NAME_HE_ON_H  contient le tableau DSPBAR_SUR_DEPBAR_HE_ON_H [i_pbar][i_nucleon]
* - FILE_NAME_HE_ON_HE contient le tableau DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon]
*
* Le nombre entier i_pbar varie de 0 a DIM_TAB_PBAR.
* Le nombre entier i_proton OU ALORS i_nucleon varie de 0 a DIM_TAB_PROTON.
*
*/

#define FILE_NAME_H_ON_H   "DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included"
#define FILE_NAME_H_ON_HE  "DSPBAR_SUR_DEPBAR_H_on_HE_high_energy"
#define FILE_NAME_HE_ON_H  "DSPBAR_SUR_DEPBAR_HE_on_H_high_energy"
#define FILE_NAME_HE_ON_HE "DSPBAR_SUR_DEPBAR_HE_on_HE_high_energy"

/********************************************************************************************/
#define E_DISC      0.1
/* La demi-epaisseur du disque mince galactique est exprimee en [kpc]. */

#define R_GAL       20.0
/* Le rayon du disque galactique est exprime en [kpc]. */

#define R_EARTH      8.5
/* Le rayon galactocentrique du systeme solaire est exprime en [kpc]. */

#define DENSITE_H_DISC 0.9
/* La densite d'hydrogene neutre diffus dans le disque mince est exprimee en [cm^{-3}]. */

#define DENSITE_HE_DISC 0.1
/* La densite d'helium neutre diffus dans le disque mince est exprimee en [cm^{-3}]. */

/********************************************************************************************/
#define PI (acos(-1.))
//#define PI 3.14159265358979

#define CELERITY_LIGHT 2.99792458e10
/* La vitesse de la lumiere est exprimee en [cm s^{-1}]. */

#define CM_PAR_PARSEC 3.0856775807e+18
/* Il s'agit du nombre de centimetres correspondant a une distance
 d'un parsec [pc] */

#define CM_PAR_KPC 3.0856775807e+21
/* Il s'agit du nombre de centimetres correspondant a une distance
 d'un kiloparsec [kpc]. */

#define SEC_PAR_MGYR (3.15581498e13)
/* Il s'agit du nombre de secondes correspondant a un laps de
 temps d'un million d'annees -- soit 1 megayear ou encore 1 MGYR -- [sec]. */

#define SEC_PAR_KYR (3.15581498e10)
/* Il s'agit du nombre de secondes correspondant a un laps de temps
 d'un millier d'annees -- soit 1 kiloyear ou encore 1 KYR -- [sec]. */

#define SEC_PAR_YR (3.15581498e7)
/* Il s'agit du nombre de secondes correspondant a un laps de temps
 d'une annee [sec]. */

#define GEV_PER_KEV 1.0e-6
#define KEV_PER_GEV 1.0e+6
#define ERG_PER_GEV 1.60217657e-3

/********************************************************************************************/

//	EINASTO PROFILE PARAMETERS

#define rhos_Ein	0.033			// [Gev cm{-3}]
#define rs_Ein		28.44			// [kpc]
#define alpha_Ein	0.17			// [NO UNIT]

/********************************************************************************************/

//	PARAMETRES DES ANTIPROTONS PRIMAIRES

//	Choix du cannal d'annihilation :

//	1	:	eL
//	2	:	eR
//	3	:	muL
//	4	:	muR
//	5	:	tauL
//	6	:	tauR
//	7	:	q
//	8	:	c
//	9	:	b
//	10	:	t
//	11	:	WL
//	12	:	WT
//	13	:	ZL
//	14	:	ZT
//	15	:	g
//	16	:	gamma

#define number_channels     16
#define channel_choice		 1


//	Choix de la masse du WIMP et de sa section efficace d'annihilation

#define mass_chi_choice		 10.0					//	[GeV]
#define sigma_v_annihilation 3.0e-26				// [cm^{3} s^{-1}]

/********************************************************************************************/

#define N_x_pbar_scan	300
#define x_pbar_scan_min	1.0e-9
#define x_pbar_scan_max	1.0

#define N_gaelle_masses	62

/********************************************************************************************/

//	For SOLAR MINIMUM, the Fisk potential is PHI_FISK_MIN =  500 MV = 0.5 GV.
//	For SOLAR MAXIMUM, the Fisk potential is PHI_FISK_MAX = 1000 MV = 1.0 GV.

#define fisk_potential 0.0

/********************************************************************************************/

#define antiproton_spectrum_file_name	"./results/primary_antiproton_spectrum_test" 

/********************************************************************************************/
/********************************************************************************************/
#endif