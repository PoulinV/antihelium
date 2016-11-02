#ifndef COMMON_H
#define COMMON_H

/********************************************************************************************/
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

#define MASSE_PROTON 0.938272013
#undef  MASSE_PROTON
#define MASSE_PROTON 0.938
/* La masse du proton est exprimee en [GeV]. */

/********************************************************************************************/

//	GALACTIC PARAMETERS

#define E_DISC      0.1
/* La demi-epaisseur du disque mince galactique est exprimee en [kpc]. */

#define R_GAL       20.0
/* Le rayon du disque galactique est exprime en [kpc]. */

//#define R_EARTH      8.5
#define R_EARTH      8.33
/* Le rayon galactocentrique du systeme solaire est exprime en [kpc]. */

#define DENSITE_H_DISC 0.9
/* La densite d'hydrogene neutre diffus dans le disque mince est exprimee en [cm^{-3}]. */

#define DENSITE_HE_DISC 0.1
/* La densite d'helium neutre diffus dans le disque mince est exprimee en [cm^{-3}]. */

/********************************************************************************************/

//	PROPAGATION PARAMETERS

//#define MIN
#define MED
//#define MAX

//	Nombre de parametres de propagation
#define nParamProp  5

//	Nombre de jeux de parametres de propagation
#define nJeuxParam  1623 
//#define nJeuxParam  10 


//	For SOLAR MINIMUM, the Fisk potential is PHI_FISK_MIN =  500 MV = 0.5 GV.
//	For SOLAR MAXIMUM, the Fisk potential is PHI_FISK_MAX = 1000 MV = 1.0 GV.

#define fisk_potential 0.724								//	[GV]
//#define fisk_potential 0.0								//	[GV]
//#define fisk_potential 0.93									//	[GV]


/********************************************************************************************/

//	PBAR SPECTRUM TABLE PARAMETERS

#define T_PBAR_MIN (0.1)
//#define T_PBAR_MIN (0.01)

//#define T_PBAR_MIN (0.2)
/* T_PBAR_MIN est la valeur MINIMALE de l'energie CINETIQUE des antiprotons
que l'on considere dans le probleme. Elle est exprimee en [GeV]. */

//#define T_PBAR_MAX (1.e4)
#define T_PBAR_MAX (1.e5)
/* T_PBAR_MAX est la valeur MAXIMALE de l'energie CINETIQUE des antiprotons
que l'on considere dans le probleme. Elle est exprimee en [GeV]. */

//#define DIM_TAB_PBAR 250
#define DIM_TAB_PBAR 300
//#define DIM_TAB_PBAR 350

/* DIM_TAB_PBAR est le nombre d'intervalles en energie CINETIQUE des antiprotons.
ATTENTION : pour pouvoir realiser une integration a la SIMPSON, il est imperatif
d'avoir ici un nombre PAIR. */
#if ((DIM_TAB_PBAR % 2) == 1)
  #error DIM_TAB_PBAR DOIT ETRE PAIR = LE MODIFIER EN CONSEQUENCE !
#endif

/********************************************************************************************/

//	PROTON SPECTRUM TABLE PARAMETERS

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

// On definit l'intervalle d'energie cinetique des protons pour afficher leur spectre.
#define DIM_TAB_PROTON_SPECTRUM	DIM_TAB_PBAR
#define T_PROTON_SPECTRUM_MIN	T_PBAR_MIN						// [GeV] 
#define T_PROTON_SPECTRUM_MAX	T_PBAR_MAX						// [GeV]

/********************************************************************************************/

//	BESSEL PARAMETERS

#define NDIM 100
//#define NDIM 170

/* Nombre de fonctions J0(alpha_i*rho) utilisees dans le developpement en serie de Bessel. */

/********************************************************************************************/

//	EXPERIMENTAL PROTON FLUX PARAMETRIZATION

//	On choisit ici l'experience qui nous donne le flux de protons.

//	New BESS data from Shikaze et al. VERSION AS OF 080903.
//#define BESS_2008_proton_Shikaze

//	New parametrization from Fiorenza and David fits to H data. VERSION AS OF 081023.
//#define Fit_2008_proton_Maurin_Donato

//	New parameterization proposed by Julien Lavalle and based on the CREAM high energy CR proton data. The F1p fit is published in arXiv:1011.3063.
//#define CREAM_2010_proton_Lavalle

//	New parameterization proposed by Julien Lavalle and based on the ATIC_2 high energy CR proton data. The F2p fit is published in arXiv:1011.3063.
//#define ATIC2_2010_proton_Lavalle

//   New parameterization proposed by Timur Delahaye and based on the PAMELA high energy CR proton data.
#define PAMELA_2012_proton_Delahaye

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Fiorenza Donato in arXiv:1402.0321.
//#define AMS02_2013_proton_Donato

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Kappl and Winkler in arXiv:1408.0299.
//#define AMS02_2013_proton_Kappl_Winkler

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Vittino in a forthcoming paper.
//#define AMS02_2013_proton_Vittino

//	Fit performed by Manuela Vecchi on data from a forthcoming paper of AMS-02.
//#define AMS02_2015_proton_manuela
	
//	Fit performed by Yoann on data from a forthcoming paper of AMS-02.
//#define AMS02_2015_proton_yoann	

/********************************************************************************************/

//	EXPERIMENTAL HELIUM FLUX PARAMETRIZATION

//	On choisit ici l'experience qui nous donne le flux d'helium.

//	New BESS data from Shikaze et al. VERSION AS OF 080903.
//#define BESS_2008_helium_Shikaze

//	New parametrization from Fiorenza and David fits to H data. VERSION AS OF 081023.
//#define Fit_2008_helium_Maurin_Donato

//	New parameterization proposed by Julien Lavalle and based on the CREAM and ATIC_2 high energy CR helium data. The F1He fit is published in arXiv:1011.3063.
//#define CREAM_ATIC2_2010_helium_Lavalle

//   New parameterization proposed by Timur Delahaye and based on the PAMELA high energy CR helium data.
#define PAMELA_2012_helium_Delahaye

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Fiorenza Donato in arXiv:1402.0321.
//#define AMS02_2013_helium_Donato

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Kappl and Winkler in arXiv:1408.0299.
//#define AMS02_2013_helium_Kappl_Winkler

//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Vittino in a forthcoming paper.
//#define AMS02_2013_helium_Vittino

//	Fit performed by Yoann on data from AMS Days
//#define AMS02_2015_helium_yoann



/********************************************************************************************/

//	SECONDARY PBAR PRODUCTION CROSS-SECTION PARAMETERS

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


// Ajouter Tan_Ng_theta_only

//#define DSPBAR_SUR_DEPBAR_H_on_H_unknown
//#define DSPBAR_SUR_DEPBAR_H_on_H_Duperray
//#define DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included
#define DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12

#ifdef DSPBAR_SUR_DEPBAR_H_on_H_unknown
	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_unknown.txt"
#elif defined DSPBAR_SUR_DEPBAR_H_on_H_Duperray
	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_Duperray.txt"
#elif defined DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included
	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included.txt"
#elif defined DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12.txt"
#else
	#error You must to specify one parametrization for H on H reaction!
#endif



#define FILE_NAME_H_ON_HE	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_HE_high_energy.txt"
#define FILE_NAME_HE_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_HE_on_H_high_energy.txt"
#define FILE_NAME_HE_ON_HE	"../sources/cross_section/DSPBAR_SUR_DEPBAR_HE_on_HE_high_energy.txt"

/********************************************************************************************/

//	PRIMARY PBAR PARAMETERS

#define RHO_CHI_SOLAR 0.3
// La densite de masse des neutralinos dans le voisinage solaire est exprimee en [GeV cm^{-3}]. 
#define RHO_CHI_0 1.0
// La valeur de reference pour la densite de masse des neutralinos est exprimee en [GeV cm^{-3}]. 
#define RC_SMBH 0.1
// Renormalization radius expressed in [kpc]


//	Dark matter profile

//#define NFW
//#define moore
#define einasto
//#define einastoB
//#define isothermal
//#define burkert


#define WIMP_annihilation
//#define WIMP_decay

#define number_channels     23

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
//	17	:	h
//	18	:	nue
//	19	:	numu
//	20	:	nutau
//	21	:	V->e
//	22	:	V->mu
//	23	:	V->tau

#define channel_choice	9


//	Choix de la masse du WIMP et de sa section efficace d'annihilation

#define mass_chi_choice		 1000.0					//	[GeV]

//	Choix de sa section efficace d'annihilation et de son taux de desintegration

#define sigma_v_annihilation 3.0e-26				//	[cm^{3} s^{-1}]
#define decay_rate 1.0e-26							//	[s^{-1}]


#define N_x_pbar_scan	300
#define x_pbar_scan_min	1.0e-9
#define x_pbar_scan_max	1.0

#define N_gaelle_masses	62


#define mass_chi_inf  5.0
#define mass_chi_sup  5000.0

#define N_mass_chi 2

/********************************************************************************************/

//	RESULTS FILE NAMES

//	On definit ici le nom des fichiers dans lesquels on stocke les resultats.

#define proton_IS_spectrum_file_name	"./results/proton_IS_spectrum.txt"
#define proton_TOA_spectrum_file_name	"./results/proton_TOA_spectrum.txt"
#define proton_exp_spectrum_file_name	"./results/proton_exp_spectrum.txt"

#define helium_exp_spectrum_file_name	"./results/helium_exp_spectrum.txt"

#define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum.txt"
#define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum.txt" 

#define pbar_over_p_IS_spectrum_file_name	"./results/pbar_over_p_IS_spectrum.txt"
#define pbar_over_p_TOA_spectrum_file_name	"./results/pbar_over_p_TOA_spectrum.txt" 

#define pbar_IS_spectra_MIN_MED_MAX_file_name	"./results/pbar_IS_spectra_MIN_MED_MAX.txt"
#define pbar_TOA_spectra_MIN_MED_MAX_file_name	"./results/pbar_TOA_spectra_MIN_MED_MAX.txt" 

#define pbar_over_p_IS_uncertainty_spectrum_file_name	"./results/pbar_over_p_IS_uncertainty_spectrum.txt"
#define pbar_over_p_TOA_uncertainty_spectrum_file_name	"./results/pbar_over_p_TOA_uncertainty_spectrum.txt" 


/********************************************************************************************/
/********************************************************************************************/
#endif