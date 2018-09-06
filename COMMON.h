#ifndef COMMON_H
#define COMMON_H

/********************************************************************************************/
/********************************************************************************************/
#define _TRUE_ 1
#define _FALSE_ 0
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

#define MASS_NEUTRON 939.565
#define MASS_HELIUM_TROIS 2809.413
/* La masse du proton est exprimee en [GeV]. */

#define SIGMA_DBARP_NOANN 4e-27     //in cm^2
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
// #define MAX

// #define BREAK_IN_DIFFUSION

#ifdef BREAK_IN_DIFFUSION
  #define _DIFFUSION_COEFF_RIGIDITY_BREAK_ 200. // [GV] Defines the rigidity break of the diffusion coefficient K. If set negative or nul, K will have no break.
  #define BREAK_PUISSANCE_COEFF_DIFF 0.33 // Defines the value of delta after the rigidity break of the diffusion coefficient K.
#else
  #define _DIFFUSION_COEFF_RIGIDITY_BREAK_ 0. // [GV] Defines the rigidity break of the diffusion coefficient K. If set negative or nul, K will have no break.
  #define BREAK_PUISSANCE_COEFF_DIFF 0. // Defines the value of delta after the rigidity break of the diffusion coefficient K.
#endif
//	Nombre de parametres de propagation
#define nParamProp  5

//	Nombre de jeux de parametres de propagation
#define nJeuxParam  1623
//#define nJeuxParam  10


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
// #define DIM_TAB_PBAR 300
#define DIM_TAB_PBAR 350

/* DIM_TAB_PBAR est le nombre d'intervalles en energie CINETIQUE des antiprotons.
ATTENTION : pour pouvoir realiser une integration a la SIMPSON, il est imperatif
d'avoir ici un nombre PAIR. */
#if ((DIM_TAB_PBAR % 2) == 1)
  #error DIM_TAB_PBAR DOIT ETRE PAIR = LE MODIFIER EN CONSEQUENCE !
#endif

// /*********************WHICH NUCLEI SHOULD BE COMPUTED?*******************************************************/

// #define  ANTIPROTON
// #define  ANTIDEUTERIUM
// #define  ANTIHELIUM3
#define  ANTIHELIUM4

///PUT TO FALSE ONCE FILES HAVE BEEN WRITTEN: ACCELERATES GREATLY CALCULATION ///
// #define WRITE_CROSS_SECTION _TRUE_
#define WRITE_CROSS_SECTION _FALSE_
//

//
#define TERTIARY_COMPUTATION _TRUE_
// #define TERTIARY_COMPUTATION _FALSE_

//AN ALTERNATIVE SCHEME FOR COALESCENCE, PUT TO FALSE TO FOLLOW 1808.08961 (Poulin++)
// #define A_LA_CHARDONNET _TRUE_
#define A_LA_CHARDONNET _FALSE_



// /**************************************************************************/
//	PROTON SPECTRUM TABLE PARAMETERS


#ifdef ANTIPROTON
  #define E_PROTON_MIN (7.*MASSE_PROTON) //for He3
  #define M_NUCLEI 0.938
  #define A_NUCLEI 1
  #define Z_NUCLEI 1
  #define BA_COAL_FID 1 //USELESS
  #define P_COALESCENCE 0.261 //USELESS

#endif
#ifdef ANTIDEUTERIUM
  #define E_PROTON_MIN (17.*MASSE_PROTON) //for He3
  #define M_NUCLEI 1.876
  #define A_NUCLEI 2
  #define Z_NUCLEI 1
  #define BA_COAL_FID 0.016 //From 1709.08522 fig 10
  #define BA_COAL 0.021 //From 1709.08522 fig 10
  #define P_COALESCENCE 0.261 //GeV
#endif
#ifdef ANTIHELIUM3
  #define E_PROTON_MIN (31.*MASSE_PROTON) //for He3
  #define M_NUCLEI 2.809
  #define A_NUCLEI 3
  #define Z_NUCLEI 2
  #define BA_COAL_FID 2e-4 //From 1709.08522 fig 11 right panel
  #define BA_COAL 3e-4  //From 1709.08522 fig 10
  #define P_COALESCENCE 0.261 //GeV
#endif
#ifdef ANTIHELIUM4
  #define E_PROTON_MIN (49.*MASSE_PROTON) //for He3
  #define M_NUCLEI 3.730
  #define A_NUCLEI 4
  #define Z_NUCLEI 2
  #define BA_COAL_FID 0.0 //NOT YET MEASURED!
  #define BA_COAL     0.0 //NOT YET MEASURED!
  #define P_COALESCENCE 0.261 //GeV
#endif
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
/********************************************************************************************/

// Valeur moyenne pour AMS-02 (2011-2013) donne par Ghelfi et al. 2015
#define Fisk_Average
//#define Fisk_Plus_1_Sigma
//#define Fisk_Plus_2_Sigma
//#define Fisk_Plus_3_Sigma
//#define Fisk_Minus_1_Sigma
//#define Fisk_Minus_2_Sigma
//#define Fisk_Minus_3_Sigma

#ifdef Fisk_Average
#define phi_fisk      0.724		//[GV]
//#define phi_fisk      0.0		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_724MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_724MV
#endif

#ifdef Fisk_Plus_1_Sigma
#define phi_fisk      0.755		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_755MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_755MV
#endif

#ifdef Fisk_Plus_2_Sigma
#define phi_fisk      0.783		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_783MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_783MV
#endif

#ifdef Fisk_Plus_3_Sigma
#define phi_fisk      0.830		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_830MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_830MV
#endif

#ifdef Fisk_Minus_1_Sigma
#define phi_fisk      0.701		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_701MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_701MV
#endif

#ifdef Fisk_Minus_2_Sigma
#define phi_fisk      0.671		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_671MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_671MV
#endif

#ifdef Fisk_Minus_3_Sigma
#define phi_fisk      0.647		//[GV]
#define AMS_CREAM_BCV_2015_H_phi_fisk_647MV
#define AMS_CREAM_BCV_2015_He_phi_fisk_647MV
#endif


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
//#define DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
#define DSPBAR_SUR_DEPBAR_H_on_H_RW

// #ifdef ANTIPROTON
  #ifdef DSPBAR_SUR_DEPBAR_H_on_H_unknown
  	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_unknown.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_Duperray
  	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_Duperray.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included
  	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
  	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_RW
  	#define FILE_NAME_H_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_H_RW.txt"
  #else
  	#error You must to specify one parametrization for H on H reaction!
  #endif
// #endif

// #ifdef ANTIPROTON
  #define FILE_NAME_H_ON_HE	"../sources/cross_section/DSPBAR_SUR_DEPBAR_H_on_HE_high_energy.txt"
  #define FILE_NAME_HE_ON_H	"../sources/cross_section/DSPBAR_SUR_DEPBAR_HE_on_H_high_energy.txt"
  #define FILE_NAME_HE_ON_HE	"../sources/cross_section/DSPBAR_SUR_DEPBAR_HE_on_HE_high_energy.txt"
// #endif
//
// #ifdef ANTIDEUTERIUM
  #if (A_LA_CHARDONNET == _TRUE_)
    #define FILE_NAME_DE_H_ON_H "../sources/cross_section/DSDe_SUR_DEDe_H_on_H_MDGS_F12_vChardonnet.txt"
  #else
    #define FILE_NAME_DE_H_ON_H "../sources/cross_section/DSDe_SUR_DEDe_H_on_H_MDGS_F12.txt"
  #endif
  #define FILE_NAME_DE_H_ON_HE "../sources/cross_section/DSDe_SUR_DEDe_H_on_HE_MDGS_F12.txt"
  #define FILE_NAME_DE_HE_ON_H "../sources/cross_section/DSDe_SUR_DEDe_HE_on_H_MDGS_F12.txt"
  #define FILE_NAME_DE_HE_ON_HE "../sources/cross_section/DSDe_SUR_DEDe_HE_on_HE_MDGS_F12.txt"
// #endif
//
  #ifdef DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
    #if (A_LA_CHARDONNET == _TRUE_)
      #define FILE_NAME_HE3_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_H_on_H_MDGS_F12_vChardonnet.txt"
    #else
      #define FILE_NAME_HE3_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_H_on_H_MDGS_F12.txt"
    #endif
    #define FILE_NAME_HE3_H_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_H_on_HE_MDGS_F12.txt"
    #define FILE_NAME_HE3_HE_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_HE_on_H_MDGS_F12.txt"
    #define FILE_NAME_HE3_HE_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_HE_on_HE_MDGS_F12.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_RW
    #define FILE_NAME_HE3_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_H_on_H_RW.txt"
    #define FILE_NAME_HE3_H_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_H_on_HE_RW.txt"
    #define FILE_NAME_HE3_HE_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_HE_on_H_RW.txt"
    #define FILE_NAME_HE3_HE_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE3BAR_HE_on_HE_RW.txt"
  #endif
  #ifdef DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
    #if (A_LA_CHARDONNET == _TRUE_)
      #define FILE_NAME_HE4_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_H_on_H_MDGS_F12_vChardonnet.txt"
    #else
      #define FILE_NAME_HE4_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_H_on_H_MDGS_F12.txt"
    #endif
    #define FILE_NAME_HE4_H_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_H_on_HE_MDGS_F12.txt"
    #define FILE_NAME_HE4_HE_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_HE_on_H_MDGS_F12.txt"
    #define FILE_NAME_HE4_HE_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_HE_on_HE_MDGS_F12.txt"
  #elif defined DSPBAR_SUR_DEPBAR_H_on_H_RW
    #define FILE_NAME_HE4_H_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_H_on_H_RW.txt"
    #define FILE_NAME_HE4_H_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_H_on_HE_RW.txt"
    #define FILE_NAME_HE4_HE_ON_H "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_HE_on_H_RW.txt"
    #define FILE_NAME_HE4_HE_ON_HE "../sources/cross_section/DSHEBAR_SUR_DEHE4BAR_HE_on_HE_RW.txt"
  #endif
// #endif
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

#if (TERTIARY_COMPUTATION == _TRUE_)
  #if (A_LA_CHARDONNET == _FALSE_)
    #ifdef MED
      #ifdef BREAK_IN_DIFFUSION
      #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_full_MED_BREAK.txt"
      #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_full_MED_BREAK.txt"
      #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_full_MED_BREAK.txt"
      #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_full_MED_BREAK.txt"
      #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_full_MED_BREAK.txt"
      #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_full_MED_BREAK.txt"
      #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_full_MED_BREAK.txt"
      #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_full_MED_BREAK.txt"
      #else
      #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_full_MED.txt"
      #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_full_MED.txt"
      #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_full_MED.txt"
      #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_full_MED.txt"
      #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_full_MED.txt"
      #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_full_MED.txt"
      #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_full_MED.txt"
      #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_full_MED.txt"
      #endif
    #endif
    #ifdef MAX
      #ifdef BREAK_IN_DIFFUSION
      #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_full_MAX_BREAK.txt"
      #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_full_MAX_BREAK.txt"
      #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_full_MAX_BREAK.txt"
      #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_full_MAX_BREAK.txt"
      #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_full_MAX_BREAK.txt"
      #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_full_MAX_BREAK.txt"
      #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_full_MAX_BREAK.txt"
      #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_full_MAX_BREAK.txt"
      #else
      #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_full_MAX.txt"
      #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_full_MAX.txt"
      #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_full_MAX.txt"
      #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_full_MAX.txt"
      #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_full_MAX.txt"
      #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_full_MAX.txt"
      #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_full_MAX.txt"
      #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_full_MAX.txt"
      #endif
    #endif

  #elif (A_LA_CHARDONNET == _TRUE_)
    #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_full_chardonnet.txt"
    #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_full_chardonnet.txt"
    #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_full_chardonnet.txt"
    #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_full_chardonnet.txt"
    #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_full_chardonnet.txt"
    #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_full_chardonnet.txt"
    #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_full_chardonnet.txt"
    #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_full_chardonnet.txt"
  #endif
#endif
// #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_ppOnly.txt"
// #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_ppOnly.txt"
// #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_ppOnly.txt"
// #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_ppOnly.txt"
// #define Hebar_IS_spectrum_file_name	"./results/Hebar_IS_spectrum_ppOnly.txt"
// #define Hebar_TOA_spectrum_file_name	"./results/Hebar_TOA_spectrum_ppOnly.txt"
#if (TERTIARY_COMPUTATION == _FALSE_)
  #if (A_LA_CHARDONNET == _FALSE_)
      #ifdef MED
        #ifdef BREAK_IN_DIFFUSION
          #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_noTertiary_MED_BREAK.txt"
          #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_noTertiary_MED_BREAK.txt"
          #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_noTertiary_MED_BREAK.txt"
          #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_noTertiary_MED_BREAK.txt"
          #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_noTertiary_MED_BREAK.txt"
          #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_noTertiary_MED_BREAK.txt"
          #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_noTertiary_MED_BREAK.txt"
          #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_noTertiary_MED_BREAK.txt"
        #else
          #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_noTertiary_MED.txt"
          #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_noTertiary_MED.txt"
          #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_noTertiary_MED.txt"
          #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_noTertiary_MED.txt"
          #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_noTertiary_MED.txt"
          #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_noTertiary_MED.txt"
          #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_noTertiary_MED.txt"
          #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_noTertiary_MED.txt"
        #endif
      #endif
      #ifdef MAX
        #ifdef BREAK_IN_DIFFUSION
          #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_noTertiary_MAX_BREAK.txt"
          #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_noTertiary_MAX_BREAK.txt"
          #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_noTertiary_MAX_BREAK.txt"
          #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_noTertiary_MAX_BREAK.txt"
          #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_noTertiary_MAX_BREAK.txt"
          #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_noTertiary_MAX_BREAK.txt"
          #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_noTertiary_MAX_BREAK.txt"
          #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_noTertiary_MAX_BREAK.txt"
        #else
          #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_noTertiary_MAX.txt"
          #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_noTertiary_MAX.txt"
          #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_noTertiary_MAX.txt"
          #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_noTertiary_MAX.txt"
          #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_noTertiary_MAX.txt"
          #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_noTertiary_MAX.txt"
          #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_noTertiary_MAX.txt"
          #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_noTertiary_MAX.txt"
        #endif
      #endif

    #elif (A_LA_CHARDONNET == _TRUE_)
      #define pbar_IS_spectrum_file_name	"./results/pbar_IS_spectrum_noTertiary_chardonnet.txt"
      #define pbar_TOA_spectrum_file_name	"./results/pbar_TOA_spectrum_noTertiary_chardonnet.txt"
      #define Debar_IS_spectrum_file_name	"./results/Debar_IS_spectrum_noTertiary_chardonnet.txt"
      #define Debar_TOA_spectrum_file_name	"./results/Debar_TOA_spectrum_noTertiary_chardonnet.txt"
      #define He3bar_IS_spectrum_file_name	"./results/He3bar_IS_spectrum_noTertiary_chardonnet.txt"
      #define He4bar_IS_spectrum_file_name	"./results/He4bar_IS_spectrum_noTertiary_chardonnet.txt"
      #define He3bar_TOA_spectrum_file_name	"./results/He3bar_TOA_spectrum_noTertiary_chardonnet.txt"
      #define He4bar_TOA_spectrum_file_name	"./results/He4bar_TOA_spectrum_noTertiary_chardonnet.txt"
    #endif
#endif

#define pbar_over_p_IS_spectrum_file_name	"./results/pbar_over_p_IS_spectrum.txt"
#define pbar_over_p_TOA_spectrum_file_name	"./results/pbar_over_p_TOA_spectrum.txt"
#define Debar_over_p_IS_spectrum_file_name	"./results/Debar_over_p_IS_spectrum.txt"
#define Debar_over_p_TOA_spectrum_file_name	"./results/Debar_over_p_TOA_spectrum.txt"
#define He3bar_over_p_IS_spectrum_file_name	"./results/He3bar_over_p_IS_spectrum.txt"
#define He4bar_over_p_IS_spectrum_file_name	"./results/He4bar_over_p_IS_spectrum.txt"
#define He3bar_over_p_TOA_spectrum_file_name	"./results/He3bar_over_p_TOA_spectrum.txt"
#define He4bar_over_p_TOA_spectrum_file_name	"./results/He4bar_over_p_TOA_spectrum.txt"

#define pbar_IS_spectra_MIN_MED_MAX_file_name	"./results/pbar_IS_spectra_MIN_MED_MAX.txt"
#define pbar_TOA_spectra_MIN_MED_MAX_file_name	"./results/pbar_TOA_spectra_MIN_MED_MAX.txt"

#define pbar_over_p_IS_uncertainty_spectrum_file_name	"./results/pbar_over_p_IS_uncertainty_spectrum.txt"
#define pbar_over_p_TOA_uncertainty_spectrum_file_name	"./results/pbar_over_p_TOA_uncertainty_spectrum.txt"
#ifdef ANTIPROTON
    #define secondary_source_term_file_name "./results/pbar_H_on_H_source_term.txt"
#endif
#ifdef ANTIDEUTERIUM
  #if (A_LA_CHARDONNET == _TRUE_)
    #define secondary_source_term_file_name "./results/Debar_H_on_H_source_term_chardonnet.txt"
  #else
    #define  secondary_source_term_file_name "./results/Debar_H_on_H_source_term.txt"
  #endif
#endif
#ifdef ANTIHELIUM3
  #if (A_LA_CHARDONNET == _TRUE_)
    #define secondary_source_term_file_name "./results/HE3bar_H_on_H_source_term_chardonnet.txt"
  #else
      #define secondary_source_term_file_name "./results/HE3bar_H_on_H_source_term.txt"
  #endif
#endif
#ifdef ANTIHELIUM4
  #if (A_LA_CHARDONNET == _TRUE_)
    #define secondary_source_term_file_name "./results/HE4bar_H_on_H_source_term_chardonnet.txt"
  #else
      #define secondary_source_term_file_name "./results/HE4bar_H_on_H_source_term.txt"
  #endif
#endif
/********************************************************************************************/
/********************************************************************************************/
#endif
