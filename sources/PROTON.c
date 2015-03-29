#include "PROTON.h"

/********************************************************************************************/
/********************************************************************************************/
/*
*  Il s'agit du flux interstellaire des protons. Ce flux est un flux differentiel par
*  rapport a l'energie. Il est exprime en unite de [cm^{-2} s^{-1} sr^{-1} GeV^{-1}].
*
*  Several options are possible.
*  The GAISSER and SCHAEFER proton spectrum is
*  resultat = NORMALISATION_FLUX_PROTONS * pow(impulsion,-DELTA) /
*  sqrt(1.+pow(P0/impulsion,2));
*
*  The BEREZINSKY proton spectrum is
*  resultat = 1.8 * (impulsion/E) * pow(E,-2.75);
*  with minimal spectral index of
*  resultat = 1.8 * (impulsion/E) * pow(E,-2.70);
*  and maximal spectral index of
*  resultat = 1.8 * (impulsion/E) * pow(E,-2.80);
*
*  A possible value for the fit is
*  resultat = N * beta_lorentz^{-1} * impulsion^{-gamma}
*  resultat = N * (E/impulsion) * pow(impulsion,-gamma)
*  with
*  I_A)  gamma = 2.78        N = 1.8250
*  I_B)  gamma = 2.88        N = 1.3700
*  I_C)  gamma = 2.69        N = 2.2800
*
*  Another parametrization is
*  resultat = N * beta_lorentz * E^{-gamma}
*  resultat = N * (impulsion/E) * pow(E,-gamma)
*  with
*  II_A)  gamma = 2.76        N = 1.5950     MEDIAN FLUX OF THE PAPER
*  II_B)  gamma = 2.85        N = 1.2300
*  II_C)  gamma = 2.67        N = 1.9600
*
*  Recently, FIORENZA, NICOLAO and SANDRO found less conservative
*  fits to the proton spectra.
*  resultat = A * beta_lorentz * E^{-gamma}
*  resultat = A * (impulsion/E) * pow(E,-gamma)
*  with
*  III_A)     A = 1.5300     gamma = 2.67
*  III_B)     A = 1.6600     gamma = 2.85
*  III_C)     A = 1.2300     gamma = 2.61     MINIMAL FLUX OF THE PAPER
*  III_D)     A = 1.9600     gamma = 2.89     MAXIMAL FLUX OF THE PAPER
*
**********************************************************************************************
*
* FEBRUARY 2000
* ANALYSIS BY FIORENZA DONATO
* Fitting both helium and proton data -- assuming the same spectral index --
* for the experiments Imax, Caprice, Mass, Leap and BESS, and AMS proton:
* 
*  N = 1.3215   gamma = 2.74   median  flux
*  N = 1.2750   gamma = 2.75   minimal flux 
*  N = 1.3680   gamma = 2.73   maximal flux
*
*  where we have chosen a parameterization N * rigidity^{- gamma}
*  Please NOTE that the correct variable is now the RIGIDITY.
*
* MARCH 2005
* FIORENZA HAS SLIGHTLY MODIFIED THE PROTON FLUX WITH
*
*  N = 1.3249   gamma = 2.72   modified median  flux
*
*  where we have chosen the NEW parameterization N * T^{- gamma}
*  Please NOTE that the correct variable is now the KINETIC ENERGY.
*
*/
double flux_proton_EXP(double E_proton)
{
	double T,impulsion,resultat;
	double A,N,gamma,p1,p2,p3,p4,Ep1,Ep2,Ep3,Ep4,R,beta;
	double phi_0,alpha_p;
  
	T = E_proton - MASSE_PROTON;
	if (T<=0.0)
	{
		return (0.0);
	}
	else
	{
		impulsion = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));

		//	OLD VERSION FRIDAY_001222
/*		N = 1.3215;
		gamma = 2.74;
		resultat = N  * pow(impulsion,-gamma);
*/
		//	NEW VERSION TUESDAY_050412
/*		N = 1.3249;
		gamma = 2.72;
		resultat = N  * pow(T,-gamma);
*/
/*		N = 1.3249;
		gamma = 2.72;
		resultat = N  * pow(T,-gamma);
*/
	

		#ifdef BESS_2008_proton_Shikaze
		//	New BESS data from Shikaze et al. VERSION AS OF 080903.
		
			A = 1.94;
			p1   = 0.7;
			p2   = 2.76;

			R    = sqrt(T*T + 2.*MASSE_PROTON*T);
			beta = R / (T+MASSE_PROTON);

			resultat = A * pow(beta,p1) * pow(R,-p2);


	
		#elif defined Fit_2008_proton_Maurin_Donato
		//	New parametrization from Fiorenza and David fits to H data. VERSION AS OF 081023.

    		if (T <= 20.0)
			{
		  		A    = 1.94;
		  		p1   = 0.7;
		 		p2   = 2.76;		  
			}
			else
			{
		  		A    = 2.4132;
		  		p1   = 0.0;
		  		p2   = 2.839;		  
			}

			R    = sqrt(T*T + 2.*MASSE_PROTON*T);
			beta = R / (T+MASSE_PROTON);

			resultat = A * pow(beta,p1) * pow(R,-p2);



		#elif defined CREAM_2010_proton_Lavalle
		//	New parameterization proposed by Julien Lavalle and based on the CREAM high energy CR proton data. The F1p fit is published in arXiv:1011.3063.

    		phi_0   =  3.09e-3;
			alpha_p =  2.8;
			Ep1     =  4.0;
			p1      =  1.05;
			Ep2     =  2.5e3;
			p2      =  0.34;
			Ep3     = 10.e3;
			p3      = -0.29;

    		resultat  = phi_0 * (1.0 - exp(-pow((T/Ep1),p1))) * pow((T/10.),(-alpha_p));
			resultat *= pow((1. + (T/Ep2)),p2) * pow((1. + (T/Ep3)),p3);
		
		
		
		#elif defined ATIC2_2010_proton_Lavalle
		//	New parameterization proposed by Julien Lavalle and based on the ATIC_2 high energy CR proton data. The F2p fit is published in arXiv:1011.3063.

    		phi_0   =  3.09e-3;
			alpha_p =  2.8;
			Ep1     =  4.0;
			p1      =  1.05;
			Ep2     =  1.5e3;
			p2      =  0.4;
			Ep3     = 10.e3;
			p3      = -0.35;

    		resultat  = phi_0 * (1.0 - exp(-pow((T/Ep1),p1))) * pow((T/10.),(-alpha_p));
			resultat *= pow((1. + (T/Ep2)),p2) * pow((1. + (T/Ep3)),p3);
			


		#elif defined PAMELA_2012_proton_Delahaye
		//	New parameterization proposed by Timur Delahaye and based on the PAMELA high energy CR proton data.

    		phi_0   =  3.53e-3;
			alpha_p =  2.5;
			Ep1     =  2.5;
			p1      =  0.9;
			Ep2     =  16.;
			p2      =  -0.5;
			Ep3     = 300.;
			p3      = 0.46;
			Ep4     = 5.e3;
			p4      = -0.21;

    		resultat  = phi_0 * (1.0 - exp(-pow((T/Ep1),p1))) * pow((T/10.),(-alpha_p));
			resultat *= pow((1. + (T/Ep2)),p2) * pow((1. + (T/Ep3)),p3) * pow((1. + (T/Ep4)),p4);
			
			
			
		#elif defined AMS02_2013_proton_Donato
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Fiorenza Donato in arXiv:1402.0321.

			A    = 2.2450;
			p1   = 2.32;
			p2   = 2.8232;

			R    = sqrt(T*T + 2.*MASSE_PROTON*T);
			beta = R / (T+MASSE_PROTON);

			resultat = A * pow(beta,p1) * pow(R,-p2);
			
			

		#elif defined AMS02_2013_proton_Kappl_Winkler
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Kappl and Winkler in arXiv:1408.0299.

			A          = 1.7407;
			gamma = 2.775;

			resultat = A * pow(T,-gamma);
			
		
		
		#elif defined AMS02_2013_proton_Vittino
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Vittino in a forthcoming paper.

			A          = 1.183042;
			gamma = 2.712;

			resultat = A * pow(T,-gamma);
		
			
			
		#else
			printf("ERROR : function 'flux_proton_EXP' \nYou must specify one proton flux parametrization in COMMON.h! \n");
			exit(0);
		#endif
			

/*	
	phi_0   =  3.53e-3;
	alpha_p =  2.5;
	Ep1     =  2.5;
	p1      =  0.9;
	Ep2     =  16.;
	p2      =  -0.5;
	Ep3     = 300.;
	p3      = 0.46;
	Ep4     = 5.e3;
	p4      = -0.21;

	resultat  = phi_0 * (1.0 - exp(-pow((T/Ep1),p1))) * pow((T/10.),(-alpha_p));
	resultat *= pow((1. + (T/Ep2)),p2) * pow((1. + (T/Ep3)),p3) * pow((1. + (T/Ep4)),p4); // [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] 

*/	  	
	return resultat;

	}
}

/********************************************************************************************/
/********************************************************************************************/
/*
* On calcule la fonction de production totale galactique de protons cosmiques.
* Elle s'exprime en [nombre de protons s^{-1} GeV^{-1}].
*
*/
double Q_proton_tot(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation)
{
  long i;
  double impulsion_proton,v_proton,K_proton;
  double production_E_proton,Si,Ai;
  
  impulsion_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
  v_proton         = CELERITY_LIGHT * impulsion_proton / E_proton;
  K_proton         =  K_space_diffusion(E_proton,MASSE_PROTON,1.0,pt_Propagation);
  
  production_E_proton = 0.0;

  for (i=1;i<=NDIM;i++)
  {
/*  Si est exprime en [kpc^{-1}].  */
    Si =
    sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_proton,2));
/*  Ai est exprime en [cm s^{-1}].  */
    Ai  = pt_Propagation->VENT_GALACTIQUE;
    Ai += 2.0*E_DISC*CM_PAR_KPC * (sigma_total_pH(E_proton) * v_proton * DENSITE_H_DISC);
    Ai += K_proton * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);
    
    production_E_proton += besselj0(alpha_i[i]*(R_EARTH/R_GAL)) * pt_Proton->q_i[i] / Ai;
/*  S'exprime en unites de [s^{+1} cm^{-3}]
*/
  }
  
  production_E_proton  = 4. * PI * flux_proton_EXP(E_proton) / v_proton /
  production_E_proton;
  return production_E_proton; /* [protons s^{-1} GeV^{-1}] */
}

/********************************************************************************************/
void calcul_method_A_BESSEL_COEF_i(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation)
{
  long i;
  double impulsion_proton,v_proton,K_proton;
  double production_E_proton,Si,Ai;
  
/*
* On remet a zero le tableau BESSEL_COEF_i[i].
*/
  for (i=0;i<=NDIM;i++)
  {
    pt_Proton->BESSEL_COEF_i[i] = 0.0;
  }

  impulsion_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
  v_proton         = CELERITY_LIGHT * impulsion_proton / E_proton;
  K_proton         =  K_space_diffusion(E_proton,MASSE_PROTON,1.0,pt_Propagation);
/*
* On calcule la fonction de production totale galactique.
*/
  production_E_proton = Q_proton_tot(E_proton, alpha_i, pt_Proton,pt_Propagation);

  for (i=1;i<=NDIM;i++)
  {
/*  Si est exprime en [kpc^{-1}].  */
    Si =
    sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_proton,2));
/*  Ai est exprime en [cm s^{-1}].  */
    Ai  = pt_Propagation->VENT_GALACTIQUE;
    Ai += 2.0*E_DISC*CM_PAR_KPC * (sigma_total_pH(E_proton) * v_proton * DENSITE_H_DISC);
    Ai += K_proton * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);
    
    pt_Proton->BESSEL_COEF_i[i] = production_E_proton * pt_Proton->q_i[i] / Ai;
/*  S'exprime en unites de [protons cm^{-3} GeV^{-1}].
*/
  }
  return;
}

/********************************************************************************************/
void calcul_method_B_BESSEL_Pi(double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation)
{
  long i;
  double impulsion_proton,v_proton,K_proton;
  double production_E_proton,Si,Ai;
  
/*
* On remet a zero le tableau BESSEL_COEF_i[i].
*/
  for (i=0;i<=NDIM;i++)
  {
    pt_Proton->BESSEL_COEF_i[i] = 0.0;
  }

  impulsion_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
  v_proton         = CELERITY_LIGHT * impulsion_proton / E_proton;
  K_proton         =  K_space_diffusion(E_proton,MASSE_PROTON,1.0,pt_Propagation);

  for (i=1;i<=NDIM;i++)
  {
/*  Si est exprime en [kpc^{-1}].  */
    Si =
    sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_proton,2));
/*  Ai est exprime en [cm s^{-1}].  */
    Ai  = pt_Propagation->VENT_GALACTIQUE;
    Ai += 2.0*E_DISC*CM_PAR_KPC * (sigma_total_pH(E_proton) * v_proton * DENSITE_H_DISC);
    Ai += K_proton * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);
    
    pt_Proton->BESSEL_COEF_i[i] = pt_Proton->q_i[i] / Ai;
/*  S'exprime pour l'instant en unites de [s^{+1} cm^{-3}].
*/
  }
  production_E_proton = flux_proton_EXP(E_proton) /
  GENERIC_FLUX(R_EARTH,0.,E_proton,MASSE_PROTON,1.,alpha_i,pt_Proton->BESSEL_COEF_i,pt_Propagation); /* [protons s^{-1} GeV^{-1}] */
  for (i=1;i<=NDIM;i++)
  {
    pt_Proton->BESSEL_COEF_i[i] *= production_E_proton;
/*  S'exprime maintenant en unites de [protons cm^{-3} GeV^{-1}].
*/
  }
  return;
}

/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Proton->BESSEL_COEF_Enuc_i
* en fonction de l'energie des protons E_protons et des coefficients de BESSEL i.
*
*/
void calculation_BESSEL_PROTON_Ep_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation)
{
  
 
  double T_proton,E_proton;
  long i_proton,i;
/*
* On remet a zero le tableau pt_Proton->BESSEL_COEF_Enuc_i.
*/
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Proton->BESSEL_COEF_Enuc_i[i_proton][i] = 0.0;
    }
  }
/*
* On le remplit maintenant.
*/
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
	  E_proton = E_PROTON_MIN *
		pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));
	//calcul_method_A_BESSEL_COEF_i(E_proton);
	  calcul_method_B_BESSEL_Pi(E_proton, alpha_i, pt_Proton,pt_Propagation);
    for (i=1;i<=NDIM;i++)
    {
      pt_Proton->BESSEL_COEF_Enuc_i[i_proton][i] = pt_Proton->BESSEL_COEF_i[i];
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
double flux_proton_TH(double r,double z,double E_proton, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation)
{
  double x,az;
  double impulsion_proton,v_proton,K_proton;
  long i;
  double Si,resultat;
  
  x = r/R_GAL;
  az = (z<0) ? -z : z;
  resultat = 0.0;
  
  if (az>=pt_Propagation->E_DIFFUS || x>=1.0)
  {
    return resultat;
  }
  else
  {
    impulsion_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
    v_proton         = CELERITY_LIGHT * impulsion_proton / E_proton;
    K_proton         =  K_space_diffusion(E_proton,MASSE_PROTON,1.0,pt_Propagation);

    for (i=1;i<=NDIM;i++)
    {
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_proton,2));
      resultat += pt_Proton->BESSEL_COEF_i[i] * besselj0(alpha_i[i]*x) *
      exp(pt_Propagation->VENT_GALACTIQUE*az*CM_PAR_KPC / (2.*K_proton)) *
      sinh((Si/2.)*(pt_Propagation->E_DIFFUS-az)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS);
    }
    resultat *= (1. / 4. / PI) * v_proton;
    return resultat; /* [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] */
  }
}

/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module calcule l'emissivite en antiprotons d'un atome d'hydrogene immerge dans
* le flux de protons flux_proton_EXP(E_proton) defini dans le fichier present PROTON.c
* Il convient de proceder a une integrale par la methode de SIMPSON.
* Cette integrale est menee de E_PROTON_MIN a E_PROTON_MAX en DIM_TAB_PROTON pas.
*
*/
double pbar_emissivity_per_H_solar(double E_pbar)
{
  double E_proton,dlog_E_proton,weight_SIMSPON;
  double resultat;
  long i_proton;

  dlog_E_proton = pow((E_PROTON_MAX/E_PROTON_MIN),(1./(double)DIM_TAB_PROTON)) - 1.;
  resultat = 0.0;
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    E_proton = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));
    
    if (i_proton==0 || i_proton==DIM_TAB_PROTON) {weight_SIMSPON = 1./3.;}
    else {weight_SIMSPON = (1. + (double)(i_proton % 2)) * 2. / 3.;}
    
    resultat += dSpbar_sur_dEpbar_SIMPSON(E_proton,E_pbar,100) * flux_proton_EXP(E_proton) *
    E_proton * weight_SIMSPON*dlog_E_proton;
  }
  resultat *= 4. * PI;
  return resultat; /* [antiprotons GeV^{-1} s^{-1} H^{-1}] */
}

/********************************************************************************************/
/********************************************************************************************/

void PROTON_IS_SPECTRUM_calculation(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1])
{
	long i_proton,i;
	double T_proton ,E_proton ,flux_proton_IS;

	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
		T_proton = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
		E_proton = T_proton + MASSE_PROTON;
	
  	  	calcul_method_B_BESSEL_Pi(E_proton, alpha_i, pt_Proton,pt_Propagation);
		
		flux_proton_IS = GENERIC_FLUX_04(R_EARTH,0.,E_proton,MASSE_PROTON,1.,alpha_i,pt_Proton->BESSEL_COEF_i, pt_Propagation);
		
		PROTON_IS_SPECTRUM[i_proton] = flux_proton_IS;																		// [#proton cm^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}
}
	
/********************************************************************************************/
/********************************************************************************************/

void PROTON_TOA_SPECTRUM_calculation(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double PROTON_TOA_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double T_PROTON_TOA[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Propagation* pt_Propagation)
{
	long i_proton;
	double T_proton_TOA,E_proton_TOA,flux_proton_TOA;
	double T_proton_IS ,E_proton_IS,flux_proton_IS;
	

	pt_Propagation->PHI_FISK = fisk_potential;
	

	for (i_proton=0;i_proton<=DIM_TAB_PBAR;i_proton++)
	{
		T_proton_IS = T_PROTON_SPECTRUM_MIN *
		pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
		E_proton_IS = T_proton_IS + MASSE_PROTON;

//		Nous modulons maintenant les spectres PBAR obtenus.

		FFA_IS_to_TOA(1.,1.,pt_Propagation->PHI_FISK,E_proton_IS,PROTON_IS_SPECTRUM[i_proton],&E_proton_TOA,&flux_proton_TOA);

		if (E_proton_TOA <= MASSE_PROTON)
		{
			T_PROTON_TOA[i_proton]        = 0.0;
			PROTON_TOA_SPECTRUM[i_proton] = 0.0;
			continue;
		}
		
		T_proton_TOA = E_proton_TOA - MASSE_PROTON;

//		Nous les stockons en memoire dans les tableaux RESULTS_T_PROTON_TOA[DIM_TAB_PBAR+1] et RESULTS_SPECTRUM_TOA_MIN_MED_MAX[DIM_TAB_PBAR+1];
	
		T_PROTON_TOA[i_proton] = T_proton_TOA;
		PROTON_TOA_SPECTRUM[i_proton] = flux_proton_TOA;																		// [#proton cm^{-3} sr^{-1} s^{-1} GeV^{-1}]
		
	}
}

/********************************************************************************************/
/********************************************************************************************/

void PROTON_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1])
{
	long i_proton;
	
	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
	    SPECTRUM[i_proton] = 0.0;
	}
	
}

/********************************************************************************************/
/********************************************************************************************/



