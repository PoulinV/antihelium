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
		
		#elif defined AMS02_2015_proton_manuela
		//	Fit performed by Manuela Vecchi on data from a forthcoming paper of AMS-02.
			
			
			//resultat = fit_proton_flux_AMS02_manuela_demodulated(E_proton);
			resultat = fit_proton_flux_AMS02_manuela_demodulated_2(E_proton);
			
		#elif defined AMS02_2015_proton_yoann
		//	Fit performed by Manuela Vecchi on data from a forthcoming paper of AMS-02.
	
			resultat = fit_proton_flux_AMS02_yoann_demodulated_2(E_proton);

		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_830MV
		  //AMS and CREAM parametrization with phi_fisk=830MV (2015)
			double C = 2.96486e+04, alpha = -0.693431, gamma = -2.90507, InvRb = 0.00290823, DeltaGamma = 0.250774, s = 0.17;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		    //AMS and CREAM parametrization with phi_fisk=783MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_783MV
			double C = 2.84802e+04, alpha = -0.597458, gamma = -2.89580, InvRb = 0.00256341, DeltaGamma = 0.250503, s = 0.17;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		    //AMS and CREAM parametrization with phi_fisk=755MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_755MV
			double C = 2.7863e+04, alpha = -0.552239, gamma = -2.89087, InvRb = 0.0024036, DeltaGamma = 0.249856, s = 0.169558;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		//AMS and CREAM parametrization with phi_fisk=724MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_724MV
			double C = 2.716e+04, alpha = -0.5115, gamma = -2.885, InvRb = 0.002357, DeltaGamma = 0.242, s = 0.1556;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		    //AMS and CREAM parametrization with phi_fisk=701MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_701MV
			double C = 2.669e+04, alpha = -0.484480, gamma = -2.88066, InvRb = 0.0023288, DeltaGamma = 0.236575, s = 0.146339;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		   //AMS and CREAM parametrization with phi_fisk=671MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_671MV
			double C = 2.6125e+04, alpha = -0.452629, gamma = -2.87576, InvRb = 0.0022972, DeltaGamma = 0.23017, s = 0.135759;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);

		   //AMS and CREAM parametrization with phi_fisk=647MV (2015)
		#elif defined AMS_CREAM_BCV_2015_H_phi_fisk_647MV
			double C = 2.57104e+04, alpha = -0.429501, gamma = -2.87218, InvRb = 0.00227701, DeltaGamma = 0.225418, s = 0.128178;
			double R = sqrt(T*T + 2.*MASSE_PROTON*T);
			double Beta = R/(T+MASSE_PROTON);
			double dR_On_dT = (T+MASSE_PROTON)/R;
			resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R*alpha)) * pow(R,gamma)*pow(1+pow((R*InvRb),DeltaGamma/s),s);
			
		#else
			printf("ERROR : function 'flux_proton_EXP' \nYou must specify one proton flux parametrization in COMMON.h! \n");
			exit(0);
		#endif
			

	  	
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
	

	pt_Propagation->PHI_FISK = phi_fisk;
	

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

double fit_proton_flux_AMS02_manuela(double R)
{
	double C, pow_alpha, pow_gamma, one_over_R_b, delta_gamma, s, Phi;
	double fisk_potential_AMS02;
	double r;
	
	fisk_potential_AMS02 = 0.62;
	C = 0.4556;
 	pow_alpha = -0.3835;
 	pow_gamma = -2.876;
 	one_over_R_b = 0.001927;
 	delta_gamma = 0.2563;
	
	
	// delta_gamma_plus
	//delta_gamma = 0.2563 + 0.1808;
	
	// delta_gamma_minus
	//delta_gamma = 0.2563 - 0.1808;
	
 	s = 0.1823;
 	Phi = 0.62;
	 
	r = C*(1.0 - exp(pow_alpha*R));
	r *= pow(R/(R+fisk_potential_AMS02), 2.0);
	r *= pow((R+fisk_potential_AMS02)/45, pow_gamma);
	r *= pow(1.0 + pow((R+fisk_potential_AMS02)*one_over_R_b, delta_gamma/s), s);
	r *= 1.0e-4;
	
	return r;
}

/********************************************************************************************/
/********************************************************************************************/

double fit_proton_flux_AMS02_yoann(double R)
{
	double C, pow_alpha, pow_gamma, one_over_R_b, delta_gamma, s, Phi;
	double fisk_potential_AMS02, beta;
	double r;
	
	fisk_potential_AMS02 = 0.62;
	C = 23566.0;
 	pow_alpha = -0.519;
 	pow_gamma = -2.849;
 	one_over_R_b = 1.0 / 355.0;
 	delta_gamma = 0.146;
	beta = 1.21;
	
	
 	s = 0.0325;
 	Phi = 0.62;
	 
	r = C*(1.0 - beta*exp(pow_alpha*R));
	r *= pow(R/(R+fisk_potential_AMS02), 2.0);
	r *= pow((R+fisk_potential_AMS02), pow_gamma);
	r *= pow(1.0 + pow((R+fisk_potential_AMS02)*one_over_R_b, delta_gamma/s), s);
	r *= 1.0e-4;
	
	return r;
}


/********************************************************************************************/
/********************************************************************************************/

double fit_proton_flux_AMS02_manuela_demodulated(double EnIS)
{
    double pnIS,pnTOA;
    double pnC,EnC,En_min,En_trans;
	double EnTOA;
	double r;
	double Z, A, PHI;
	
	
	Z = 1.0;
	A = 1.0;
  	PHI = 0.62;
  
    pnC = Z * RIGIDITY_MS_C / A; /* Dans nos unites, la charge de l'electron vaut 1. */
    EnC = sqrt(pow(MASSE_PROTON,2) + pow(pnC,2));
  
    En_trans = EnC + Z*PHI/A;
    En_min   = En_trans + pnC*log(MASSE_PROTON/(EnC+pnC));
	
	//printf("En_min = %.5e \t En_trans = %.5e \t  EnC = %.5e \n", En_min, En_trans, EnC);
	
	
   	if (EnIS <= En_trans)
    {
      EnTOA = MASSE_PROTON * cosh((EnIS - En_min)/pnC);
      pnIS  = sqrt(pow(EnIS,2) - pow(MASSE_PROTON,2));
      pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
	  
	  r = pow((pnIS/pnTOA),2.) * fit_proton_flux_AMS02_manuela(pnTOA);
	 // printf("EnTOA = %.5e \t pnIS = %.5e \t pnTOA = %.5e \n", EnTOA, pnIS, pnTOA);

    }
    else
    {
      EnTOA = EnIS - Z*PHI/A;
      pnIS  = sqrt(pow(EnIS,2) - pow(MASSE_PROTON,2));
      pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));

	  r = pow((pnIS/pnTOA),2.) * fit_proton_flux_AMS02_manuela(pnTOA);
	  //printf("EnTOA = %.5e \t pnIS = %.5e \t pnTOA = %.5e \n", EnTOA, pnIS, pnTOA);
	  
    }
	
	r *= EnIS / sqrt(pow(EnIS,2) - pow(MASSE_PROTON,2));
	

	//printf("r = %.5e \n", r);
	
	return r;
}

/********************************************************************************************/
/********************************************************************************************/

double fit_proton_flux_AMS02_manuela_demodulated_2(double EnIS)
{
    double pnIS,pnTOA;
    double pnC,EnC,En_min,En_trans;
	double EnTOA;
	double r;
	double Z, A, PHI;
	//double TnIS;
	
	double flux_IS, flux_TOA;
	
	flux_IS = 0.0;
	
	
	FFA_IS_to_TOA_modified(1.0, 1.0, 0.62, EnIS, flux_IS, &EnTOA, &flux_TOA);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
	
	flux_TOA = fit_proton_flux_AMS02_manuela(pnTOA);
	flux_TOA *= EnTOA / pnTOA;
	//printf("EnTOA = %.5e \t flux_TOA = %.5e  \n", EnTOA, flux_TOA);
	
	
	FFA_TOA_to_IS(1.0, 1.0, 0.62, EnTOA, flux_TOA, &EnIS, &flux_IS);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
	r = flux_IS;
	
	
	
	//printf("r = %.5e \n", r);
	
	return r;
}


/********************************************************************************************/
/********************************************************************************************/

double fit_proton_flux_AMS02_yoann_demodulated_2(double EnIS)
{
    double pnIS,pnTOA;
    double pnC,EnC,En_min,En_trans;
	double EnTOA;
	double r;
	double Z, A, PHI;
	//double TnIS;
	
	double flux_IS, flux_TOA;
	
	flux_IS = 0.0;
	
	
	FFA_IS_to_TOA_modified(1.0, 1.0, 0.62, EnIS, flux_IS, &EnTOA, &flux_TOA);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
	
	flux_TOA = fit_proton_flux_AMS02_yoann(pnTOA);
	flux_TOA *= EnTOA / pnTOA;
	
	//printf("EnTOA = %.5e \t flux_TOA = %.5e  \n", EnTOA, flux_TOA);
	
	
	FFA_TOA_to_IS(1.0, 1.0, 0.62, EnTOA, flux_TOA, &EnIS, &flux_IS);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
	r = flux_IS;
	
	
	
	//printf("r = %.5e \n", r);
	
	return r;
}


/********************************************************************************************/
/********************************************************************************************/

void FFA_IS_to_TOA_modified(double A,double Z,double PHI,
double EnIS,double flux_IS,double *EnTOA,double *flux_TOA)
{
  double pnIS,pnTOA;
  double pnC,EnC,En_min,En_trans;
  
  pnC = Z * RIGIDITY_MS_C / A; /* Dans nos unites, la charge de l'electron vaut 1. */
  EnC = sqrt(pow(MASSE_PROTON,2) + pow(pnC,2));
  
  En_trans = EnC + Z*PHI/A;
  En_min   = En_trans + pnC*log(MASSE_PROTON/(EnC+pnC));

/*  
  if (EnIS <= En_min)
  {
    printf(
    " TON NUCLEON EST TROP MOU : IL N'ARRIVERA JAMAIS SUR TERRE !\n"
    " SON ENERGIE INTER_STELLAIRE EST INFERIEURE A LA VALEUR CRITIQUE = %.5e [GEV]\n"
    " PRENDS UNE ENERGIE SUPERIEURE \n",
    En_min);
    *EnTOA = MASSE_PROTON;
    *flux_TOA = 0.0;
    return;
  }
 */
 /* else*/ if (EnIS <= En_trans)
  {
    *EnTOA = MASSE_PROTON * cosh((EnIS - En_min)/pnC);
    pnIS  = sqrt(pow(EnIS,2) - pow(MASSE_PROTON,2));
    pnTOA = sqrt(pow(*EnTOA,2) - pow(MASSE_PROTON,2));
    *flux_TOA = pow((pnTOA/pnIS),2) * flux_IS;
    return;
  }
  else
  {
    *EnTOA = EnIS - Z*PHI/A;
    pnIS  = sqrt(pow(EnIS,2) - pow(MASSE_PROTON,2));
    pnTOA = sqrt(pow(*EnTOA,2) - pow(MASSE_PROTON,2));
    *flux_TOA = pow((pnTOA/pnIS),2) * flux_IS;
    return;
  }
}


/********************************************************************************************/
/********************************************************************************************/









