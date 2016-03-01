#include "HELIUM.h"

/********************************************************************************************/
/********************************************************************************************/
/*
*  Il s'agit du flux interstellaire des noyaux d'helium. Ce flux est un flux
* differentiel par rapport a l'ENERGIE PAR NUCLEON. Il est donc exprime en unite
* de [cm^{-2} s^{-1} sr^{-1} (GeV/nucleon)^{-1}].
*
* FEBRUARY 2000
* ANALYSIS BY FIORENZA DONATO
* Fitting both helium and proton data -- assuming the same spectral index --
* for the experiments Imax, Caprice, Mass, Leap and BESS, and AMS proton:
* 
*  N = 0.0645   gamma = 2.70   median  flux
*
*  where we have chosen a parameterization N * T_nucleon^{- gamma}
*  Please NOTE that the correct variable is now the KINETIC ENERGY PER NUCLEON.
*
* MARCH 2005
* FIORENZA HAS SLIGHTLY MODIFIED THE HELIUM FLUX WITH
*
*  N = 0.0721   gamma = 2.74   modified median  flux
*
*/
double flux_helium_EXP(double E_nucleon)
{
 	double T_nucleon,P_nucleon,rigidity,resultat;
	double A,N,gamma;
	double beta,R,Rp1,Rp2,Rp3,Rp4,p1,p2,p3,p4,Ttot,Mtot;
  
	T_nucleon = E_nucleon - MASSE_PROTON;
  	if (T_nucleon<=0.0)
  	{
		return (0.0);
  	}
  	else
  	{
    	P_nucleon = sqrt(pow(E_nucleon,2) - pow(MASSE_PROTON,2));

		//   NEW VERSION TUESDAY_050412
/*
    	rigidity  = (A_HE / Z_HE) * P_nucleon;
    	N = 0.0721;
    	gamma = 2.74;
    	resultat = N  * pow(T_nucleon,-gamma);

    	rigidity  = (A_HE / Z_HE) * P_nucleon;
    	N = 0.0721;
    	gamma = 2.74;
    	resultat = N  * pow(T_nucleon,-gamma);
*/
		
		
		#ifdef BESS_2008_helium_Shikaze
		//	New BESS data from Shikaze et al. VERSION AS OF 080903.

   	 		A  = 0.71;
			p1 = 0.5;
			p2 = 2.78;

			Ttot = 4. * T_nucleon;
			Mtot = 4. * MASSE_PROTON;
			R    = sqrt(Ttot * (Ttot + 2.*Mtot)) / 2.;
			beta = 2. * R / (Ttot + Mtot);

			resultat = A * pow(beta,p1) * pow(R,-p2);
			
			
		
		#elif defined Fit_2008_helium_Maurin_Donato
		//	New parametrization from Fiorenza and David fits to HE data. VERSION AS OF 081023.

    		if (T_nucleon <= 20.0)
			{
		  	  	A  = 0.71;
		  		p1 = 0.5;
		  	  	p2 = 2.78;		  
			}
			else
			{
		  	  	A    = 0.8866;
		  	  	p1   = 0.0;
		  	  	p2   = 2.85;		  
			}

			Ttot = 4. * T_nucleon;
			Mtot = 4. * MASSE_PROTON;
			R    = sqrt(Ttot * (Ttot + 2.*Mtot)) / 2.;
			beta = 2. * R / (Ttot + Mtot);

			resultat = A * pow(beta,p1) * pow(R,-p2);



		#elif defined CREAM_ATIC2_2010_helium_Lavalle
		//	New parameterization proposed by Julien Lavalle and based on the CREAM and ATIC_2 high energy CR helium data. The F1He fit is published in arXiv:1011.3063.

    		A   =  0.71;
			p1  =  0.5;
			p2  =  2.78;
			p3  =  0.5;
			Rp3 =  1.e3;
			p4  = -0.5;
			Rp4 = 10.e3;

			Ttot = 4. * T_nucleon;
			Mtot = 4. * MASSE_PROTON;
			R    = sqrt(Ttot * (Ttot + 2.*Mtot)) / 2.;
			beta = 2. * R / (Ttot + Mtot);

			resultat  = A * pow(beta,p1) * pow(R,-p2);
			resultat *= pow((1. + (R/Rp3)),p3) * pow((1. + (R/Rp4)),p4);

			
			
		#elif defined PAMELA_2012_helium_Delahaye
		//	New parameterization proposed by Timur Delahaye and based on the PAMELA high energy CR helium data.

    		A       = 1.5e-5;
			Rp1     = 50.;
			p1      = 2.7;
			Rp2     = 250.;
			p2      = -1.3;
			Rp3     = 1.e3;
			p3      = 5.4;
			Rp4     = 2.e3;
			p4      = -4.15;

			Ttot = 4. * T_nucleon;
			Mtot = 4. * MASSE_PROTON;
			R    = sqrt(Ttot * (Ttot + 2.*Mtot)) / 2.;
			beta = 2. * R / (Ttot + Mtot);

    		resultat  = A * pow((R/Rp1),-p1);
			resultat *= pow((1. + (R/Rp2)),p2) * pow((1. + (R/Rp3)),p3) * pow((1. + (R/Rp4)),p4);



		#elif defined AMS02_2013_helium_Donato
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Fiorenza Donato in arXiv:1402.0321.

			A  = 0.5220;
			p1 = 1.34;
			p2 = 2.6905;

			Ttot = 4. * T_nucleon;
			Mtot = 4. * MASSE_PROTON;
			R    = sqrt(Ttot * (Ttot + 2.*Mtot)) / 2.;
			beta = 2. * R / (Ttot + Mtot);

			resultat = A * pow(beta,p1) * pow(R,-p2);



		#elif defined AMS02_2013_helium_Kappl_Winkler
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Kappl and Winkler in arXiv:1408.0299.

			A          = 0.05972;
			gamma = 2.630;

			resultat = A * pow(T_nucleon,-gamma);
			
			
			
		#elif defined AMS02_2013_helium_Vittino
		//	New AMS02 data presented at ICRC 2013 in Rio de Janeiro Parameterized by Vittino in a forthcoming paper.

			A          = 0.06576578;
			gamma = 2.667;

			resultat = A * pow(T_nucleon,-gamma);
		
		#elif defined AMS02_2015_helium_yoann
		//	Fit performed by Yoann on data from AMS days
	
			resultat = fit_helium_flux_AMS02_yoann_demodulated_2(E_nucleon);
	
		
		#else
			printf("ERROR : function 'flux_helium_EXP' \nYou must specify one helium flux parametrization in COMMON.h! \n");
			exit(0);
		#endif



		return resultat;
	}
}

/********************************************************************************************/
/********************************************************************************************/
void calcul_method_B_BESSEL_HEi(double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation)
{
  long i;
  double impulsion_nucleon,v_helium,K_helium;
  double production_E_helium,Si,Ai;
  
/*
* On remet a zero le tableau BESSEL_COEF_i[i].
*/
  for (i=0;i<=NDIM;i++)
  {
    pt_Helium->BESSEL_COEF_i[i] = 0.0;
  }

  impulsion_nucleon = sqrt(pow(E_nucleon,2) - pow(MASSE_PROTON,2));
  v_helium          = CELERITY_LIGHT * impulsion_nucleon / E_nucleon;
/*
  K_helium          = K_space_diffusion(E_nucleon,MASSE_PROTON,(Z_HE / A_HE));
*/
  K_helium          = K_space_diffusion((A_HE*E_nucleon),(A_HE*MASSE_PROTON),Z_HE,pt_Propagation);

  for (i=1;i<=NDIM;i++)
  {
/*  Si est exprime en [kpc^{-1}].  */
    Si =
    sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_helium,2));
/*  Ai est exprime en [cm s^{-1}].  */
    Ai  = pt_Propagation->VENT_GALACTIQUE;
    Ai += 2.0*E_DISC*CM_PAR_KPC *
    (pow(A_HE,2.2/3.) * sigma_total_pH(E_nucleon) * v_helium * DENSITE_H_DISC);
    Ai += K_helium * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);
    
    pt_Helium->BESSEL_COEF_i[i] = pt_Helium->q_i[i] / Ai;
/*  S'exprime pour l'instant en unites de [s^{+1} cm^{-3}].
*/
  }
/*
  production_E_helium = flux_helium_EXP(E_nucleon) /
  GENERIC_FLUX(R_EARTH,0.,E_nucleon,MASSE_PROTON,(Z_HE / A_HE),pt_Helium->BESSEL_COEF_i);
*/
  production_E_helium = flux_helium_EXP(E_nucleon) /
  GENERIC_FLUX(R_EARTH,0.,(A_HE*E_nucleon),(A_HE*MASSE_PROTON),Z_HE,alpha_i,pt_Helium->BESSEL_COEF_i,pt_Propagation);
/*
* S'exprime en [helions s^{-1} (GeV/nucleon)^{-1}]
*/
  for (i=1;i<=NDIM;i++)
  {
    pt_Helium->BESSEL_COEF_i[i] *= production_E_helium;
/*  S'exprime maintenant en unites de [helions cm^{-3} (GeV/nucleon)^{-1}].
*/
  }
  return;
}

/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Helium->BESSEL_COEF_Enuc_i
* en fonction de l'energie par nucleon E_nucleon des noyaux d'helium  et des
* coefficients de BESSEL i.
*
*/
void calculation_BESSEL_HELIUM_Ep_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation)
{

  
  double T_nucleon,E_nucleon;
  long i_nucleon,i;
/*
* On remet a zero le tableau pt_Helium->BESSEL_COEF_Enuc_i.
*/
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Helium->BESSEL_COEF_Enuc_i[i_nucleon][i] = 0.0;
    }
  }
/*
* On le remplit maintenant.
*/
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
	  E_nucleon = E_PROTON_MIN *
		pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_nucleon/(double)DIM_TAB_PROTON));
    calcul_method_B_BESSEL_HEi(E_nucleon, alpha_i, pt_Helium,pt_Propagation);
    for (i=1;i<=NDIM;i++)
    {
     pt_Helium->BESSEL_COEF_Enuc_i[i_nucleon][i] = pt_Helium->BESSEL_COEF_i[i];
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
double flux_helium_TH(double r,double z,double E_nucleon, double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Helium, struct Structure_Propagation* pt_Propagation)
{
  double x,az;
  double impulsion_nucleon,v_helium,K_helium;
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
    impulsion_nucleon = sqrt(pow(E_nucleon,2) - pow(MASSE_PROTON,2));
    v_helium          = CELERITY_LIGHT * impulsion_nucleon / E_nucleon;
/*
    K_helium          = K_space_diffusion(E_nucleon,MASSE_PROTON,(Z_HE / A_HE));
*/
    K_helium          = K_space_diffusion((A_HE*E_nucleon),(A_HE*MASSE_PROTON),Z_HE,pt_Propagation);

    for (i=1;i<=NDIM;i++)
    {
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_helium,2));
      resultat += pt_Helium->BESSEL_COEF_i[i] * besselj0(alpha_i[i]*x) *
      exp(pt_Propagation->VENT_GALACTIQUE*az*CM_PAR_KPC / (2.*K_helium)) *
      sinh((Si/2.)*(pt_Propagation->E_DIFFUS-az)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS);
    }
    resultat *= (1. / 4. / PI) * v_helium;
    return resultat; /* [cm^{-2} s^{-1} sr^{-1} (GeV/nucleon)^{-1}] */
  }
}

/********************************************************************************************/
/********************************************************************************************/

double fit_helium_flux_AMS02_yoann(double R)
{
	double C, pow_alpha, pow_gamma, one_over_R_b, delta_gamma, s, Phi;
	double fisk_potential_AMS02, beta;
	double r;
	
	fisk_potential_AMS02 = 0.62;
	C = 4075.0;
 	pow_alpha = -0.163;
	beta = 0.41;
 	pow_gamma = -2.795;
 	one_over_R_b = 1.0/284.0;
 	delta_gamma = 0.162;
	
	// delta_gamma_plus
	//delta_gamma = 0.2563 + 0.1808;
	
	// delta_gamma_minus
	//delta_gamma = 0.2563 - 0.1808;
	
 	s = 0.078;
	 
	r = C*(1.0 - beta*exp(pow_alpha*R));
	r *= pow(R/(R+fisk_potential_AMS02), 2.0);
	r *= pow((R+fisk_potential_AMS02), pow_gamma);
	r *= pow(1.0 + pow((R+fisk_potential_AMS02)*one_over_R_b, delta_gamma/s), s);
	r *= 1.0e-4;
	
	return r;
}

/********************************************************************************************/
/********************************************************************************************/

double fit_helium_flux_AMS02_yoann_demodulated_2(double EnIS)
{
    double pnIS,pnTOA;
    double pnC,EnC,En_min,En_trans;
	double EnTOA;
	double r;
	double Z, A, PHI;
	//double TnIS;
	
	double flux_IS, flux_TOA;
	
	flux_IS = 0.0;
	
	
	FFA_IS_to_TOA_modified(4.0, 2.0, 0.62, EnIS, flux_IS, &EnTOA, &flux_TOA);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
	
	flux_TOA = fit_helium_flux_AMS02_yoann(2.0*pnTOA);
	flux_TOA *= 2.0*EnTOA / pnTOA;
	
	
	FFA_TOA_to_IS(4.0, 2.0, 0.62, EnTOA, flux_TOA, &EnIS, &flux_IS);
	
	//printf("EnIS = %.5e \t EnTOA = %.5e \n", EnIS, EnTOA);
	
	r = flux_IS;
		
	
	//printf("r = %.5e \n", r);
	
	return r;
}

/********************************************************************************************/
/********************************************************************************************/


