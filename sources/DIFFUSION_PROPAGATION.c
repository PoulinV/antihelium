#include "DIFFUSION_PROPAGATION.h"

/********************************************************************************************/
/********************************************************************************************/
/*
*
* K = DIFFUSION_0_GV * beta * pow(rigidity,PUISSANCE_COEFF_DIFF);
*
*/
double K_space_diffusion(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation)
{
  double momentum,rigidity,beta,K;

  momentum = sqrt(pow(energy,2.) - pow(mass,2.)); /* [GeV]        */
  rigidity = momentum / Z_em;                     /* [GV]         */
  beta     = momentum / energy;                   /* NO DIMENSION */

/*
* Coefficient de diffusion spatiale a la Strong & al, ApJ 509, 212, 1998.
* Il est exprime en [cm^{2} s^{-1}].
*/
  K        = pt_Propagation->DIFFUSION_0_GV * beta * pow(rigidity,pt_Propagation->PUISSANCE_COEFF_DIFF);
  return K;
}

/********************************************************************************************/
/*
* This function computes the diffusion coefficient D_EE in energy space that
* is responsible for the diffusive reacceleration. The cosmic ray spectrum gets
* smeared by this process and the average particle energy tends to increase.
*
* Note that both the TOTAL energy and the TOTAL mass of the cosmic ray species
* are considered here. Whenever one deals with energies per nucleon in the rest
* of the code, conversion into TOTAL energy must be properly performed as this
* function is called.
*
* Energies and masses are expressed in [GeV].
*
*/
double D_energy_diffusion(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation)
{
  double momentum,rigidity,beta,D_EE;

  momentum = sqrt(pow(energy,2) - pow(mass,2)); /* [GeV]        */
  beta     = momentum / energy;                 /* NO DIMENSION */

  D_EE     = (2./9.) * pow(pt_Propagation->V_ALFEN,2) * pow(beta,4) * pow(energy,2) /
  K_space_diffusion(energy,mass,Z_em,pt_Propagation); /* [GeV^{2} s^{-1}] */

  return D_EE;
}

/********************************************************************************************/
/*
* This module computes the energy LOSSES of a given cosmic ray species
* that result from IONIZATION as well as COULOMB FRICTION on the gas of
* the galactic ridge. As the function describes the variation in time of
* the TOTAL energy, note the presence of a MINUS sign.
*
* Note that both the TOTAL energy and the TOTAL mass of the cosmic ray
* species are considered here. Whenever one deals with energies per nucleon
* in the rest of the code, conversion into TOTAL energy must be properly
* performed whenever this function is called.
*
* Energies and masses are expressed in [GeV].
*
*/
double b_energy_losses(double energy,double mass,double Z_em, struct Structure_Propagation* pt_Propagation)
{
  double momentum,beta,gamma_lorentz;
  double Q_max,B1,B2,ionisation;
  double x_m,ln_lambda,coulomb,adiabatic;
  double dE_dt;
/*
* La perte d'energie dE_dt est exprimee en [GeV sec^{-1}].
*/
  momentum      = sqrt(pow(energy,2) - pow(mass,2)); /* [GeV]        */
  beta          = momentum / energy;                 /* NO DIMENSION */
  gamma_lorentz = energy / mass;                     /* NO DIMENSION */

/*
* IONISATION CONTRIBUTION
* ***********************
* L'energie maximale transferee est notee Q_max et est exprimee en [GeV].
*/
  Q_max = 2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) /
  (1.0  +  (2. * gamma_lorentz * MASSE_ELECTRON / mass));
/*
* Les coefficients B1 et B2 sont sans dimension.
*/
  B1  = log(2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) * Q_max / pow(V_ION_H ,2.));
  B1 -= 2. * pow(beta,2.);
  B2  = log(2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) * Q_max / pow(V_ION_HE,2.));
  B2 -= 2. * pow(beta,2.);
/*
* La perte d'energie par ionisation est exprimee en [GeV s^{-1}].
*/
  ionisation =
  (-1.) * 2.*PI*pow(RADIUS_ELECTRON,2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em,2) *
  ((DENSITE_H_DISC * B1) + (DENSITE_HE_DISC * B2)) / beta; /* [GeV s^{-1}] */

/*
* COULOMB CONTRIBUTION
* ********************
*/
  x_m = pow((3.*sqrt(PI)/4.),(1./3.)) *
  sqrt(2. * K_BOLTZMANN * T_ELECTRONIC / MASSE_ELECTRON);

  ln_lambda  = mass / (mass + 2.*gamma_lorentz*MASSE_ELECTRON);
  ln_lambda *= pow((MASSE_ELECTRON/H_BAR/CELERITY_LIGHT),2.) / RADIUS_ELECTRON; /* [cm^{-3}] */
  ln_lambda *= pow(gamma_lorentz,2.) * pow(beta,4.) / PI / DENSITE_FREE_ELECTRON;
  ln_lambda  = 0.5 * log(ln_lambda);

  coulomb =
  (-1.) * 4.*PI*pow(RADIUS_ELECTRON,2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em,2) *
  DENSITE_FREE_ELECTRON *
  ln_lambda * pow(beta,2.) / ( pow(x_m,3.) + pow(beta,3.) ); /* [GeV s^{-1}] */

/*
* ADIABATIC CONTRIBUTION
* **********************
*
* ATTENTION ! Pour pouvoir ecrire ce terme sous une forme exactement similaire
* aux pertes precedentes, il faut factoriser le terme 2h, d'ou un terme
* \beq
* b_adiab \; = \; - \, E \, \frac{V_{c}}{3h} \;\; .
* \eeq
*/
  adiabatic =
	(-1.) * (momentum*momentum/energy) *
	(pt_Propagation->VENT_GALACTIQUE/(3.*E_DISC*CM_PAR_KPC)); /* [GeV s^{-1}] */

/*
* WE SUM UP BOTH CONTRIBUTIONS
* ****************************
* dE_dt = ionisation + coulomb + adiabatic;
* dE_dt = ionisation + coulomb;
*
*/
  dE_dt = ionisation + coulomb + adiabatic;
  return dE_dt;
}

/********************************************************************************************/
/*
* This module computes the energy LOSSES of a given cosmic ray species
* that result from IONIZATION, COULOMB FRICTION on the gas of
* the galactic ridge and ADIABATIC CONVECTION. As these processes result
* into a LOSS of energy, there should be a MINUS sign.
*
* This module returns also the DRIFT term that comes into play when
* we start with the phase space distribution function f and not from
* the energy distribution function $\psi$.
*
* Energies and masses are expressed in [GeV].
*
*/
double b_energy_losses_for_wei_calculations(double energy,double mass,double Z_em,double *ionisation,double *coulomb,double *adiabatic,double *drift, struct Structure_Propagation* pt_Propagation)
{
  double momentum,beta,gamma_lorentz;
  double Q_max,B1,B2;
  double x_m,ln_lambda;
  double dE_dt;
/*
* La perte d'energie dE_dt est exprimee en [GeV sec^{-1}].
*/
  momentum      = sqrt(pow(energy,2) - pow(mass,2)); /* [GeV]        */
  beta          = momentum / energy;                 /* NO DIMENSION */
  gamma_lorentz = energy / mass;                     /* NO DIMENSION */

/*
* IONISATION CONTRIBUTION
* ***********************
* L'energie maximale transferee est notee Q_max et est exprimee en [GeV].
*/
  Q_max = 2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) /
  (1.0  +  (2. * gamma_lorentz * MASSE_ELECTRON / mass));
/*
* Les coefficients B1 et B2 sont sans dimension.
*/
  B1  = log(2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) * Q_max / pow(V_ION_H ,2.));
  B1 -= 2. * pow(beta,2.);
  B2  = log(2. * MASSE_ELECTRON * pow((beta*gamma_lorentz),2.) * Q_max / pow(V_ION_HE,2.));
  B2 -= 2. * pow(beta,2.);
/*
* La perte d'energie par ionisation est exprimee en [GeV s^{-1}].
*/
  *ionisation =
  (-1.) * 2.*PI*pow(RADIUS_ELECTRON,2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em,2) *
  ((DENSITE_H_DISC * B1) + (DENSITE_HE_DISC * B2)) / beta; /* [GeV s^{-1}] */

/*
* COULOMB CONTRIBUTION
* ********************
*/
  x_m = pow((3.*sqrt(PI)/4.),(1./3.)) *
  sqrt(2. * K_BOLTZMANN * T_ELECTRONIC / MASSE_ELECTRON);

  ln_lambda  = mass / (mass + 2.*gamma_lorentz*MASSE_ELECTRON);
  ln_lambda *= pow((MASSE_ELECTRON/H_BAR/CELERITY_LIGHT),2.) / RADIUS_ELECTRON; /* [cm^{-3}] */
  ln_lambda *= pow(gamma_lorentz,2.) * pow(beta,4.) / PI / DENSITE_FREE_ELECTRON;
  ln_lambda  = 0.5 * log(ln_lambda);

  *coulomb =
  (-1.) * 4.*PI*pow(RADIUS_ELECTRON,2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em,2) *
  DENSITE_FREE_ELECTRON *
  ln_lambda * pow(beta,2.) / ( pow(x_m,3.) + pow(beta,3.) ); /* [GeV s^{-1}] */

/*
* ADIABATIC CONTRIBUTION
* **********************
*
* ATTENTION ! Pour pouvoir ecrire ce terme sous une forme exactement similaire
* aux pertes precedentes, il faut factoriser le terme 2h, d'ou un terme
* \beq
* b_adiab \; = \; - \, E \, \frac{V_{c}}{3h} \;\; .
* \eeq
*/
  *adiabatic =
	(-1.) * (momentum*momentum/energy) *
	(pt_Propagation->VENT_GALACTIQUE/(3.*E_DISC*CM_PAR_KPC)); /* [GeV s^{-1}] */

/*
* DRIFT TERM
* **********
*
* This term comes from the formalism based on the phase space distribution function.
* It may be expressed as
*
* \beq
* b_drift \; = \; \frac{1 \, + \, \beta^{2}}{\beta^{2}} \times \frac{K_EE}{E} \;\; .
* \eeq
*/
  *drift =
	(+1.) * (1. + beta*beta) * D_energy_diffusion(energy,mass,Z_em,pt_Propagation) / energy / (beta*beta); /* [GeV s^{-1}] */

/*
* WE SUM UP ALL THE LOSS CONTRIBUTIONS
* ************************************
* dE_dt = ionisation + coulomb + adiabatic;
*
*/
  dE_dt = *ionisation + *coulomb + *adiabatic;
  return dE_dt;
}

/********************************************************************************************/
/********************************************************************************************/
double GENERIC_FLUX(double r,double z,
double energy,double mass,double Z_em, double alpha_i[NDIM+1],double BESSEL_COEFFICIENTi[NDIM+1], struct Structure_Propagation* pt_Propagation)
{
  double x,az;
  double momentum,velocity,K;
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
    momentum = sqrt(pow(energy,2) - pow(mass,2));
    velocity = CELERITY_LIGHT * momentum / energy;
    K        = K_space_diffusion(energy,mass,Z_em,pt_Propagation);
    
    for (i=1;i<=NDIM;i++)
    {
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K,2));
      resultat += BESSEL_COEFFICIENTi[i] * besselj0(alpha_i[i]*x) *
      exp(pt_Propagation->VENT_GALACTIQUE*az*CM_PAR_KPC / (2.*K)) *
      sinh((Si/2.)*(pt_Propagation->E_DIFFUS-az)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS);
    }
    resultat *= (1. / 4. / PI) * velocity;
    return resultat; /* [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] */
  }
}

/********************************************************************************************/
/********************************************************************************************/

double GENERIC_FLUX_04(double r,double z,
double energy,double mass,double Z_em,double alpha_i[NDIM+1],double BESSEL_COEFFICIENTi[NDIM+1], struct Structure_Propagation* pt_Propagation)
{
  double x,az;
  double momentum,velocity,K;
  long i;
  double Si,resultat;
  double coefficient_torsten = 1.0;
  
  x = r/R_GAL;
  az = (z<0) ? -z : z;
  resultat = 0.0;
  
  if (az>=pt_Propagation->E_DIFFUS || x>=1.0)
  {
    return resultat;
  }
  else
  {
    momentum = sqrt(pow(energy,2) - pow(mass,2));
    velocity = CELERITY_LIGHT * momentum / energy;
    K        = K_space_diffusion(energy,mass,Z_em,pt_Propagation);
	
    
    for (i=1;i<=NDIM;i++)
    {
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K,2));
	

      if (i==(NDIM-15)){coefficient_torsten = 624./625.;}
			if (i==(NDIM-14)){coefficient_torsten = 620./625.;}
			if (i==(NDIM-13)){coefficient_torsten = 610./625.;}
			if (i==(NDIM-12)){coefficient_torsten = 590./625.;}
			if (i==(NDIM-11)){coefficient_torsten = 555./625.;}
			if (i==(NDIM-10)){coefficient_torsten = 503./625.;}
			if (i==(NDIM-9)) {coefficient_torsten = 435./625.;}
			if (i==(NDIM-8)) {coefficient_torsten = 355./625.;}
			if (i==(NDIM-7)) {coefficient_torsten = 270./625.;}
			if (i==(NDIM-6)) {coefficient_torsten = 190./625.;}
			if (i==(NDIM-5)) {coefficient_torsten = 122./625.;}
			if (i==(NDIM-4)) {coefficient_torsten =  70./625.;}
			if (i==(NDIM-3)) {coefficient_torsten =  35./625.;}
			if (i==(NDIM-2)) {coefficient_torsten =  15./625.;}
			if (i==(NDIM-1)) {coefficient_torsten =   5./625.;}
			if (i==(NDIM))   {coefficient_torsten =   1./625.;}

      resultat += coefficient_torsten * BESSEL_COEFFICIENTi[i] * besselj0(alpha_i[i]*x) *
      exp(pt_Propagation->VENT_GALACTIQUE*az*CM_PAR_KPC / (2.*K)) *
	  sinh((Si/2.)*(pt_Propagation->E_DIFFUS-az)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS);
    }
	resultat *= (1. / 4. / PI) * velocity;
    
	return resultat; /* [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] */
  }
}

/********************************************************************************************/
/********************************************************************************************/

void MIN_MED_MAX_loading(struct Structure_Propagation* pt_Propagation)
{
			// COEFFICIENTS RELATIFS AU MODELE DE DIFFUSION_PROPAGATION MAXIMUM 
	
	#ifdef MAX

		pt_Propagation->DIFFUSION_0_GV		 = 0.0765 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR;   												// [cm^{2} s^{-1}]
		pt_Propagation->PUISSANCE_COEFF_DIFF = 0.46;                                                										// [NO UNIT]
		pt_Propagation->E_DIFFUS         	 = 15.0;                                           												// [kpc]
		pt_Propagation->VENT_GALACTIQUE  	 = 5.0  * 1.0e5;                               													// [cm s^{-1}]
		pt_Propagation->V_ALFEN          	 = 117.6  * 1.0e5;                             													// [cm s^{-1}]

	#endif

			//COEFFICIENTS RELATIFS AU MODELE DE DIFFUSION_PROPAGATION MEDIUM.

	#ifdef MED

		pt_Propagation->DIFFUSION_0_GV		 = 0.0112 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR;     												// [cm^{2} s^{-1}]
		pt_Propagation->PUISSANCE_COEFF_DIFF = 0.7;                                                 										// [NO UNIT]
		pt_Propagation->E_DIFFUS         	 = 4.0;                                           												// [kpc]
		pt_Propagation->VENT_GALACTIQUE  	 = 12.0  * 1.0e5;                                												// [cm s^{-1}]              
		pt_Propagation->V_ALFEN          	 = 52.9  * 1.0e5;                                												// [cm s^{-1}]

	#endif

			// COEFFICIENTS RELATIFS AU MODELE DE DIFFUSION_PROPAGATION MINIMUM 

	#ifdef MIN

		pt_Propagation->DIFFUSION_0_GV       = 0.0016 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR;     												// [cm^{2} s^{-1}]
		pt_Propagation->PUISSANCE_COEFF_DIFF = 0.85;                                                 										// [NO UNIT]
		pt_Propagation->E_DIFFUS         	 = 1.0;                                           												// [kpc]
		pt_Propagation->VENT_GALACTIQUE  	 = 13.5  * 1.0e5;                                												// [cm s^{-1}]
		pt_Propagation->V_ALFEN          	 = 22.4  * 1.0e5;                                												// [cm s^{-1}]

	#endif
	
}

/********************************************************************************************/
/********************************************************************************************/

//	Nous imprimons les coefficients de diffusion_propagation choisis dans le calcul.

void print_propagation_parameters(struct Structure_Propagation* pt_Propagation)
{
	printf(" DELTA           = %.5e [NO UNIT]\n",pt_Propagation->PUISSANCE_COEFF_DIFF);
	printf(" DIFFUSION_0_GV  = %.5e [cm^{2} s^{-1}]\n",pt_Propagation->DIFFUSION_0_GV);
	printf(" E_DIFFUS        = %.5e [kpc]\n",pt_Propagation->E_DIFFUS);
	printf(" VENT_GALACTIQUE = %.5e [cm s^{-1}]\n",pt_Propagation->VENT_GALACTIQUE);
	printf(" V_ALFEN         = %.5e [cm s^{-1}]\n\n",pt_Propagation->V_ALFEN);
}

/********************************************************************************************/
/********************************************************************************************/
