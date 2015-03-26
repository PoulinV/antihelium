#include "ANTI_PROTON.h"

#define NRANSI

/********************************************************************************************/
/* Nous definissons a ce niveau les variables externes du programme ANTI_PROTON.c
* exclusivement.
*/
//double TABLE_Abar_i[DIM_TAB_PBAR+1][NDIM+1];

/********************************************************************************************/
/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Pbar->BESSEL_PBAR_SEC_Epbar_i
* en fonction de l'energie CINETIQUE des antiprotons T_pbar et des coefficients de BESSEL i.
*
*/
void calculation_BESSEL_PBAR_SECONDARY_Epbar_i(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium, struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation)
{
  long i_pbar,i,i_proton;
  double T_pbar,E_pbar,impulsion_pbar,v_pbar,K_pbar;
  double Si,Abar_i;
  double dlog_E_proton,E_proton;
  double impulsion_proton[DIM_TAB_PROTON+1],weight_SIMSPON[DIM_TAB_PROTON+1];
/*
* On remet a zero les tableaux pt_Pbar->BESSEL_PBAR_SEC_Epbar_i et TABLE_Abar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] = 0.0;
      pt_Pbar->TABLE_Abar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On s'occupe un peu des protons. Il serait bon en effet d'economiser le temps de calcul.
* Nous bouclons maintenant sur les protons afin d'evaluer E_proton, v_proton et aussi
* weight_SIMSPON, c'est-a-dire le poids de SIMPSON de chaque point.
*/
  dlog_E_proton = pow((E_PROTON_MAX/E_PROTON_MIN),(1./(double)DIM_TAB_PROTON)) - 1.;
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    E_proton = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));
    impulsion_proton[i_proton] = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
    
    if (i_proton==0 || i_proton==DIM_TAB_PROTON) {weight_SIMSPON[i_proton] = 1./3.;}
    else {weight_SIMSPON[i_proton] = (1. + (double)(i_proton % 2)) * 2. / 3.;}
  }
/*
* On remplit maintenant le tableau BESSEL_PBAR_SEC_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    impulsion_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));
    v_pbar         = CELERITY_LIGHT * impulsion_pbar / E_pbar;
    K_pbar         = K_space_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
    
    for (i=1;i<=NDIM;i++)
    {
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_pbar,2));
/*    Abar_i est exprime en [cm s^{-1}].  */
      Abar_i  = pt_Propagation->VENT_GALACTIQUE;
      Abar_i += 2.0*E_DISC*CM_PAR_KPC *
      ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar)
      - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) * v_pbar *
			(DENSITE_H_DISC + pow(4.,(2./3.))*1.0*DENSITE_HE_DISC));
      Abar_i += K_pbar * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);
      pt_Pbar->TABLE_Abar_i[i_pbar][i] = Abar_i;
	  
      
      for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
      {
/*
*       CONTRIBUTION H ON H
*/
        pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] += pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] *
        pt_Proton->BESSEL_COEF_Enuc_i[i_proton][i] * DENSITE_H_DISC *
        impulsion_proton[i_proton] * CELERITY_LIGHT *
        weight_SIMSPON[i_proton] * dlog_E_proton; /* [antiprotons GeV^{-1} s^{-1} cm^{-3}] */
/*
*       CONTRIBUTION H ON HE
*/
        pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] += pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton] *
        pt_Proton->BESSEL_COEF_Enuc_i[i_proton][i] * DENSITE_HE_DISC *
        impulsion_proton[i_proton] * CELERITY_LIGHT *
        weight_SIMSPON[i_proton] * dlog_E_proton; /* [antiprotons GeV^{-1} s^{-1} cm^{-3}] */
/*
*       CONTRIBUTION HE ON H
*/
        pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] += pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_proton] *
        pt_Helium->BESSEL_COEF_Enuc_i[i_proton][i] * DENSITE_H_DISC *
        impulsion_proton[i_proton] * CELERITY_LIGHT *
        weight_SIMSPON[i_proton] * dlog_E_proton; /* [antiprotons GeV^{-1} s^{-1} cm^{-3}] */
/*
*       CONTRIBUTION HE ON HE
*/
        pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] += pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_proton] *
        pt_Helium->BESSEL_COEF_Enuc_i[i_proton][i] * DENSITE_HE_DISC *
        impulsion_proton[i_proton] * CELERITY_LIGHT *
        weight_SIMSPON[i_proton] * dlog_E_proton; /* [antiprotons GeV^{-1} s^{-1} cm^{-3}] */
      }
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] *= 2. * E_DISC*CM_PAR_KPC / Abar_i;
/*    S'exprime maintenant en unites de [antiprotons cm^{-3} GeV^{-1}].
*/
	}
  }
  return;
}

/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Pbar->BESSEL_PBAR_TER_Epbar_i
* en fonction de l'energie CINETIQUE des antiprotons T_pbar et des coefficients de BESSEL i.
*
*/
void calculation_BESSEL_PBAR_TERTIARY_Epbar_i(double alpha_i[NDIM+1], struct Structure_Pbar* pt_Pbar)
{
  long i_pbar,i;
  double dlog_T_pbar,T_pbar,E_pbar,impulsion_pbar,v_pbar;
  static double S_inel_NOANN_fois_v_pbar[DIM_TAB_PBAR+1];
  double Abar_i,SUM;
/*
* On s'occupe un peu des antiprotons. Il serait bon en effet d'economiser le temps
* de calcul. Nous bouclons maintenant sur les antiprotons afin d'evaluer le produit
* sigma_inelastic_NOANN_pbarH_TAN_and_NG * v_pbar afin de le reutiliser dans les
* boucles suivantes.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar         = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar         = T_pbar + MASSE_PROTON;
    impulsion_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));
    v_pbar         = CELERITY_LIGHT * impulsion_pbar / E_pbar;
    
    S_inel_NOANN_fois_v_pbar[i_pbar] = sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)
    * v_pbar;
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TER_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TER_Epbar_i.
*/
  dlog_T_pbar = pow((T_PBAR_MAX/T_PBAR_MIN),(1./(double)DIM_TAB_PBAR)) - 1.;
  for (i=1;i<=NDIM;i++)
  {
    SUM = 0.0;
    for (i_pbar=DIM_TAB_PBAR;i_pbar>=0;i_pbar--)
    {
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].  */
	  
      SUM += dlog_T_pbar * S_inel_NOANN_fois_v_pbar[i_pbar] *
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i];
      
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = SUM -
      S_inel_NOANN_fois_v_pbar[i_pbar] * pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i];
/*    S'exprime en unites de [antiprotons GeV^{-1} s^{-1}].
*/
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] *=
      2. * E_DISC*CM_PAR_KPC * (DENSITE_H_DISC + pow(4.,(2./3.))*1.0*DENSITE_HE_DISC) / Abar_i;
/*    S'exprime maintenant en unites de [antiprotons cm^{-3} GeV^{-1}].
*/
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i
* en sommant les contributions primaires, secondaires et tertiaires en fonction de
* l'energie CINETIQUE des antiprotons T_pbar et des coefficients de BESSEL i.
*
*/
void calculation_BESSEL_PBAR_SUM_123_Epbar_i(struct Structure_Pbar* pt_Pbar)
{
  long i_pbar,i;
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=1;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] + pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
	}
  }
  
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX. We define
* x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* Then, we define u[i_pbar] = BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] for a given BESSEL
* order i. By making it discontinuous, we modify the differential equation to be
* solved into an algebraic relation [A] * [u] = [r]. The matrix [A] is tridiagonal
* so that inversion is a straightforward BUT SOMETIMES HAZARDOUS process !
*
*/
void calculation_BESSEL_PBAR_TOT_direct_inversion_A(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_pbar,i;
  double DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  static double vec_a[DIM_TAB_PBAR+1];
  static double vec_b[DIM_TAB_PBAR+1];
  static double vec_c[DIM_TAB_PBAR+1];
  static double vec_r[DIM_TAB_PBAR+1];
  static double vec_u[DIM_TAB_PBAR+1];
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en inversant une matrice tridiagonale.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  for (i=1;i<=NDIM;i++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      vec_r[i_pbar] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
/*    Nous chargeons les contributions secondaires et tertiaires de la production
*     d'antiprotons.
*/      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*     Nous remplissons les coefficients de la matrice tridiagonale a inverser.
*/
      if (i_pbar == 0)
      {
        vec_a[0]  = 0.0;
        vec_b[0]  = 1.0;
        vec_b[0] -= grand_C_cal * b_coeff[0] / DELTA_x;
        vec_b[0] += grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
        vec_c[0]  = grand_C_cal * b_coeff[1] / DELTA_x;
        vec_c[0] -= grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        vec_a[DIM_TAB_PBAR]  = 0.0;
        vec_b[DIM_TAB_PBAR]  = 1.0;
        vec_c[DIM_TAB_PBAR]  = 0.0;
      }
      else
      {
        vec_a[i_pbar]  = - grand_C_cal * b_coeff[i_pbar-1] / 2. / DELTA_x;
        vec_a[i_pbar] -= grand_C_cal * a_coeff[i_pbar] / pow(DELTA_x,2);
        vec_b[i_pbar]  = 1.0;
        vec_b[i_pbar] += grand_C_cal * (a_coeff[i_pbar] + a_coeff[i_pbar+1]) /
        pow(DELTA_x,2);
        vec_c[i_pbar]  = grand_C_cal * b_coeff[i_pbar+1] / 2. / DELTA_x;
        vec_c[i_pbar] -= grand_C_cal * a_coeff[i_pbar+1] / pow(DELTA_x,2);
      }
/*
      vec_b[i_pbar] = 1.0;
      vec_a[i_pbar] = 0.0;
      vec_c[i_pbar] = 0.0;
*/
    }
    inversion_tridiagonal(vec_a,vec_b,vec_c,vec_r,vec_u);
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = vec_u[i_pbar];
    }
  }
  return;
}

/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX. We define
* x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* Then, we define u[i_pbar] = BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] for a given BESSEL
* order i. By making it discontinuous, we modify the differential equation to be
* solved into an algebraic relation [A] * [u] = [r]. The matrix [A] is tridiagonal
* so that inversion is a straightforward BUT SOMETIMES HAZARDOUS process !
*
*/
void calculation_BESSEL_PBAR_TOT_direct_inversion_B(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_pbar,i;
  double DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  static double vec_a[DIM_TAB_PBAR+1];
  static double vec_b[DIM_TAB_PBAR+1];
  static double vec_c[DIM_TAB_PBAR+1];
  static double vec_r[DIM_TAB_PBAR+1];
  static double vec_u[DIM_TAB_PBAR+1];
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en inversant une matrice tridiagonale.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  for (i=1;i<=NDIM;i++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      vec_r[DIM_TAB_PBAR - i_pbar] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
/*    Nous chargeons les contributions secondaires et tertiaires de la production
*     d'antiprotons.
*/      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*     Nous remplissons les coefficients de la matrice tridiagonale a inverser.
*/
      if (i_pbar == 0)
      {
        vec_c[DIM_TAB_PBAR]  = 0.0;
        vec_b[DIM_TAB_PBAR]  = 1.0;
        vec_b[DIM_TAB_PBAR] -= grand_C_cal * b_coeff[0] / DELTA_x;
        vec_b[DIM_TAB_PBAR] += grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
        vec_a[DIM_TAB_PBAR]  = grand_C_cal * b_coeff[1] / DELTA_x;
        vec_a[DIM_TAB_PBAR] -= grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        vec_c[0]  = 0.0;
        vec_b[0]  = 1.0;
        vec_a[0]  = 0.0;
      }
      else
      {
        vec_c[DIM_TAB_PBAR - i_pbar]  = - grand_C_cal * b_coeff[i_pbar-1] / 2. / DELTA_x;
        vec_c[DIM_TAB_PBAR - i_pbar] -= grand_C_cal * a_coeff[i_pbar] / pow(DELTA_x,2);
        vec_b[DIM_TAB_PBAR - i_pbar]  = 1.0;
        vec_b[DIM_TAB_PBAR - i_pbar] += grand_C_cal *
        (a_coeff[i_pbar] + a_coeff[i_pbar+1]) / pow(DELTA_x,2);
        vec_a[DIM_TAB_PBAR - i_pbar]  = grand_C_cal * b_coeff[i_pbar+1] / 2. / DELTA_x;
        vec_a[DIM_TAB_PBAR - i_pbar] -= grand_C_cal * a_coeff[i_pbar+1] / pow(DELTA_x,2);
      }
/*
      vec_b[i_pbar] = 1.0;
      vec_a[i_pbar] = 0.0;
      vec_c[i_pbar] = 0.0;
*/
    }
    inversion_tridiagonal(vec_a,vec_b,vec_c,vec_r,vec_u);
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[DIM_TAB_PBAR - i_pbar][i] = vec_u[i_pbar];
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX. We define
* x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* Then, we define u[i_pbar] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] for a given BESSEL
* order i. By making it discontinuous, we modify the differential equation to be
* solved into an algebraic relation [A] * [u] = [r]. The matrix [A] is tridiagonal
* so that inversion is a straightforward BUT SOMETIMES HAZARDOUS process !
*
*/
void calculation_BESSEL_PBAR_TOT_direct_inversion_TD_NR(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_pbar,i;
  double DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  float *vec_a, *vec_b, *vec_c, *vec_r, *vec_u;
  vec_a = vector(1,DIM_TAB_PBAR+1);
  vec_b = vector(1,DIM_TAB_PBAR+1);
  vec_c = vector(1,DIM_TAB_PBAR+1);
  vec_r = vector(1,DIM_TAB_PBAR+1);
  vec_u = vector(1,DIM_TAB_PBAR+1);
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en inversant une matrice tridiagonale.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  for (i=1;i<=NDIM;i++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      vec_r[i_pbar + 1] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
/*    Nous chargeons les contributions secondaires et tertiaires de la production
*     d'antiprotons.
*/      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*     Nous remplissons les coefficients de la matrice tridiagonale a inverser.
*/
      if (i_pbar == 0)
      {
        vec_a[1]  = 0.0;
        vec_b[1]  = 1.0;
        vec_b[1] -= grand_C_cal * b_coeff[0] / DELTA_x;
        vec_b[1] += grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
        vec_c[1]  = grand_C_cal * b_coeff[1] / DELTA_x;
        vec_c[1] -= grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2);
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        vec_a[DIM_TAB_PBAR + 1]  = 0.0;
        vec_b[DIM_TAB_PBAR + 1]  = 1.0;
        vec_c[DIM_TAB_PBAR + 1]  = 0.0;
      }
      else
      {
        vec_a[i_pbar + 1]  = - grand_C_cal * b_coeff[i_pbar-1] / 2. / DELTA_x;
        vec_a[i_pbar + 1] -= grand_C_cal * a_coeff[i_pbar] / pow(DELTA_x,2);
        vec_b[i_pbar + 1]  = 1.0;
        vec_b[i_pbar + 1] += grand_C_cal * (a_coeff[i_pbar] + a_coeff[i_pbar+1]) /
        pow(DELTA_x,2);
        vec_c[i_pbar + 1]  = grand_C_cal * b_coeff[i_pbar+1] / 2. / DELTA_x;
        vec_c[i_pbar + 1] -= grand_C_cal * a_coeff[i_pbar+1] / pow(DELTA_x,2);
      }
/*
      vec_b[i_pbar + 1] = 1.0;
      vec_a[i_pbar + 1] = 0.0;
      vec_c[i_pbar + 1] = 0.0;
*/
    }
    tridag(vec_a,vec_b,vec_c,vec_r,vec_u,(DIM_TAB_PBAR + 1));
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = vec_u[i_pbar + 1];
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX. We define
* x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* Then, we define u[i_pbar] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] for a given BESSEL
* order i. By making it discontinuous, we modify the differential equation to be
* solved into an algebraic relation [A] * [u] = [r]. The matrix [A] is tridiagonal
* so that inversion is a straightforward BUT SOMETIMES HAZARDOUS process !
*
*/
void calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_pbar,j_pbar,i;
  double DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  float **matrix_A, **matrix_U;
  matrix_A = matrix(1,(DIM_TAB_PBAR+1),1,(DIM_TAB_PBAR+1));
  matrix_U = matrix(1,(DIM_TAB_PBAR+1),1,1);
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en inversant la matrice
* matrix_A[i_pbar + 1][j_pbar + 1] via la methode de GAUSS--JORDAN.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  for (i=1;i<=NDIM;i++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      matrix_U[i_pbar + 1][1] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
/*    Nous chargeons les contributions secondaires et tertiaires de la production
*     d'antiprotons.
*/      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*     Nous remplissons les coefficients de la matrice tridiagonale a inverser.
*/
      for (j_pbar=0;j_pbar<=DIM_TAB_PBAR;j_pbar++)
      {
        matrix_A[i_pbar + 1][j_pbar + 1] = 0.0;
      }
      if (i_pbar == 0)
      {
        matrix_A[1][1] = 1.0  -  (grand_C_cal * b_coeff[0] / DELTA_x)  +
        (grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2));
        matrix_A[1][2] = (grand_C_cal * b_coeff[1] / DELTA_x)  -
        (grand_C_cal * (a_coeff[1] - a_coeff[0]) / pow(DELTA_x,2));
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        matrix_A[DIM_TAB_PBAR + 1][DIM_TAB_PBAR + 1] = 1.0;
      }
      else
      {
        matrix_A[i_pbar + 1][i_pbar] = - (grand_C_cal*b_coeff[i_pbar-1]/2./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar]/pow(DELTA_x,2));
        matrix_A[i_pbar + 1][i_pbar + 1] = 1.0 +
        (grand_C_cal*(a_coeff[i_pbar] + a_coeff[i_pbar+1])/pow(DELTA_x,2));
        matrix_A[i_pbar + 1][i_pbar + 2] = (grand_C_cal*b_coeff[i_pbar+1]/2./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar+1]/pow(DELTA_x,2));
      }
    }
    gaussj(matrix_A,(DIM_TAB_PBAR + 1),matrix_U,1);
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = matrix_U[i_pbar + 1][1];
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX by evolving
* in time an initial distribution
*
* u^{0} \; = \;
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter} \;\; ,
*
* through the diffusion equation
*
* \frac{\partial u}{\partial t} \; + \; u \; + \;
* \left( \frac{2 h}{\bar{A}_{i} \, T} \right) \, \frac{\partial}{\partial x}
* \left\{ b u \, - \, a \frac{partial u}{\partial x} \right\} \; = \; 0 \;\; .
*
* u \equiv \bar{\cal{P}}_{i}(E) which may be written as
* u[i_pbar] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i]
*
* x = log(T/T_PBAR_MIN) and x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that
* x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* grand_C_cal = \left( \frac{2 h}{\bar{A}_{i} \, T} \right)
*
* We then evolve in time u(x,t) by using a CRANK & NICHOLSON scheme.
* There is always a little problem at the boundaries of the domain that is for
* T_PBAR_MIN and T_PBAR_MAX where the time evolution is partially explicit.
*
*/
void calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(double time_max,long n_time,struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_time,i_pbar,i;
  double DELTA_t,DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  static double vec_a[DIM_TAB_PBAR+1];
  static double vec_b[DIM_TAB_PBAR+1];
  static double vec_c[DIM_TAB_PBAR+1];
  static double vec_r[DIM_TAB_PBAR+1];
  static double vec_u[DIM_TAB_PBAR+1];
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en faisant evoluer dans le temps une
* distribution diffusant en energie et s'attenuant peu a peu.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  DELTA_t = time_max / (double)n_time;
  for (i=1;i<=NDIM;i++)
  {
/*
*   INITIALISATION = PHASE 1 !
******************************
*/
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
/*
*     On charge la distribution initiale u^{0}.
*/
      vec_u[i_pbar] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*   On calcule maintenant les coefficients de la matrice tridiagonale
*   d'evolution temporelle.
*/
      if (i_pbar == 0)
      {
        vec_a[0] = 0.0;
        
        vec_b[0] = 1.0 + (DELTA_t/2.) - (grand_C_cal*b_coeff[0]*DELTA_t/2./DELTA_x) +
        (grand_C_cal*(a_coeff[1] - a_coeff[0])*DELTA_t/2./pow(DELTA_x,2));
        
        vec_c[0] = (grand_C_cal*b_coeff[1]*DELTA_t/2./DELTA_x) -
        (grand_C_cal*(a_coeff[1] - a_coeff[0])*DELTA_t/2./pow(DELTA_x,2));
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        vec_a[DIM_TAB_PBAR] = 0.0;
        vec_b[DIM_TAB_PBAR] = 1.0 + (DELTA_t/2.);
        vec_c[DIM_TAB_PBAR] = 0.0;
      }
      else
      {
        vec_a[i_pbar] = - (grand_C_cal*b_coeff[i_pbar-1]*DELTA_t/4./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar]*DELTA_t/2./pow(DELTA_x,2));
        
        vec_b[i_pbar] = 1 + (DELTA_t/2.) +
        (grand_C_cal*(a_coeff[i_pbar] + a_coeff[i_pbar+1])*DELTA_t/2./pow(DELTA_x,2));
        
        vec_c[i_pbar] = (grand_C_cal*b_coeff[i_pbar+1]*DELTA_t/4./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar+1]*DELTA_t/2./pow(DELTA_x,2));
      }
    }
/*
*   TIME EVOLUTION = PHASE 2 !
******************************
*/
    for (i_time=1;i_time<=n_time;i_time++)
    {
      for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
      {
        pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] += DELTA_t * vec_u[i_pbar];
        
        Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*      Abar_i est exprime en [cm s^{-1}].
*/
        T_pbar = T_PBAR_MIN *
        pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
        grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
        
        if (i_pbar == 0)
        {
          vec_r[0] = (2.0 - vec_b[0])*vec_u[0]  -  vec_c[0]*vec_u[1];
        }
        else if (i_pbar == DIM_TAB_PBAR)
        {
          vec_r[DIM_TAB_PBAR] = (2.0 - vec_b[DIM_TAB_PBAR])*vec_u[DIM_TAB_PBAR];
        }
        else
        {
          vec_r[i_pbar] = - vec_a[i_pbar]*vec_u[i_pbar-1] +
          (2.0 - vec_b[i_pbar])*vec_u[i_pbar] - vec_c[i_pbar]*vec_u[i_pbar+1];
        }
      }
/*
*     On calcule maintenant le pas temporel suivant.
*/
      inversion_tridiagonal(vec_a,vec_b,vec_c,vec_r,vec_u);
    }
/*
* On passe au nouveau coefficient de BESSEL.
*/
  }
  return;
}

/********************************************************************************************/
/*
* We need to solve the differential equation that describes the energy behaviour
* of the BESSEL transforms \bar{\cal{P}}_{i} :
*
* \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
* b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
* \bar{\cal{P}}_{i} \right\} \; = \;
* \bar{A}_{i} \, \left\{
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
* \right\} \;\; .
*
* This is a mere diffusion equation. This routine tries to solve that equation
* on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX by evolving
* in time an initial distribution
*
* u^{0} \; = \;
* \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter} \;\; ,
*
* through the diffusion equation
*
* \frac{\partial u}{\partial t} \; + \; u \; + \;
* \left( \frac{2 h}{\bar{A}_{i} \, T} \right) \, \frac{\partial}{\partial x}
* \left\{ b u \, - \, a \frac{partial u}{\partial x} \right\} \; = \; 0 \;\; .
*
* u \equiv \bar{\cal{P}}_{i}(E) which may be written as
* u[i_pbar] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i]
*
* x = log(T/T_PBAR_MIN) and x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that
* x[i_pbar] = DELTA_x * (double)i_pbar
* with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR.
*
* grand_C_cal = \left( \frac{2 h}{\bar{A}_{i} \, T} \right)
*
* We then evolve in time u(x,t) by using a CRANK & NICHOLSON scheme.
* There is always a little problem at the boundaries of the domain that is for
* T_PBAR_MIN and T_PBAR_MAX where the time evolution is partially explicit.
*
*/
void calculation_BESSEL_PBAR_TOT_diffusion_soluce_B(double time_max,long n_time,struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation)
{
  long i_time,i_pbar,i;
  double DELTA_t,DELTA_x,T_pbar,E_pbar,Abar_i,grand_C_cal;
  static double a_coeff[DIM_TAB_PBAR+2],b_coeff[DIM_TAB_PBAR+1];

  static double vec_a[DIM_TAB_PBAR+1];
  static double vec_b[DIM_TAB_PBAR+1];
  static double vec_c[DIM_TAB_PBAR+1];
  static double vec_r[DIM_TAB_PBAR+1];
  static double vec_u[DIM_TAB_PBAR+1];
/*
* We compute the coefficients a_{j-1/2} = a_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR+1;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((-0.5+(double)i_pbar)/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    a_coeff[i_pbar] = D_energy_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation) / T_pbar;
/*  Ce coefficient s'exprime en unites de [GeV s^{-1}].
*/
  }
/*
* We compute the coefficients b_{j} = b_coeff[j] with j = i_pbar.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    T_pbar = T_PBAR_MIN *
    pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
    E_pbar = T_pbar + MASSE_PROTON;
    b_coeff[i_pbar] = b_energy_losses(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
/*  Ce coefficient s'exprime egalement en unites de [GeV s^{-1}].
*/
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_TOT_Epbar_i a partir
* de pt_Pbar->BESSEL_PBAR_PRI_Epbar_i, de pt_Pbar->BESSEL_PBAR_SEC_Epbar_i ainsi que
* de pt_Pbar->BESSEL_PBAR_TER_Epbar_i en faisant evoluer dans le temps une
* distribution diffusant en energie et s'attenuant peu a peu.
*/
  DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (double)DIM_TAB_PBAR; /* [sans unite] */
  DELTA_t = time_max / (double)n_time;
  for (i=1;i<=NDIM;i++)
  {
/*
*   INITIALISATION = PHASE 1 !
******************************
*/
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
/*
*     On charge la distribution initiale u^{0}.
*/
      vec_u[DIM_TAB_PBAR - i_pbar] =
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] +
      pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i];
      
      Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*    Abar_i est exprime en [cm s^{-1}].
*/
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
/*
*   On calcule maintenant les coefficients de la matrice tridiagonale
*   d'evolution temporelle.
*/
      if (i_pbar == 0)
      {
        vec_c[DIM_TAB_PBAR] = 0.0;
        
        vec_b[DIM_TAB_PBAR] = 1.0 + (DELTA_t/2.) -
        (grand_C_cal*b_coeff[0]*DELTA_t/2./DELTA_x) +
        (grand_C_cal*(a_coeff[1] - a_coeff[0])*DELTA_t/2./pow(DELTA_x,2));
        
        vec_a[DIM_TAB_PBAR] = (grand_C_cal*b_coeff[1]*DELTA_t/2./DELTA_x) -
        (grand_C_cal*(a_coeff[1] - a_coeff[0])*DELTA_t/2./pow(DELTA_x,2));
      }
      else if (i_pbar == DIM_TAB_PBAR)
      {
        vec_c[0] = 0.0;
        vec_b[0] = 1.0 + (DELTA_t/2.);
        vec_a[0] = 0.0;
      }
      else
      {
        vec_c[DIM_TAB_PBAR - i_pbar] =
        - (grand_C_cal*b_coeff[i_pbar-1]*DELTA_t/4./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar]*DELTA_t/2./pow(DELTA_x,2));
        
        vec_b[DIM_TAB_PBAR - i_pbar] = 1 + (DELTA_t/2.) +
        (grand_C_cal*(a_coeff[i_pbar] + a_coeff[i_pbar+1])*DELTA_t/2./pow(DELTA_x,2));
        
        vec_a[DIM_TAB_PBAR - i_pbar] =
        (grand_C_cal*b_coeff[i_pbar+1]*DELTA_t/4./DELTA_x) -
        (grand_C_cal*a_coeff[i_pbar+1]*DELTA_t/2./pow(DELTA_x,2));
      }
    }
/*
*   TIME EVOLUTION = PHASE 2 !
******************************
*/
    for (i_time=1;i_time<=n_time;i_time++)
    {
      for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
      {
        pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i] += DELTA_t * vec_u[DIM_TAB_PBAR - i_pbar];
        
        Abar_i = pt_Pbar->TABLE_Abar_i[i_pbar][i];
/*      Abar_i est exprime en [cm s^{-1}].
*/
        T_pbar = T_PBAR_MIN *
        pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
        grand_C_cal = 2. * E_DISC*CM_PAR_KPC / T_pbar / Abar_i;
        
        if (i_pbar == 0)
        {
          vec_r[DIM_TAB_PBAR] = (2.0 - vec_b[DIM_TAB_PBAR])*vec_u[DIM_TAB_PBAR] -
          vec_a[DIM_TAB_PBAR]*vec_u[DIM_TAB_PBAR-1];
        }
        else if (i_pbar == DIM_TAB_PBAR)
        {
          vec_r[0] = (2.0 - vec_b[0])*vec_u[0];
        }
        else
        {
          vec_r[i_pbar] = - vec_a[i_pbar]*vec_u[i_pbar-1] +
          (2.0 - vec_b[i_pbar])*vec_u[i_pbar] - vec_c[i_pbar]*vec_u[i_pbar+1];
        }
      }
/*
*     On calcule maintenant le pas temporel suivant.
*/
      inversion_tridiagonal(vec_a,vec_b,vec_c,vec_r,vec_u);
    }
/*
* On passe au nouveau coefficient de BESSEL.
*/
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
* This routine solves the matrix equation [A] * [u] = [r]
* in the case of a tridiagonal matrix [A]. Pivoting is not used.
* PROCEED AT YOUR OWN RISK !
*
*/
void inversion_tridiagonal(double a[DIM_TAB_PBAR+1],double b[DIM_TAB_PBAR+1],
double c[DIM_TAB_PBAR+1],double r[DIM_TAB_PBAR+1],double u[DIM_TAB_PBAR+1])
{
  extern FILE *probleme;
  long j;
  double bet,gam[DIM_TAB_PBAR+1];

  if (b[0] == 0.0)
  {
    fprintf(probleme,
    " PROBLEME DANS L'INVERSION DE LA MATRICE TRIDIAGONALE \n"
    " b[0] = %.5e \n",b[0]);
    return;
  }

  bet  = b[0];
  u[0] = r[0] / bet;
  for (j=1;j<=DIM_TAB_PBAR;j++)
  {
    gam[j] = c[j-1] / bet;
    bet = b[j] - a[j]*gam[j];
    if (bet == 0.0)
    {
      fprintf(probleme,
      " PROBLEME DANS L'INVERSION DE LA MATRICE TRIDIAGONALE \n"
      " j = %.ld , beta_j_j = %.5e \n",j,bet);
      return;
    }
    u[j] = (r[j] - a[j]*u[j-1]) / bet;
  }

  for (j=DIM_TAB_PBAR-1;j>=0;j--)
  {
    u[j] = u[j] - gam[j+1]*u[j+1];
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/

void PBAR_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PBAR+1])
{
	long i_pbar;
	
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
	    SPECTRUM[i_pbar] = 0.0;
	}
	
}


//	On remet a zero les tableaux Pbar.BESSEL_PBAR_SEC_Epbar_i et Pbar.BESSEL_PBAR_TER_Epbar_i.

void PBAR_BESSEL_TABLES_123_initialization(struct Structure_Pbar* pt_Pbar)
{
	long i_pbar,i;
	
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		for (i=0;i<=NDIM;i++)
		{
			pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
			pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] = 0.0;
			pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = 0.0;
		}
	}
}

/********************************************************************************************/
/********************************************************************************************/

void PBAR_IS_SPECTRUM_calculation(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1])
{
	long i_pbar,i;
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;
		for (i=1;i<=NDIM;i++)
		{
			pt_Pbar->BESSEL_PBARi[i] = pt_Pbar->BESSEL_PBAR_TOT_Epbar_i[i_pbar][i];
		}
		flux_antiproton_IS = GENERIC_FLUX_04(R_EARTH,0.,E_pbar_IS,MASSE_PROTON,1.,alpha_i,pt_Pbar->BESSEL_PBARi, pt_Propagation);
		//flux_antiproton_IS = GENERIC_FLUX(R_EARTH,0.,E_pbar_IS,MASSE_PROTON,1.,alpha_i,pt_Pbar->BESSEL_PBARi, &Propagation);

		
		PBAR_IS_SPECTRUM[i_pbar] = flux_antiproton_IS;	
	}
}

/********************************************************************************************/
/********************************************************************************************/

//	Modulation des spectres IS ===> TOA.
//	Cette fonction doit etre precedee de la fonction 'PBAR_IS_SPECTRUM_calculation' qui remplit le tableau PBAR_IS_SPECTRUM necessaire au calcule du flux d'antiprotons TOA.


void PBAR_TOA_SPECTRUM_calculation(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double PBAR_TOA_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_TOA[DIM_TAB_PBAR+1], struct Structure_Propagation* pt_Propagation)
{
	long i_pbar;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA,flux_proton_TOA;
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double flux_pbar, flux_pbar_TOA;
	


	pt_Propagation->PHI_FISK = fisk_potential;
	

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

//		Nous modulons maintenant les spectres PBAR obtenus.

		FFA_IS_to_TOA(1.,1.,pt_Propagation->PHI_FISK,E_pbar_IS,PBAR_IS_SPECTRUM[i_pbar],&E_pbar_TOA,&flux_pbar_TOA);
		
		if (E_pbar_TOA <= MASSE_PROTON)
		{
			T_PBAR_TOA[i_pbar]        = 0.0;
			PBAR_TOA_SPECTRUM[i_pbar] = 0.0;
			continue;
		}
		
		T_pbar_TOA = E_pbar_TOA - MASSE_PROTON;

//		Nous les stockons en memoire dans les tableaux RESULTS_T_PBAR_TOA[DIM_TAB_PBAR+1] et RESULTS_SPECTRUM_TOA_MIN_MED_MAX[DIM_TAB_PBAR+1];
	
		T_PBAR_TOA[i_pbar] = T_pbar_TOA;
		PBAR_TOA_SPECTRUM[i_pbar] = flux_pbar_TOA;										// [#pbar cm^{-3} sr^{-1} s^{-1} GeV^{-1}]
		
	}
}

/********************************************************************************************/
/********************************************************************************************/

void tertiary_component_effect_calculation(struct Structure_Pbar* pt_Pbar, double alpha_i[NDIM+1])
{
	long i_iteration;

	for (i_iteration=1;i_iteration<=10;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, pt_Pbar);
		calculation_BESSEL_PBAR_SUM_123_Epbar_i(pt_Pbar);
	}
}

/********************************************************************************************/
/********************************************************************************************/

void ELDR_effect_calculation(struct Structure_Propagation* pt_Propagation, struct Structure_Pbar* pt_Pbar, double alpha_i[NDIM+1])
{
	long i_iteration;

	for (i_iteration=1;i_iteration<=15;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, pt_Pbar);
		calculation_BESSEL_PBAR_TOT_direct_inversion_A(pt_Pbar, pt_Propagation);
		//calculation_BESSEL_PBAR_TOT_direct_inversion_B(&Pbar, &Propagation);
		//calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(&Pbar, &Propagation);
		//calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(15., 1200., &Pbar, &Propagation);
	}
}

/********************************************************************************************/
/********************************************************************************************/
	
#undef NRANSI
