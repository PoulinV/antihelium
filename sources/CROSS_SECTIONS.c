#include "CROSS_SECTIONS.h"

/*
* CE MODULE RENVOIE EN MILLIBARN LA SECTION EFFICACE DES PROTONS, DES ANTIPROTONS ET
* DES ANTIDEUTERONS SUR L'HYDROGENE, SECTION EFFICACE PARAMETREE PAR LES DONNEES DU CERN
* (http://pdg.lbl.gov)
*
* THE PARAMETRIZATION IS
* sigma(p) = A + B*p^n + C*ln^2(p) + D*ln(p)
* WHERE SIGMA IS IN MILLIBARN AND p IS IN GEV/C. THE BEST-FIT COEFFICIENTS A, B, C AND D,
* AND THE EXPONENT n ARE TABULATED BELOW. ALSO GIVEN IS THE RANGE OF MOMENTUM OVER WHICH
* THE FIT WAS DONE.
*/
/********************************************************************************************/
/********************************************************************************************/
/*
* sigma_total_pH est la section efficace totale de collision d'un proton cosmique
* avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
*
*/
double sigma_total_pH(double E_proton)
{
  double p;
  double resultat;

//p = sqrt(pow(E_proton,2.) - pow(m_p*uma_Gev ,2.));
  p = sqrt(pow(E_proton,2.) - pow(MASSE_PROTON,2.));

  if (p<=p_seuil_pp_tot)
  {
    resultat = P1_pp_tot + P2_pp_tot*log(p) + P3_pp_tot*pow(log(p),2.) +
               P4_pp_tot*pow(log(p),3.) + P5_pp_tot*pow(log(p),4.) +
               P6_pp_tot*pow(log(p),5.) + P7_pp_tot*pow(log(p),6.);
  }
  else
  {
    resultat = A_pp_tot + B_pp_tot*pow(p,n_pp_tot) + C_pp_tot*pow(log(p),2.) +
               D_pp_tot*log(p);
  }
  return (resultat*mb_cm2);
}

/********************************************************************************************/
/*
* sigma_inelastic_pH_TAN_and_NG est la section efficace inelastique de
* collision d'un proton cosmique avec un atome d'hydrogene au repos.
*
* Cette section efficace est exprimee en [cm^{2}].
* TAN et NG la notent \sigma^{i}_{p}.
*
*/
double sigma_inelastic_pH_TAN_and_NG(double E_proton)
{
  double EK_proton,U,Cp;
  double resultat;

  EK_proton = E_proton - MASSE_PROTON;
  resultat = 0.0;
  if (EK_proton<=0.0)
  {
    return resultat;
  }
  else if (EK_proton<0.3)
  {
    return resultat;
/*
*   JE NE CROIS PAS A LEUR VALEUR LORSQUE L'ENERGIE CINETIQUE EST INFERIEURE
*   A 300 MEV ! TO BE CHECKED !
*/
  }
  else if (EK_proton<3.0)
  {
    U = log(E_proton/200.0);
    resultat  = 1. + 0.0273*U;
    resultat *= 32.2e-27;
    Cp = 17.9 + 13.8*log(EK_proton) + 4.41*log(EK_proton)*log(EK_proton);
    resultat /= 1. + 2.62e-3*pow(EK_proton,-Cp);
    return resultat;
  }
  else if (E_proton<200.0)
  {
    U = log(E_proton/200.0);
    resultat = 1. + 0.0273*U;
    resultat *= 32.2e-27;
    return resultat;
  }
  else
  {
    U = log(E_proton/200.0);
    resultat = 1. + 0.0273*U + 0.01*U*U;
    resultat *= 32.2e-27;
    return resultat;
  }
}

/********************************************************************************************/
/********************************************************************************************/
/*
* sigma_total_pbarH est la section efficace totale de collision d'un antiproton cosmique
* avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
* Remarquons que :
* sigma_total_pbarH = sigma_elastic_pbarH + sigma_inelastic_pbarH.
*
*/
double sigma_total_pbarH(double E_pbar)
{
  double p;
  double resultat;

  p = sqrt(pow(E_pbar,2.) - pow(m_p*uma_Gev,2.));
  if (p<=p_seuil_pbarp_tot)
  {
    resultat = P1_pbarp_tot + P2_pbarp_tot*log(p) + P3_pbarp_tot*pow(log(p),2.) +
               P4_pbarp_tot*pow(log(p),3.) + P5_pbarp_tot*pow(log(p),4.) +
               P6_pbarp_tot*pow(log(p),5.) + P7_pbarp_tot*pow(log(p),6.);
  }
  else
  {
    resultat = A_pbarp_tot + B_pbarp_tot*pow(p,n_pbarp_tot) + C_pbarp_tot*pow(log(p),2.) +
               D_pbarp_tot*log(p);
  }
  return (resultat*mb_cm2);
}

/********************************************************************************************/
/*
* sigma_elastic_pbarH est la section efficace elastique de collision d'un antiproton
* cosmique avec un atome d'hydrogene au repos.
* Cette section efficace est exprimee en [cm^{2}]. Elle correspond au processus
*
* pbar + p -----> pbar + p.
*
* En diffusant elastiquement sur des protons interstellaires au repos, les antiprotons
* ne perdent pas d'energie dans la mesure ou l'angle de deviation est estime tres faible.
*
* Cependant, ils peuvent entrer en interaction avec des protons au repos en les excitant
* de maniere inelastique sans toutefois s'annihiler avec eux. Cette interaction inelastique
* mais non-annihilante engendre une forte deperdition d'energie pour les antiprotons.
*
*/
double sigma_elastic_pbarH(double E_pbar)
{
  double p;
  double resultat;

  p = sqrt(pow(E_pbar,2.) - pow(m_p*uma_Gev,2.));
  if (p<=p_seuil_pbarp_el)
  {
    resultat = P1_pbarp_el + P2_pbarp_el*log(p) + P3_pbarp_el*pow(log(p),2.) +
               P4_pbarp_el*pow(log(p),3.) + P5_pbarp_el*pow(log(p),4.) +
               P6_pbarp_el*pow(log(p),5.) + P7_pbarp_el*pow(log(p),6.);
  }
  else
  {
    resultat = A_pbarp_el + B_pbarp_el*pow(p,n_pbarp_el) + C_pbarp_el*pow(log(p),2.) +
               D_pbarp_el*log(p);
  }
  return (resultat*mb_cm2);
}

/********************************************************************************************/
/*
* sigma_inelastic_pbarH_TAN_and_NG est la section efficace inelastique de
* collision d'un antiproton cosmique avec un atome d'hydrogene au repos.
* Cette section efficace est exprimee en [cm^{2}]. Remarquons que :
*
* sigma_inelastic_pbarH_TAN_and_NG =
* sigma_ANN_pbarH_TAN_and_NG + sigma_inelastic_NOANN_pbarH_TAN_and_NG.
*
* TAN et NG notent cette relation sous la forme
* \sigma^{inel}_{pbar} \; = \; \sigma^{an}_{pbar} + \sigma^{i}_{pbar} \;\; .
*
* sigma_inelastic_pbarH_TAN_and_NG       est designee par \sigma^{inel}_{pbar}
* sigma_ANN_pbarH_TAN_and_NG             est designee par \sigma^{an}_{pbar}
* sigma_inelastic_NOANN_pbarH_TAN_and_NG est designee par \sigma^{i}_{pbar}
*
* A basse energie, l'annihilation domine. Par contre, a haute energie, l'interaction
* inelastique d'un antiproton avec un proton devient non-annihilante.
*
* A la limite des hautes energies, on a de surcroit le comportement asymptotique
* \sigma^{i}_{pbar} \simeq \sigma^{i}_{p} \;\; .
*
*/
double sigma_inelastic_pbarH_TAN_and_NG(double E_pbar)
{
  extern FILE *probleme;
  double EK_pbar;
  double resultat;

  EK_pbar = E_pbar - MASSE_PROTON;
  resultat = 0.0;
  if (EK_pbar<=0.0)
  {
    return resultat;
  }
  else if (EK_pbar<0.05)
  {
    fprintf(probleme,
    " WARNING : sigma_inelastic_pbarH_TAN_and_NG NOT DEFINED ! \n");
    return resultat;
  }
  else
  {
    resultat  = 1. + 0.584*pow(EK_pbar,-0.115) + 0.856*pow(EK_pbar,-0.566);
    resultat *= 24.7e-27;
    return resultat;
  }
}

/********************************************************************************************/
/*
* sigma_inelastic_NOANN_pbarH_TAN_and_NG est la section efficace inelastique
* MAIS NON-ANNIHILANTE d'un antiproton cosmique avec un atome d'hydrogene au
* repos. Cette section efficace est exprimee en [cm^{2}].
*
* En diffusant inelastiquement sur des protons interstellaires au repos,
* les antiprotons peuvent parfois ne pas s'annihiler. Ils perdent alors
* une grande quantite d'energie en rechappant de l'annihilation, energie
* qu'ils transferent aux cibles excitees.
*
*/
double sigma_inelastic_NOANN_pbarH_TAN_and_NG(double E_pbar)
{
  extern FILE *probleme;
  double EK_pbar,sigma_ANN_pbarH_TAN_and_NG;
  double resultat;
/*
* Nous calculons tout d'abord la section efficace d'annihilation
* sigma_ANN_pbarH_TAN_and_NG.
*/
  EK_pbar = E_pbar - MASSE_PROTON;
  if (EK_pbar<=0.0)
  {
    sigma_ANN_pbarH_TAN_and_NG = 0.0;
  }
  else if (EK_pbar<0.05)
  {
    fprintf(probleme,
    " WARNING : sigma_ANN_pbarH_TAN_and_NG NOT DEFINED ! \n");
    sigma_ANN_pbarH_TAN_and_NG = 0.0;
  }
  else
  {
    sigma_ANN_pbarH_TAN_and_NG  =
    1. + 0.0115*pow(EK_pbar,-0.774) - 0.948*pow(EK_pbar,0.0151);
    sigma_ANN_pbarH_TAN_and_NG *= 661.0e-27;
  }

  resultat = 0.0;
  if (EK_pbar<=0.0)
  {
    return resultat;
  }
  else if (EK_pbar<0.05)
  {
    fprintf(probleme,
    " WARNING : sigma_inelastic_NOANN_pbarH_TAN_and_NG NOT DEFINED ! \n");
    return resultat;
  }
  else if (EK_pbar<=13.3)
  {
    resultat = sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_ANN_pbarH_TAN_and_NG;
    return resultat;
  }
  else
  {
    resultat = sigma_inelastic_pH_TAN_and_NG(E_pbar);
    return resultat;
  }
}

/********************************************************************************************/
/********************************************************************************************/
/*
* sigma_total_dbarH est la section efficace totale de collision d'un antideuteron cosmique
* avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
* Voir Physical Review D, 45 (1992) page III.83 et III.98. Voir plus recemment The
* Review of Particles Physics, The European Physical Journal C, 208 (1998).
*
* Par conjugaison de charge, elle est egale a la section efficace d'interaction
* d'un antiproton sur un cible de deuterium.
*
*/
double sigma_total_dbarH(double E_dbar)
{
  double p;
  double resultat;

  resultat = 0.0;

  p = sqrt(pow(E_dbar,2.) - pow(MASSE_DEUT,2.)) / 2.; /* Impulsion du proton */

  if (p<=p_seuil_min_pbard_tot)
  {
    resultat = P1_pbard_min_tot + P2_pbard_min_tot*(p) + P3_pbard_min_tot*pow((p),2.) +
               P4_pbard_min_tot*pow((p),3.) + P5_pbard_min_tot*pow((p),4.) +
               P6_pbard_min_tot*pow((p),5.) + P7_pbard_min_tot*pow((p),6.) +
               P8_pbard_min_tot*pow((p),7.) + P9_pbard_min_tot*pow((p),8.);
  }
  else if ((p_seuil_min_pbard_tot<p) && (p<=p_seuil_max_pbard_tot))
  {
    resultat = exp(P1_pbard_max_tot + P2_pbard_max_tot*log(p));
  }
  else if (p_seuil_max_pbard_tot<p)
  {
    resultat = A_pbard_tot + B_pbard_tot*pow(p,n_pbard_tot) + C_pbard_tot*pow(log(p),2.) +
               D_pbard_tot*log(p);
  }

  return (resultat*mb_cm2);
}

/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton + proton -----> antiproton + X
* consideree dans le referentiel du centre de masse CMF.
* Cette LI section efficace est exprimee en unites de [millibarns GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_ECM(double E_CMF,double pLstar,double pTstar)
{
  double Q,Estar,E_MAX_star,XR,f_XR,A,B;
  double XR_min,DXR,LE_sur_RS;
  double resultat;

  Q = E_CMF - 4.*MASSE_PROTON;
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) /
  (2.*E_CMF);
  Estar = sqrt(pow(MASSE_PROTON,2) + pow(pLstar,2) + pow(pTstar,2));
  XR = Estar / E_MAX_star;
/*
* SOME CHECKS ARE IN ORDER.
*/
  resultat = 0.0;
  if (Q<=0.0) goto LA_FIN;
  if (XR>=1.0) goto LA_FIN;

/*
* f_XR is expressed in [millibarns GeV^{-2}].
*/
  f_XR = 2.10 * pow((1.-XR),7.80);
  if (XR<=0.5)
  {
    f_XR += 3.34 * exp(-17.6*XR);
  }
  A = 3.95 * exp(-2.76*XR);                          /* [GeV^{-1}] */
  B = 40.5 * exp(-3.21*XR) * pow(XR,2.13);           /* [GeV^{-2}] */
  resultat = f_XR * exp(-A*pTstar-B*pTstar*pTstar);  /* [millibarns GeV^{-2}] */

  XR_min = sqrt(pow(pTstar,2) + pow(MASSE_PROTON,2)) / E_MAX_star;
  DXR = XR - XR_min;
  LE_sur_RS = (6.25e-3) * (exp(-0.592*Q) + 493.*exp(-5.40*Q)) *
  (exp(6.08 + 2.57*DXR + 7.95*DXR*DXR) - 1.) * exp(3.00*DXR*(3.09-Q));
  LE_sur_RS += 1.;
  resultat *= LE_sur_RS; /* [millibarns GeV^{-2}] */

LA_FIN :
  return resultat;
}

/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton + proton -----> antiproton + X
* consideree dans le referentiel du laboratoire dans lequel un des protons est au repos.
* Cette LI section efficace est exprimee en unites de [millibarns GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_LAB(double E_proton,double pL,double pT,
double *E_CMF,double *pLstar,double *pTstar)
{
  double gamma,beta,E_pbar,resultat;

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_pbar = sqrt(pow(MASSE_PROTON,2) + pow(pL,2) + pow(pT,2));
  *E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  *pLstar = gamma * (pL - beta*E_pbar);
  *pTstar = pT;

  resultat = E_d3S_on_d3P_PBAR_ECM(*E_CMF,*pLstar,*pTstar); /* [millibarns GeV^{-2}] */
  return resultat;
}

/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + proton (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee de maniere DIRECTE.
*
*/
double dSpbar_sur_dEpbar_DIRECTE(double E_proton,double E_pbar,long n_step_theta)
{
  double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long i_theta;
  double cos_theta,sin_theta,cos_theta_min,cos_theta_max,dcos_theta,resultat;

  resultat = 0.0;

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON)      goto LA_FIN;
  if (E_pbar>=(E_proton-2.*MASSE_PROTON)) goto LA_FIN;
  if (E_pbar<=MASSE_PROTON)         goto LA_FIN;
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma      = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta       = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) /
  (2.*E_CMF);

  cos_theta_max = 1.0;
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
  if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
  if (cos_theta_min>=+1.0) {goto LA_FIN;}

  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step_theta;
  cos_theta  = cos_theta_min - dcos_theta/2.;
  for (i_theta=1;i_theta<=n_step_theta;i_theta++)
  {
    cos_theta += dcos_theta;
    sin_theta  = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    resultat += E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar) *
    dcos_theta;
  }
  resultat *= 2.* PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de TAN et NG ne concerne que la production d'antiprotons
* lors d'une collision proton sur proton. Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */

LA_FIN :
  return resultat;
}

/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + proton (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON(double E_proton,double E_pbar,long n_step_theta)
{
  double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long i_theta;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
  double h,x1,x2,x3,f1,f2,f3,resultat;

  resultat = 0.0;

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON)      goto LA_FIN;
  if (E_pbar>=(E_proton-2.*MASSE_PROTON)) goto LA_FIN;
  if (E_pbar<=MASSE_PROTON)         goto LA_FIN;
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) /
  (2.*E_CMF);

  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
  if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
  if (cos_theta_min>=+1.0) {goto LA_FIN;}

  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step_theta;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

  for (i_theta=1;i_theta<=n_step_theta;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de TAN et NG ne concerne que la production d'antiprotons
* lors d'une collision proton sur proton. Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */

LA_FIN :
  return resultat;
}

/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module ecrit dans le fichier FILE_NAME_H_ON_H le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs de la section efficace differentielle du processus
* proton (E_proton) + proton (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'energie TOTALE des protons E_proton varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_H_ON_H_write_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei, double Z_nuclei, double M_nuclei)
{
  long i_proton,i_pbar;
  double E_proton,T_pbar,E_pbar;
  FILE *p_to_pbar_file;
/*
* On remet a zero le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton].
*/
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = 0.0;
    }
  }
/*
* Et maintenant, on le remplit !
*/
  if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_H_ON_H,"w");
  else if(A_nuclei == 2){
    p_to_pbar_file = fopen(FILE_NAME_DE_H_ON_H,"w");
  }
  else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_H_ON_H,"w");

  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    E_proton = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));
	  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      	T_pbar = T_PBAR_MIN*A_nuclei *
      	pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      	E_pbar = T_pbar + A_nuclei * MASSE_PROTON;


    if(A_nuclei == 1){
      #ifdef	DSPBAR_SUR_DEPBAR_H_on_H_unknown
        pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = dSpbar_sur_dEpbar_SIMPSON(E_proton,E_pbar,1000);
      #elif defined DSPBAR_SUR_DEPBAR_H_on_H_Duperray
        pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = dSpbar_sur_dEpbar_SIMPSON_H_ON_H_duperray(E_proton,E_pbar,1000);
      #elif defined DSPBAR_SUR_DEPBAR_H_on_H_high_energy_TAN_NG_mass_T_included
        pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = dSpbar_sur_dEpbar_SIMPSON_H_ON_H_tan_ng_mass_T(E_proton,E_pbar,1000);
      #elif defined DSPBAR_SUR_DEPBAR_H_on_H_MDGS_F12
        pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = dSpbar_sur_dEpbar_SIMPSON_H_ON_H_MDGS_F12(E_proton,E_pbar,1000);
      #else
        printf("ERROR : function 'DSPBAR_SUR_DEPBAR_H_ON_H_write_file' \n You must to specify one parametrization for H on H reaction! \n");
        exit(0);
      #endif
    }
    else if(A_nuclei == 2 || A_nuclei == 3){
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = dS3He_sur_sE3He_H_on_H_MDGS_F12_coalescence(E_proton,E_pbar,pt_Cross_Section->P_coal,1000,A_nuclei,Z_nuclei,M_nuclei);

    }
	  fprintf(p_to_pbar_file," %.6e \n",pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]);

	  // printf(" E_proton = %.5e -- T_pbar = %.5e -- pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H  = %.6e \n",
			// E_proton,T_pbar,pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/*
* Ce module nettoie -- remet a zero -- tous les elements des fichiers ou seront
* stockees les valeurs de la section efficace differentielle de production des
* antiprotons \frac{d \sigma_{\pbar}}{d E_{\pbar}} au cours de la reaction generique
*
* {P ou ALPHA} + {H ou HE}_{milieu interstellaire au repos} -----> PBAR + X
*
* Liste des tableaux remis a zero :
* - pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H  [i_pbar][i_proton]
* - pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE [i_pbar][i_proton]
* - pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H [i_pbar][i_nucleon]
* - pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon]
*
*/
void CLEANING_ALL_THE_DSPBAR_SUR_DEPBAR(struct Structure_Cross_Section* pt_Cross_Section)
{
  long i_proton,i_pbar;
/*
* On remet a zero les tableaux pt_Cross_Section->DSPBAR_SUR_DEPBAR_A_ON_B[i_pbar][i_proton].
*/
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton] = 0.0;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton] = 0.0;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_proton] = 0.0;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_proton] = 0.0;
    }
  }
  return;
}

/********************************************************************************************/
/*
* Ce module lit dans le fichier FILE_NAME_H_ON_H le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs -- exprimees en [cm^{2} GeV^{-1}] -- de la section
* efficace \frac{d \sigma_{\pbar}}{d E_{\pbar}} du processus
*
* PROTON INCIDENT + HYDROGENE AU REPOS -----> PBAR + X.
*
* L'energie TOTALE des protons E_proton varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_H_ON_H_read_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_proton,i_pbar;
  FILE *p_to_pbar_file;
/*
* On remplit le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton].
*/
  if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_H_ON_H,"r");
  else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_H_ON_H,"r");
  else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_H_ON_H,"r");
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      fscanf(p_to_pbar_file, " %lf ",&pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}



/********************************************************************************************/
/*
* Ce module lit dans le fichier FILE_NAME_H_ON_HE le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs -- exprimees en [cm^{2} GeV^{-1}] -- de la section
* efficace \frac{d \sigma_{\pbar}}{d E_{\pbar}} du processus
*
* PROTON INCIDENT + HELIUM AU REPOS -----> PBAR + X.
*
* L'energie TOTALE des protons E_proton varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_H_ON_HE_read_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_proton,i_pbar;
  FILE *p_to_pbar_file;
/*
* On remplit le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton].
*/
if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_H_ON_HE,"r");
else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_H_ON_HE,"r");
else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_H_ON_HE,"r");
for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      fscanf(p_to_pbar_file, " %lf ",&pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/*
* Ce module lit dans le fichier FILE_NAME_HE_ON_H le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs -- exprimees en [cm^{2} GeV^{-1}] -- de la section
* efficace \frac{d \sigma_{\pbar}}{d E_{\pbar}} du processus
*
* ALPHA INCIDENT + HYDROGENE AU REPOS -----> PBAR + X.
*
* L'energie TOTALE PAR NUCLEON E_nucleon des particules alpha incidentes varie
* de E_PROTON_MIN a E_PROTON_MAX en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_HE_ON_H_read_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_nucleon,i_pbar;
  FILE *p_to_pbar_file;
/*
* On remplit le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon].
*/
if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_HE_ON_H,"r");
else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_HE_ON_H,"r");
else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_HE_ON_H,"r");

  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      fscanf(p_to_pbar_file, " %lf ",&pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/*
* Ce module lit dans le fichier FILE_NAME_HE_ON_HE le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs -- exprimees en [cm^{2} GeV^{-1}] -- de la section
* efficace \frac{d \sigma_{\pbar}}{d E_{\pbar}} du processus
*
* ALPHA INCIDENT + HELIUM AU REPOS -----> PBAR + X.
*
* L'energie TOTALE PAR NUCLEON E_nucleon des particules alpha incidentes varie
* de E_PROTON_MIN a E_PROTON_MAX en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_HE_ON_HE_read_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_nucleon,i_pbar;
  FILE *p_to_pbar_file;
/*
* On remplit le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon].
*/
if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_HE_ON_HE,"r");
else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_HE_ON_HE,"r");
else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_HE_ON_HE,"r");

  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      fscanf(p_to_pbar_file, " %lf ",&pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/* CODE DE LAURENT DEROME POUR LES CALCULS DE LA REACTION H ON NUCLEUS WITH ATOMIC NUMBER A
*
*  retourne la multiplicite invariante 1/sig_reac * Edsig/d3p : p + A ---> ap + X
*  en unites de [GeV^{-2}] avec
*
*  At : Target atomic number
*  p1 : projectile momentum (GeV/c) (target reference system)
*  y  : produced particle rapidity  (target reference system) y = atanh(pl/E) = atanh(beta_CM) + atanh(pl_star /E_star)
*  mt : produced particle transverse mass (GeV/c2) = sqrt(m^2+pt^2)
*
*  #define PMASS 0.93827231 en unites de [GeV]
*
*/
double GetInvarMul_pAapX(double At,double p1,double y,double mt)
{
  double Ecm,pt,xp,xt,sig,phi;
  double E,E1,p,ptmax;
  double s,Betacm;
  double m = PMASS;
  double omega =  2.81484;
  double C1,C2,C3,C4,C5,C6,mu2,b0,D1,bb;
  double EcmMax;
  double par[] = {0.009205, 2.269, 3.708, 3.360, 2.485, 0.06394, 10.28, 0.1699, 0.4812, 0.1824};

  m = PMASS;
  omega =3.*m;
  C1 = par[0]; //   C1 =0.009205;
  C2 = par[1]; //   C2 =2.269;
  C3 = par[2]; //   C3 =3.708;
  C4 = par[3]; //   C4 =3.360;
  C5 = par[4]; //   C5 =2.485;
  C6 = par[5]; //   C6 =0.06394;
  D1 = par[6]; //   D1 = 10.28;
  b0 = par[7]; //   b0 =0.1699;
  bb = par[8]; //   bb =0.4812;
  mu2= par[9]; //   mu2=0.1824;
  sig = 0;

  E1 = sqrt(p1*p1+PMASS*PMASS); /*E1 = energie totale par nucleon */

  s = 2*PMASS*PMASS+2*PMASS*E1;
  Betacm = p1/(E1+PMASS);


  EcmMax =  (s + m*m-omega*omega)/(2.0*sqrt(s));
  //    EcmMax = Emaxferm;
  Ecm = mt*cosh(y-atanh(Betacm));
  xp  = Ecm/EcmMax;

  if(mt<m) return 0;
  pt  = sqrt(mt*mt-m*m);
  E   = mt*cosh(y);
  p   = sqrt(E*E - m*m);

  ptmax = pow(EcmMax/cosh(y-atanh(Betacm)),2.) - m*m;
  if(ptmax<0) return 0;
  ptmax = sqrt(ptmax);
  xt = pt/ptmax;

  if(xp<1.0 && xt>=0 && xt<1 && xp > 0) {
    phi = pow(sqrt(s),bb)*C1*exp(-C4*pt) + C6*exp(-C5*pow(pt,2.))/pow(sqrt(s),mu2);
    sig = pow(At,b0*log(sqrt(s)/D1)*pt)* pow((1.0-xp),C2*0.5*log(s))*exp(-C3*xp)*phi;
  }
  else sig=0;
  return sig;
}

/********************************************************************************************/
/* CODE DE LAURENT DEROME POUR LES CALCULS DE LA REACTION H ON NUCLEUS WITH ATOMIC NUMBER A
*
*  Get the total reaction cross section (mb)
*  At : Target atomic number
*  Pproj : projectile momentum (GeV/c) (target reference system)
*
*/
double GetReactionCrossSection(double At,double Pproj)
{
  if(At>1)
	{
/*
*   Section efficace totale de reaction p + Noyau avec dependance en masse A et energie labo Elab
*   selon Letaw (Ap.J.Suppl.,51,271) en unites de [millibarns].
*
*/
    double sig0=45.*pow(At,0.695)*(1.+0.016*sin(5.3-2.63*log(At)));

    double Tproj = sqrt(Pproj*Pproj+PMASS*PMASS) - PMASS;
    return sig0*(1.-0.62*exp(-Tproj*1000/200.)*sin(10.9*pow(Tproj*1000,-0.28)));
  }
	if(At==1)
	{
/*
*   from data @ 100 GeV (PDG http://pdg.lbl.gov/2004/hadronic-xsections/hadron.html)
*/
    return 38.46-7;
  }
	return 0;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton (E_proton) + helium (repos) -----> antiproton (E_pbar) + X.
* consideree dans le referentiel du laboratoire dans lequel le noyau d'helium est au repos.
* Cette LI section efficace est exprimee en unites de [millibarns GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_H_ON_HE_LAB(double E_proton,double pL,double pT,double A_nuclei,double P_coal)
{
  double P_proton,E_pbar,mt_pbar,y_pbar,resultat,Ed3S_on_d3P_pbar,Ed3S_on_d3P_nbar,B;

  P_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
	E_pbar   = sqrt(pow(MASSE_PROTON,2) + pT*pT + pL*pL);
	mt_pbar  = sqrt(pow(MASSE_PROTON,2) + pT*pT);
	y_pbar   = atanh(pL/E_pbar);
  Ed3S_on_d3P_pbar = GetInvarMul_pAapX(4.0,P_proton,y_pbar,mt_pbar) * GetReactionCrossSection(4.0,P_proton);
  if(A_nuclei == 1) resultat = Ed3S_on_d3P_pbar;
  if(A_nuclei == 2){
    B = MASSE_DEUT/(MASSE_PROTON*1e3*MASS_NEUTRON)*pow(4*PI/3*pow(P_coal,3)/8,2); //MeV^4
    B *= 1e-12; // conversion en GeV^4
    B /= 2*sigma_total_pH(E_proton)/mb_cm2; //convert cross section in mb, use a simple scaling factor 2 to be updated
    resultat = B*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*1.3;
  }
  if(A_nuclei == 3){
    B = MASS_HELIUM_TROIS/(MASSE_PROTON*1e3*MASSE_PROTON*1e3*MASS_NEUTRON)*pow(4*PI/3*pow(P_coal,3)/8,2); //MeV^4
    B *= 1e-12; //GeV^4
    B /= pow(2*sigma_total_pH(E_proton)/mb_cm2,2); //GeV^4 Simple scaling factor 2 to be updated
    resultat = B*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*1.3;
  }
  return resultat; /* [millibarns GeV^{-2}] */
}

/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + helium (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON_H_ON_HE(double E_proton,double E_pbar,long n_step,double A_nuclei,double P_coal)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_proton-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON){return resultat;}

  if (E_pbar<=MASSE_PROTON){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		fprintf(probleme," more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
	if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de LAURENT DEROME ne concerne que la production d'antiprotons
* lors d'une collision proton sur helium (He). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte. De meme pour les antiTritium si l'on calcule le flux d'antiHe3.
*/
  if(A_nuclei == 1 || A_nuclei == 3)resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE L'ANTIPROTON FINAL.
*/
  // if(A_nuclei == 1){
    mt_min = MASSE_PROTON;
    mt_max = pow(MASSE_PROTON,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
  // }
  // else if(A_nuclei == 2){
  //   mt_min = MASSE_DEUT;
  //   mt_max = pow(MASSE_DEUT,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
  // }
  // else if(A_nuclei == 3){
  //   mt_min = MASS_HELIUM_TROIS*1e-3;
  //   mt_max = pow(MASS_HELIUM_TROIS*1e-3,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
  // }

	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt);
  if(A_nuclei == 1){
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
  }
  else if(A_nuclei == 2){
    pT = sqrt(mt*mt - pow(MASSE_DEUT,2.0));
  }
  else if(A_nuclei == 3){
    pT = sqrt(mt*mt - pow(MASS_HELIUM_TROIS*1e-3,2.0));
  }

  f3     = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal) * mt / pL;

  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt);
    // if(A_nuclei == 1){
      pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    // }
    // else if(A_nuclei == 2){
    //   pT = sqrt(mt*mt - pow(MASSE_DEUT,2.0));
    // }
    // else if(A_nuclei == 3){
    //   pT = sqrt(mt*mt - pow(MASS_HELIUM_TROIS*1e-3,2.0));
    // }
    f2 = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal) * mt / pL;

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt);
    // if(A_nuclei == 1){
      pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    // }
    // else if(A_nuclei == 2){
    //   pT = sqrt(mt*mt - pow(MASSE_DEUT,2.0));
    // }
    // else if(A_nuclei == 3){
    //   pT = sqrt(mt*mt - pow(MASS_HELIUM_TROIS*1e-3,2.0));
    // }
    f3 = E_d3S_on_d3P_PBAR_H_ON_HE_LAB(E_proton,pL,pT,A_nuclei,P_coal) * mt / pL;

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de LAURENT DEROME ne concerne que la production d'antiprotons
* lors d'une collision proton sur helium (He). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte. De meme pour les antiTritium si l'on calcule le flux d'antiHe3.
*/
  if(A_nuclei == 1 || A_nuclei == 3)resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

/********************************************************************************************/
/*
* Ce module ecrit dans le fichier FILE_NAME_H_ON_HE le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs de la section efficace differentielle du processus
* proton (E_proton) + helium (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'energie TOTALE des protons E_proton varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_H_ON_HE_write_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_proton,i_pbar;
  double E_proton,T_pbar,E_pbar;
  FILE *p_to_pbar_file;
/*
* On remet a zero le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton].
*/
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton] = 0.0;
    }
  }
/*
* Et maintenant, on le remplit !
*/
  if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_H_ON_HE,"w");
  else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_H_ON_HE,"w");
  else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_H_ON_HE,"w");
  for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
  {
    E_proton = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      T_pbar = T_PBAR_MIN*A_nuclei *
      pow((T_PBAR_MAX/T_PBAR_MIN*A_nuclei),((double)i_pbar/(double)DIM_TAB_PBAR));
      E_pbar = T_pbar + MASSE_PROTON;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton] =
			dSpbar_sur_dEpbar_SIMPSON_H_ON_HE(E_proton,E_pbar,1000,A_nuclei,pt_Cross_Section->P_coal);
      fprintf(p_to_pbar_file," %.6e \n",pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_HE[i_pbar][i_proton]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* helium (E_nucleon) + proton (repos) -----> antiproton (E_pbar) + X.
* consideree dans le referentiel du laboratoire dans lequel le proton est au repos.
* Cette LI section efficace est exprimee en unites de [millibarns GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_HE_ON_H_LAB(double E_nucleon,double pL,double pT, double A_nuclei, double P_coal)
{
  double P_nucleon,E_pbar,mt_pbar,y_pbar,resultat,Ed3S_on_d3P_pbar,Ed3S_on_d3P_nbar,B;

  P_nucleon = sqrt(pow(E_nucleon,2) - pow(MASSE_PROTON,2));
	E_pbar    = sqrt(pow(MASSE_PROTON,2) + pT*pT + pL*pL);
	mt_pbar   = sqrt(pow(MASSE_PROTON,2) + pT*pT);
	y_pbar    = acosh(E_nucleon/MASSE_PROTON) - atanh(pL/E_pbar);
  Ed3S_on_d3P_pbar = GetInvarMul_pAapX(4.0,P_nucleon,y_pbar,mt_pbar) * GetReactionCrossSection(4.0,P_nucleon);

  if(A_nuclei == 1) resultat = Ed3S_on_d3P_pbar;
  if(A_nuclei == 2){
    B = MASSE_DEUT/(MASSE_PROTON*1e3*MASS_NEUTRON)*pow(4*PI/3*pow(P_coal,3)/8,2); //MeV^4
    B *= 1e-12; // conversion en GeV^4
    B /= pow(4,1./3)*sigma_total_pH(E_nucleon)/mb_cm2; //GeV^4 Simple scaling factor 2 to be updated
    resultat = B*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*1.3;
  }
  if(A_nuclei == 3){
    B = MASS_HELIUM_TROIS/(MASSE_PROTON*1e3*MASSE_PROTON*1e3*MASS_NEUTRON)*pow(4*PI/3*pow(P_coal,3)/8,2); //MeV^4
    B *= 1e-12; //GeV^4
    B /= pow(pow(4,1./3)*sigma_total_pH(E_nucleon)/mb_cm2,2); //GeV^4 Simple volume factor, to be updated
    resultat = B*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*Ed3S_on_d3P_pbar*1.3;
  }

  return resultat; /* [millibarns GeV^{-2}] */


}

/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* helium (E_nucleon) + proton (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du nucleon incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON_HE_ON_H(double E_nucleon,double E_pbar,long n_step, double A_nuclei, double P_coal)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_nucleon-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_nucleon+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON){return resultat;}

  if (E_pbar<=MASSE_PROTON){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_nucleon+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_nucleon-MASSE_PROTON)/(E_nucleon+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_nucleon = %.5e -- E_pbar = %.5e\n",E_nucleon,E_pbar);
		fprintf(probleme," more information ! E_nucleon = %.5e -- E_pbar = %.5e\n",E_nucleon,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
  if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}
//if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.00)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de LAURENT DEROME ne concerne que la production d'antiprotons
* lors d'une collision d'helium (He) sur proton. Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* d'helions de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE L'ANTIPROTON FINAL.
*/
  mt_min = MASSE_PROTON;
	mt_max = pow(MASSE_PROTON,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt);
  pT     = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
  f3     = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal) * mt / pL;

  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f2 = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal) * mt / pL;

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f3 = E_d3S_on_d3P_PBAR_HE_ON_H_LAB(E_nucleon,pL,pT,A_nuclei,P_coal) * mt / pL;

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de LAURENT DEROME ne concerne que la production d'antiprotons
* lors d'une collision d'helium (He) sur proton. Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* d'helions de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

/********************************************************************************************/
/*
* Ce module ecrit dans le fichier FILE_NAME_HE_ON_H le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs de la section efficace differentielle du processus
* helium (E_nucleon) + proton (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'energie TOTALE PAR NUCLEON E_nucleon de l'helium incident varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_HE_ON_H_write_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_nucleon,i_pbar;
  double E_nucleon,T_pbar,E_pbar;
  FILE *p_to_pbar_file;
/*
* On remet a zero le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon].
*/
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon] = 0.0;
    }
  }
/*
* Et maintenant, on le remplit !
*/
  if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_HE_ON_H,"w");
  else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_HE_ON_H,"w");
  else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_HE_ON_H,"w");
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    E_nucleon = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_nucleon/(double)DIM_TAB_PROTON));
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      T_pbar = T_PBAR_MIN*A_nuclei *
      pow((T_PBAR_MAX/T_PBAR_MIN*A_nuclei),((double)i_pbar/(double)DIM_TAB_PBAR));
      E_pbar = T_pbar + MASSE_PROTON;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon] =
			dSpbar_sur_dEpbar_SIMPSON_HE_ON_H(E_nucleon,E_pbar,1000,A_nuclei,pt_Cross_Section->P_coal);
      fprintf(p_to_pbar_file," %.6e \n",pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_H[i_pbar][i_nucleon]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module ecrit dans le fichier FILE_NAME_HE_ON_HE le tableau denote
* pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[DIM_TAB_PBAR+1][DIM_TAB_PROTON+1] dans lequel sont
* enregistrees les valeurs de la section efficace differentielle du processus
* helium (E_nucleon) + helium (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'energie TOTALE PAR NUCLEON E_nucleon de l'helium incident varie de E_PROTON_MIN a E_PROTON_MAX
* en prenant (DIM_TAB_PROTON + 1) valeurs differentes.
*
* L'energie CINETIQUE des antiprotons T_pbar varie de T_PBAR_MIN a T_PBAR_MAX
* en prenant (DIM_TAB_PBAR + 1) valeurs differentes.
*
*/
void DSPBAR_SUR_DEPBAR_HE_ON_HE_write_file(struct Structure_Cross_Section* pt_Cross_Section, double A_nuclei)
{
  long i_nucleon,i_pbar;
  double E_nucleon,T_pbar,E_pbar;
  FILE *p_to_pbar_file;
/*
* On remet a zero le tableau pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon].
*/
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon] = 0.0;
    }
  }
/*
* Et maintenant, on le remplit !
*/
  if(A_nuclei == 1)p_to_pbar_file = fopen(FILE_NAME_HE_ON_HE,"w");
  else if(A_nuclei == 2)p_to_pbar_file = fopen(FILE_NAME_DE_HE_ON_HE,"w");
  else if(A_nuclei == 3)p_to_pbar_file = fopen(FILE_NAME_HE3_HE_ON_HE,"w");
  for (i_nucleon=0;i_nucleon<=DIM_TAB_PROTON;i_nucleon++)
  {
    E_nucleon = E_PROTON_MIN *
    pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_nucleon/(double)DIM_TAB_PROTON));
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      T_pbar = T_PBAR_MIN*A_nuclei *
      pow((T_PBAR_MAX/T_PBAR_MIN*A_nuclei),((double)i_pbar/(double)DIM_TAB_PBAR));
      E_pbar = T_pbar + MASSE_PROTON;
      pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon] = pow((4.0),(1.0/3.0)) *
			dSpbar_sur_dEpbar_SIMPSON_H_ON_HE(E_nucleon,E_pbar,1000,A_nuclei,pt_Cross_Section->P_coal);
      fprintf(p_to_pbar_file," %.6e \n",pt_Cross_Section->DSPBAR_SUR_DEPBAR_HE_ON_HE[i_pbar][i_nucleon]);
    }
  }
  fclose(p_to_pbar_file);
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/* CODE DE REMI DUPERRAY et al. POUR LE CALCUL DE LA REACTION H ON HYDROGEN
*
* R.~P.~Duperray, C.-Y.~Huang, K.~V.~Protasov and M.~Bu\'enerd,
* ``Parameterization of the antiproton inclusive production cross section on nuclei,''
* Phys.\ Rev.\ D {\bf 68}, 094017 (2003)
* [arXiv:astro-ph/0305274].
*
*  Ce code retourne la multiplicite invariante 1/sig_reac * Edsig/d3p : p + H ---> ap + X
*  en unites de [GeV^{-2}] avec
*
*  p1 : projectile momentum (GeV/c) (target reference system)
*  y  : produced particle rapidity  (target reference system) y = atanh(pl/E) = atanh(beta_CM) + atanh(pl_star/E_star)
*  mt : produced particle transverse mass (GeV/c2) = sqrt(m^2+pt^2)
*
*  #define PMASS 0.93827231 en unites de [GeV]
*
*/
double invariant_multiplicity_pH_apX(double p1,double y,double mt)
{
  double Ecm,pt,xp,xt,sig,phi;
  double E,E1,p,ptmax;
  double s,Betacm;
  double m     = PMASS;
  double omega =  3.*m;
  double D1    = 3.4610;
	double D2    = 4.340;
	double D3    = 0.007855;
	double D4    = 0.5121;
	double D5    = 3.6620;
	double D6    = 0.023070;
	double D7    = 3.2540;
  double EcmMax;

  sig = 0;

  E1 = sqrt(p1*p1+PMASS*PMASS); /*E1 = energie totale par nucleon */

  s = 2*PMASS*PMASS+2*PMASS*E1;
  Betacm = p1/(E1+PMASS);


  EcmMax =  (s + m*m - omega*omega)/(2.0*sqrt(s));
  //    EcmMax = Emaxferm;
  Ecm = mt*cosh(y-atanh(Betacm));
  xp  = Ecm/EcmMax;

  if(mt<m) return 0;
  pt  = sqrt(mt*mt-m*m);
  E   = mt*cosh(y);
  p   = sqrt(E*E - m*m);

  ptmax = pow(EcmMax/cosh(y-atanh(Betacm)),2.) - m*m;
  if(ptmax<0) return 0;
  ptmax = sqrt(ptmax);
  xt = pt/ptmax;

  if(xp<1.0 && xt>=0 && xt<1 && xp > 0) {
    phi = D3 * pow(sqrt(s),D4) * exp(-D5*pt)  +  D6 * exp(-D7*pow(pt,2.));
		sig = phi * pow((1.0-xp),D1) * exp(-D2*xp);
  }
  else sig=0;
  return sig;
}

/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton (E_proton) + hydrogene (repos) -----> antiproton (E_pbar) + X.
* consideree dans le referentiel du laboratoire dans lequel l'hydrogene est au repos.
* Cette LI section efficace est exprimee en unites de [millibarns GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(double E_proton,double pL,double pT)
{
  double P_proton,E_pbar,mt_pbar,y_pbar,resultat;

  P_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
	E_pbar   = sqrt(pow(MASSE_PROTON,2) + pT*pT + pL*pL);
	mt_pbar  = sqrt(pow(MASSE_PROTON,2) + pT*pT);
	y_pbar   = atanh(pL/E_pbar);

//resultat = invariant_multiplicity_pH_apX(P_proton,y_pbar,mt_pbar) * GetReactionCrossSection(1.0,P_proton);
	resultat = invariant_multiplicity_pH_apX(P_proton,y_pbar,mt_pbar) * 1.0e27 * sigma_inelastic_pH_TAN_and_NG(E_proton);
	return resultat; /* [millibarns GeV^{-2}] */
}

/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + hydrogene (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_duperray(double E_proton,double E_pbar,long n_step)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_proton-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON){return resultat;}

  if (E_pbar<=MASSE_PROTON){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		fprintf(probleme," more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
	if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de REMI DUPERRAY ne concerne que la production d'antiprotons
* lors d'une collision proton sur hydrogene (H). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE L'ANTIPROTON FINAL.
*/
  mt_min = MASSE_PROTON;
	mt_max = pow(MASSE_PROTON,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt);
  pT     = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
  f3     = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT) * mt / pL;

  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f2 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT) * mt / pL;

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_duperray(E_proton,pL,pT) * mt / pL;

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de REMI DUPERRAY ne concerne que la production d'antiprotons
* lors d'une collision proton sur hydrogene (H). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + hydrogene (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_tan_ng_mass_T(double E_proton,double E_pbar,long n_step)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_proton-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON){return resultat;}

  if (E_pbar<=MASSE_PROTON){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		fprintf(probleme," more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
  if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de TAN & NG ne concerne que la production d'antiprotons
* lors d'une collision proton sur hydrogene (H). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE L'ANTIPROTON FINAL.
*/
  mt_min = MASSE_PROTON + 1.0e-6;
	mt_max = pow(MASSE_PROTON,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt);
  pT     = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
  f3     = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar) * mt / pL;

  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f2 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar) * mt / pL;

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f3 = E_d3S_on_d3P_PBAR_LAB(E_proton,pL,pT,&E_CMF,&pLstar,&pTstar) * mt / pL;

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * (1.e-27); /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de TAN & NG ne concerne que la production d'antiprotons
* lors d'une collision proton sur hydrogene (H). Dans notre contexte astrophysique,
* il y a autant d'antineutrons produits que d'antiprotons lors de la spallation
* de protons de haute energie sur du gaz interstellaire. Les antineutrons ne
* sont pas a priori detectes dans les experiences de haute energie au niveau
* des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un
* facteur 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
/*
* sigma_total_pH_MDGS est la section efficace totale de collision d'un proton cosmique
* avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
* Le fit est extrait de l'article arXiv:1408.0288 dont le titre est
* A new evaluation of the antiproton production cross section for cosmic ray studies.
*
*/
double sigma_total_pH_MDGS(double E_proton)
{
	double Zpp  = 33.44;
	double M    = 2.06;
	double Bpp  = PI * pow((2.99792458e10*6.58211915e-25/M),2.) / mb_cm2;
	double Y1pp = 13.53;
	double Y2pp = 6.38;
	double eta1 = 0.324;
	double eta2 = 0.324;

	double s    = 2.*MASSE_PROTON*(E_proton+MASSE_PROTON);
	double sM   = pow((2.*MASSE_PROTON + M),2.);

	double sig_tot_pp;
	sig_tot_pp = Zpp + Bpp*pow(log(s/sM),2.) + Y1pp*pow((sM/s),eta1) - Y2pp*pow((sM/s),eta2);
	return (sig_tot_pp*mb_cm2);
}

/********************************************************************************************/
/*
* sigma_elastic_pH_MDGS est la section efficace elastique de collision d'un proton cosmique
* avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
* Le fit est extrait de l'article arXiv:1408.0288 dont le titre est
* A new evaluation of the antiproton production cross section for cosmic ray studies.
*
*/
double sigma_elastic_pH_MDGS(double E_proton)
{
	double Zpp  = 144.98;
	double M    = 3.06;
	double Bpp  = PI * pow((2.99792458e10*6.58211915e-25/M),2.) / mb_cm2;
	double Y1pp = 2.64;
	double Y2pp = 137.27;
	double eta1 = 1.57;
	double eta2 = -(4.65e-3);

	double s    = 2.*MASSE_PROTON*(E_proton+MASSE_PROTON);
	double sM   = pow((2.*MASSE_PROTON + M),2.);

	double sig_el_pp;
	sig_el_pp = Zpp + Bpp*pow(log(s/sM),2.) + Y1pp*pow((sM/s),eta1) - Y2pp*pow((sM/s),eta2);
	return (sig_el_pp*mb_cm2);
}

/********************************************************************************************/
/*
* sigma_inelastic_pH_MDGS est la section efficace inelastique de collision d'un proton
* cosmique avec un atome d'hydrogene au repos. Cette section efficace est exprimee en [cm^{2}].
* Le fit est extrait de l'article arXiv:1408.0288 dont le titre est
* A new evaluation of the antiproton production cross section for cosmic ray studies.
*
*/
double sigma_inelastic_pH_MDGS(double E_proton)
{
	return (sigma_total_pH_MDGS(E_proton) - sigma_elastic_pH_MDGS(E_proton));
}

/********************************************************************************************/
/*
* Ces résultats sont extraits de l'article arXiv:1408.0288 dont le titre est
* A new evaluation of the antiproton production cross section for cosmic ray studies.
*
*  Ce code retourne la multiplicite invariante 1/sig_reac * Edsig/d3p : p + H ---> ap + X
*  en unites de [GeV^{-2}] avec
*
*  p1 : projectile momentum (GeV/c) (target reference system)
*  y  : produced particle rapidity  (target reference system) y = atanh(pl/E) = atanh(beta_CM) + atanh(pl_star/E_star)
*  mt : produced particle transverse mass (GeV/c2) = sqrt(m^2+pt^2)
*
*  #define PMASS 0.93827231 en unites de [GeV]
*
*/
double invariant_multiplicity_pH_apX_MDGS_F12(double p1,double y,double mt)
{
	double Ecm,pt,xp,xt,sig,phi;
	double E,E1,p,ptmax;
	double s,Betacm;
	double m     = PMASS;
	double omega =  3.*m;
	double EcmMax;

	double C1 = 4.499;
	double C2 = 3.41;
	double C3 = 0.00942;
	double C4 = 0.445;
	double C5 = 3.502;
	double C6 = 0.0622;
	double C7 = -0.247;
	double C8 = 2.576;

	sig    = 0;

	E1     = sqrt(p1*p1+PMASS*PMASS); /*E1 = energie totale par nucleon */
	s      = 2*PMASS*PMASS+2*PMASS*E1;
	Betacm = p1/(E1+PMASS);

	EcmMax =  (s + m*m - omega*omega)/(2.0*sqrt(s));
	Ecm    = mt*cosh(y-atanh(Betacm));
	xp     = Ecm/EcmMax;

	if(mt<m) return 0;
	pt  = sqrt(mt*mt-m*m);
	E   = mt*cosh(y);
	p   = sqrt(E*E - m*m);

	ptmax = pow(EcmMax/cosh(y-atanh(Betacm)),2.) - m*m;
	if(ptmax<0) return 0;
	ptmax = sqrt(ptmax);
	xt = pt/ptmax;

	if(xp<1.0 && xt>=0 && xt<1 && xp > 0)
	{
		phi = C3 * pow(sqrt(s),C4) * exp(-C5*pt)  +  C6 * pow(sqrt(s),C7) * exp(-C8*pow(pt,2.));
		sig = phi * pow((1.0-xp),C1) * exp(-C2*xp);
	}
	else sig=0;
	return sig;
}

/********************************************************************************************/
/*
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton (E_proton) + hydrogene (repos) -----> antiproton (E_pbar) + X.
* consideree dans le referentiel du laboratoire dans lequel l'hydrogene est au repos.
* Cette LI section efficace est exprimee en unites de [cm^{2} GeV^{-2}].
*
*/
double E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(double E_proton,double pL,double pT)
{
	double P_proton,E_pbar,mt_pbar,y_pbar,resultat;

	P_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
	E_pbar   = sqrt(pow(MASSE_PROTON,2) + pT*pT + pL*pL);
	mt_pbar  = sqrt(pow(MASSE_PROTON,2) + pT*pT);
	y_pbar   = atanh(pL/E_pbar);

	resultat = invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
	return resultat; /* [cm^{2} GeV^{-2}] */
}
/********************************************************************************************/
/* VP
* Ce module calcule la section efficace invariante de Lorentz (LI) de la reaction
* proton (E_proton) + hydrogene (repos) -----> antiHe3 (E_pbar) + X.
* consideree dans le referentiel du laboratoire dans lequel l'hydrogene est au repos.
* Cette LI section efficace est exprimee en unites de [cm^{2} GeV^{-2}].
* Nous utilisons le scenario de coalescence comme decrit dans 1711.08465.
*/
double E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(double E_proton,double pL,double pT, double P_coal, double A_nuclei, double Z_nuclei){

  double P_proton,P_proton_CMF,P_proton_reste,P_proton_reste_2,E_pbar,E_reste,E_reste_2,mt_pbar,y_pbar,resultat,Ed3S_on_d3P_pbar,Ed3S_on_d3P_nbar,B;
  double E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  double Ed3S_on_d3P_at_E_CMF=0,Ed3S_on_d3P_at_E_reste=0,Ed3S_on_d3P_at_E_reste_2=0;
  P_proton_CMF = sqrt(pow(E_CMF,2) - pow(MASSE_PROTON,2));
  // printf("P_proton_CMF %e \n",P_proton_CMF );
  P_proton = sqrt(pow(E_proton,2) - pow(MASSE_PROTON,2));
  // printf("P_proton%e \n",P_proton);


  E_pbar   = sqrt(pow(MASSE_PROTON,2) + pT*pT + pL*pL);
	mt_pbar  = sqrt(pow(MASSE_PROTON,2) + pT*pT);
	y_pbar   = atanh(pL/E_pbar);
  E_reste = E_CMF-2*E_pbar;
  E_reste_2 = E_CMF-4*E_pbar;
  if(E_reste<0)E_reste =0;
  if(E_reste_2<0)E_reste_2 =0;
  if(E_reste<MASSE_PROTON)P_proton_reste = 0;
  else P_proton_reste = sqrt(pow(E_reste,2) - pow(MASSE_PROTON,2));
  if(E_reste_2<MASSE_PROTON)P_proton_reste_2 = 0;
  else P_proton_reste_2 = sqrt(pow(E_reste_2,2) - pow(MASSE_PROTON,2));

  // if(isnan(P_proton_reste)==1)P_proton_reste=0;
  // printf("E_proton_reste %e P_proton_reste %e\n",E_reste,P_proton_reste);
  // printf("E_CMF %e P_proton_CMF %e\n",E_CMF,P_proton_CMF);
  // printf("E_proton %e P_proton %e\n",E_proton,P_proton);

  // B = MASS_HELIUM_TROIS/(MASSE_PROTON*1e3*MASSE_PROTON*1e3*MASS_NEUTRON)*pow(4*PI/3*pow(P_coal,3)/8,2); //MeV^4
  B = M_NUCLEI/(pow(MASSE_PROTON,Z_nuclei)*pow(MASS_NEUTRON*1e-3,A_nuclei-Z_nuclei))*pow(4*PI/3*pow(P_coal,3)/8,A_nuclei-1); //MeV^4
  // if(A_nuclei == 2)B *= 1e-6; //GeV^4
  // if(A_nuclei == 3)B *= 1e-12; //GeV^4
  // printf("B %e\n", B);
  B /= pow(sigma_total_pH_MDGS(E_proton),Z_nuclei); //GeV^4
  if(A_nuclei == 1)resultat = invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
  else if(A_nuclei == 2){

    if(A_LA_CHARDONNET == _FALSE_){
      Ed3S_on_d3P_pbar = invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
      Ed3S_on_d3P_nbar = 1.3*invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
      resultat = B*pow(Ed3S_on_d3P_pbar,Z_nuclei)*pow(Ed3S_on_d3P_nbar,A_nuclei-Z_nuclei); //TBC
    }

    /*** Computation as advocated by Chardonnet et al; astro-ph/9005110 */
    else if(A_LA_CHARDONNET == _TRUE_){
      Ed3S_on_d3P_at_E_CMF = invariant_multiplicity_pH_apX_MDGS_F12(P_proton_CMF,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_CMF);
      Ed3S_on_d3P_at_E_reste = invariant_multiplicity_pH_apX_MDGS_F12(P_proton_reste,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_reste);
      resultat = 0.5*B*(
                  Ed3S_on_d3P_at_E_CMF
                * 1.3*Ed3S_on_d3P_at_E_reste
                + Ed3S_on_d3P_at_E_reste
                * 1.3*Ed3S_on_d3P_at_E_CMF);
    }
  }
  else if(A_nuclei == 3){
    if(A_LA_CHARDONNET == _FALSE_){
      Ed3S_on_d3P_pbar = invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
    	Ed3S_on_d3P_nbar = 1.3*invariant_multiplicity_pH_apX_MDGS_F12(P_proton,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_proton);
      resultat = B*pow(Ed3S_on_d3P_pbar,Z_nuclei)*pow(Ed3S_on_d3P_nbar,A_nuclei-Z_nuclei); //TBC
    }
    else if(A_LA_CHARDONNET == _TRUE_){
      /*** Computation as advocated by Chardonnet et al; astro-ph/9005110 */
      Ed3S_on_d3P_at_E_CMF = invariant_multiplicity_pH_apX_MDGS_F12(P_proton_CMF,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_CMF);
      Ed3S_on_d3P_at_E_reste = invariant_multiplicity_pH_apX_MDGS_F12(P_proton_reste,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_reste);
      Ed3S_on_d3P_at_E_reste_2 = invariant_multiplicity_pH_apX_MDGS_F12(P_proton_reste_2,y_pbar,mt_pbar) * sigma_inelastic_pH_MDGS(E_reste_2);
      resultat = B*Ed3S_on_d3P_at_E_CMF*Ed3S_on_d3P_at_E_reste*Ed3S_on_d3P_at_E_reste_2;
    }
  }

  // if(resultat!=0)printf("resultat %e  B %e Ed3S_on_d3P_at_E_CMF %e Ed3S_on_d3P_at_E_reste %e Ed3S_on_d3P_at_E_reste_2 %e\n", resultat, B,Ed3S_on_d3P_at_E_CMF,Ed3S_on_d3P_at_E_reste,Ed3S_on_d3P_at_E_reste_2);
  return resultat; /* [cm^{2} GeV^{-2}] */
}

/********************************************************************************************/
/*
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + hydrogene (repos) -----> antiproton (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dSpbar_sur_dEpbar_SIMPSON_H_ON_H_MDGS_F12(double E_proton,double E_pbar,long n_step)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_proton-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  if (E_CMF<=4.0*MASSE_PROTON){return resultat;}

  if (E_pbar<=MASSE_PROTON){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		fprintf(probleme," more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
	if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar;
  pT = sin_theta * P_pbar;
  f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f2 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar;
    pT = sin_theta * P_pbar;
    f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar; /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de MDGS ne concerne que la production d'antiprotons lors d'une collision proton sur hydrogene (H).
* Dans notre contexte astrophysique, il y a autant, voire plus, d'antineutrons produits que d'antiprotons lors de la
* spallation de protons de haute energie sur du gaz interstellaire. Les antineutrons ne sont pas a priori detectes dans
* les experiences de haute energie au niveau des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un facteur >= 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE L'ANTIPROTON FINAL.
*/
  mt_min = MASSE_PROTON + 1.0e-6;
	mt_max = pow(MASSE_PROTON,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt);
  pT     = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
  f3     = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT) * mt / pL;

  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f2 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT) * mt / pL;

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt);
    pT = sqrt(mt*mt - pow(MASSE_PROTON,2.0));
    f3 = E_d3S_on_d3P_PBAR_H_ON_H_LAB_MDGS_F12(E_proton,pL,pT) * mt / pL;

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI; /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de MDGS ne concerne que la production d'antiprotons lors d'une collision proton sur hydrogene (H).
* Dans notre contexte astrophysique, il y a autant, voire plus, d'antineutrons produits que d'antiprotons lors de la
* spallation de protons de haute energie sur du gaz interstellaire. Les antineutrons ne sont pas a priori detectes dans
* les experiences de haute energie au niveau des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un facteur >= 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/

/********************************************************************************************/
/* VP
* Ce module evalue la section efficace differentielle du processus
* proton (E_proton) + hydrogene (repos) -----> antiHe3 (E_pbar) + X.
* La quantite \frac{d \sigma_{\pbar}}{d E_{\pbar}} est exprimee en [cm^{2} GeV^{-1}].
*
* L'integrale sur le cosinus de l'angle \theta que font les impulsions
* de l'antiproton final et du proton incident est menee par la methode de SIMPSON.
*
* Lorsque cet angle est petit, une integrale sur la masse transverse mt de l'antiproton final
* est alors menee par la methode de SIMPSON.
*
*/
double dS3He_sur_sE3He_H_on_H_MDGS_F12_coalescence(double E_proton,double E_pbar, double P_coal, long n_step, double A_nuclei, double Z_nuclei, double M_nuclei)
{
  extern FILE *probleme;
	double P_pbar,E_CMF,E_CMF_SEUIL,E_MAX_star,gamma,beta,pL,pT,pLstar,pTstar;
  long   i_theta,i_mt;
  double cos_theta,sin_theta,dcos_theta,cos_theta_min,cos_theta_max;
	double mt,dmt,mt_min,mt_max;
  double h,x1,x2,x3,f1,f2,f3;

  double resultat = 0.0;

  if (E_pbar>=(E_proton-2.*MASSE_PROTON)){return resultat;}

  E_CMF = sqrt(2.*MASSE_PROTON*(E_proton+MASSE_PROTON));
  E_CMF_SEUIL = sqrt(2.*MASSE_PROTON*(E_PROTON_MIN+MASSE_PROTON));
  // if (E_CMF<=4.0*MASSE_PROTON){return resultat;}
  // printf("E_proton %e E_CMF %e GeV seuil %e \n",E_proton, E_CMF,8*MASSE_PROTON);
  if (E_CMF<=E_CMF_SEUIL){return resultat;}//the threshold for De and He3 production is higher

  if (E_pbar<=M_nuclei){return resultat;}
  P_pbar = sqrt(pow(E_pbar,2) - pow(M_nuclei,2));

  gamma = sqrt((E_proton+MASSE_PROTON)/(2.*MASSE_PROTON));
  beta  = sqrt((E_proton-MASSE_PROTON)/(E_proton+MASSE_PROTON));
  E_MAX_star = (pow(E_CMF,2) - pow(3.*MASSE_PROTON,2) + pow(MASSE_PROTON,2)) / (2.*E_CMF);


  cos_theta_max = 1. - 1.e-12; /* La valeur 1 entraine des problemes avec le sin_theta */
  cos_theta_min = (E_pbar - (E_MAX_star/gamma)) / beta / P_pbar;
	if (cos_theta_min>=cos_theta_max){return resultat;}
	if (cos_theta_min<=0.0)
	{
	  printf(" cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		fprintf(probleme," cos_theta_min = %.5e and NEGATIVE !\n",cos_theta_min);
		printf(" more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		fprintf(probleme," more information ! E_proton = %.5e -- E_pbar = %.5e\n",E_proton,E_pbar);
		return resultat;
	}
//if (cos_theta_min<=-1.0) {cos_theta_min = -1.0;}
	if (cos_theta_min<cos_theta_max && cos_theta_min>=(0.99)){goto TRANSVERSE_MASS;}

/*
* INTEGRATION SUR LE COSINUS DE L'ANGLE THETA.
*/
  dcos_theta = (cos_theta_max - cos_theta_min) / (double) n_step;
  h = dcos_theta/2.;
  x1 = cos_theta_min - dcos_theta;
  x2 = x1 + h;
  x3 = x2 + h;

  cos_theta = x3;
  sin_theta = sqrt(1. - pow(cos_theta,2));
  pL = cos_theta * P_pbar / A_nuclei; //We evaluate the individual nucleus production CS at k_A/A
  pT = sin_theta * P_pbar / A_nuclei;
  f3 = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei);

  for (i_theta=1;i_theta<=n_step;i_theta++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

    cos_theta = x2;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar / A_nuclei;
    pT = sin_theta * P_pbar / A_nuclei;
    f2 = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei);

    cos_theta = x3;
    sin_theta = sqrt(1. - pow(cos_theta,2));
    pL = cos_theta * P_pbar / A_nuclei;
    pT = sin_theta * P_pbar / A_nuclei;
    f3 = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei);

    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI * P_pbar; /* [cm^{2} GeV^{-1}] */
/*
* La section efficace de MDGS ne concerne que la production d'antiprotons lors d'une collision proton sur hydrogene (H).
* Dans notre contexte astrophysique, il y a autant, voire plus, d'antineutrons produits que d'antiprotons lors de la
* spallation de protons de haute energie sur du gaz interstellaire. Les antineutrons ne sont pas a priori detectes dans
* les experiences de haute energie au niveau des accelerateurs de particules alors qu'ils se desintegrent au bout de
* \gamma \times 10 minutes dans la galaxie. Il nous faut donc rajouter un facteur >= 2 pour les prendre en compte.
*/
  resultat *= 2.; /* [cm^{2} GeV^{-1}] */
  // printf("resultat integration %e\n", resultat);
	return resultat;


TRANSVERSE_MASS :
/*
* INTEGRATION SUR LA MASSE TRANSVERSE M_T DE LA PARTICULE (pbar, Debar, Hebar) FINALE.
*/
  mt_min = M_nuclei + 1.0e-6;
	mt_max = pow( M_nuclei,2.0) + pow(P_pbar,2.0)*(1.0 - cos_theta_min*cos_theta_min);
	mt_max = sqrt(mt_max);
	dmt    = (mt_max - mt_min) / (double) n_step;
  h      = dmt/2.;
  x1     = mt_min - dmt;
  x2     = x1 + h;
  x3     = x2 + h;

  mt     = x3;
  pL     = sqrt(E_pbar*E_pbar - mt*mt)/A_nuclei;
  pT     = sqrt(mt*mt - pow(M_nuclei,2.0))/A_nuclei;
  f3     = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei) * mt / (pL*A_nuclei);
  for (i_mt=1;i_mt<=n_step;i_mt++)
  {
    x1 = x3;
    x2 = x1 + h;
    x3 = x2 + h;

    f1 = f3;

		mt = x2;
    pL = sqrt(E_pbar*E_pbar - mt*mt)/A_nuclei;
    pT = sqrt(mt*mt - pow(M_nuclei,2.0))/A_nuclei;
    f2 = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei) * mt /(pL*A_nuclei);

    mt = x3;
		pL = sqrt(E_pbar*E_pbar - mt*mt)/A_nuclei;
    pT = sqrt(mt*mt - pow(M_nuclei,2.0))/A_nuclei;
    f3 = E_d3S_on_d3P_3HE_H_ON_H_LAB_MDGS_F12_coalescence(E_proton,pL,pT,P_coal,A_nuclei,Z_nuclei) * mt /(pL*A_nuclei);

    // resultat +=0;
    resultat += h * (f1/3. + 4.*f2/3. + f3/3.);
  }
  resultat *= 2. * PI; /* [cm^{2} GeV^{-1}] */
/*
* Nous prenons maintenant en compte la production d'antineutron qui se desintegrera en antiproton en multipliant simplement par un facteur 2.3 (tient compte de la brisure de symmetrie d'isospin)
*/
  if(A_nuclei == 1)resultat *= 2.3; /* [cm^{2} GeV^{-1}] */
	return resultat;
}

// /********************************************************************************************/
// /********************************************************************************************/
// /********************************************************************************************/
