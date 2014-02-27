#include "BESSEL_PRELIM.h"

/********************************************************************************************/
/********************************************************************************************/
/*
* Ce module cherche les NDIM premiers zéros alpha_i[i] de la fonction de Bessel J0.
* Puis le calcul des NDIM coefficients q_i[i] est effectué. Les résultats sont stockés
* dans le fichier BESSEL_ALPHA_QI
*/
void bessel_preliminary_write_file(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium)
{
 
  long i;
  FILE *bessel_file;
  
  bessel_file = fopen("BESSEL_ALPHAi_Qi","w");
  search_for_zero(alpha_i);
  production(alpha_i, pt_Proton);
  production(alpha_i, pt_Helium);
  for(i=1;i<=NDIM;i++)
  {
	//printf(" i = %4d , alpha_i = %.12e , q_i_proton = %.12e , q_i_helium = %.12e \n",i,alpha_i[i],pt_Proton->q_i[i],pt_Helium->q_i[i]);
    fprintf(bessel_file," %4d\t %.12e\t %.12e\t %.12e\t \n",i,alpha_i[i],pt_Proton->q_i[i],pt_Helium->q_i[i]);
  }
  fclose(bessel_file);
  return;
}

/********************************************************************************************/
/*
* Ce module lit les résultats du fichier BESSEL_ALPHA_QI calculés auparavent.
* Il charge donc les NDIM premiers zéros alpha_i[i] de la fonction de Bessel J0
* ainsi que les NDIM coefficients q_i[i]. Cela permet de gagner du temps de calcul.
*/
void bessel_preliminary_read_file(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium)
{
 
  long i,j;
  FILE *bessel_file;
  
  bessel_file = fopen("BESSEL_ALPHAi_Qi","r");
  for(i=1;i<=NDIM;i++)
  {
	fscanf(bessel_file," %4ld %lf %lf %lf",&j,&alpha_i[i],&pt_Proton->q_i[i],&pt_Helium->q_i[i]);
    //printf(" i = %4d , alpha_i = %.12e , q_i_proton = %.12e , q_i_helium = %.12e \n",j,alpha_i[i],pt_Proton->q_i[i],pt_Helium->q_i[i]);
  }
  fclose(bessel_file);
  return;
}

/********************************************************************************************/
/********************************************************************************************/
void search_for_zero(double alpha_i[NDIM+1])
{
  int i,nzero;
  long iz;
  double z,dz,zero;
/*
* On remet à zero le tableau alpha_i[i].
*/
  for (i=0;i<=NDIM;i++)
  {
    alpha_i[i] = 0.0;
  }
  
  nzero = 1;
  z = 2.0;
  dz = 0.01;
  for (iz=1;iz<=200000 && nzero<=NDIM;iz++)
  {
    z += dz;
    if ( (besselj0(z)*besselj0(z+dz)) <= 0 )
    {
/*    Le zéro de la fonction de Bessel J0 se situe entre z et z+dz    */
      search(z,z+dz,&zero);
      alpha_i[nzero] = zero;
      nzero += 1;
    }
  }
  return;
}

/********************************************************************************************/
void search(double z1,double z2,double *zero)
{
  extern FILE *probleme;
  short i_recur;
  double zinf,zsup,zmed;
  
  if (z1>z2)
  {
    fprintf(probleme,
    " ERREUR DANS LES ARGUMENTS DE LA FONCTION search \n"
    " z1 = %f , z2 = %f \n",z1,z2);
    return;
  }
  
  zinf = z1;
  zsup = z2;
  for (i_recur=1;i_recur<=20;i_recur++)
//for (i_recur=1;i_recur<=50;i_recur++)
  {
    zmed = (zinf+zsup)/2.;
    if ( (besselj0(zinf)*besselj0(zmed)) <= 0 )
      zsup = zmed;
    if ( (besselj0(zsup)*besselj0(zmed)) <= 0 )
      zinf = zmed;
  }
  *zero = zmed;
  return;
}

/********************************************************************************************/
void production(double alpha_i[NDIM+1], struct Structure_Nuclei* pt_Nuclei)
{

  short i;
  long i_int;
  double u,du,sum,coefficient,weight_SIMPSON;
/*
* On remet à zéro le tableau q_i[i].
*/
  for (i=0;i<=NDIM;i++)
  {
    pt_Nuclei->q_i[i] = 0.0;
  }
  
  sum = 0.0;
  du  = 1. / (double) NINT_PRODUCTION;
  u   = 0.0;
  for (i_int=0;i_int<=NINT_PRODUCTION;i_int++)
  {
	  if (i_int==0 || i_int==NINT_PRODUCTION) {weight_SIMPSON = 1./3.;}
		else {weight_SIMPSON = (1. + (double)(i_int % 2)) * 2. / 3.;}

    sum += weight_SIMPSON*f_PSRD(u)*u*du;
    u += du;
  }
  coefficient = (1.0)/(PI*pow(CM_PAR_KPC*R_GAL,2))/sum; /* [cm^{-2}] */
  
/*
* On calcule maintenant les intégrales q_i de la transformée de
* Bessel de la fonction f_PSRD.
*/
  for (i=1;i<=NDIM;i++)
  {
    sum = 0.0;
    du  = 1. / (double) NINT_PRODUCTION;
    u   = 0.0;
    for (i_int=0;i_int<=NINT_PRODUCTION;i_int++)
    {
		  if (i_int==0 || i_int==NINT_PRODUCTION) {weight_SIMPSON = 1./3.;}
			else {weight_SIMPSON = (1. + (double)(i_int % 2)) * 2. / 3.;}

      sum += weight_SIMPSON*f_PSRD(u)*besselj0(alpha_i[i]*u)*u*du;
      u   += du;
    }
    pt_Nuclei->q_i[i] = coefficient*sum/pow(besselj1(alpha_i[i]),2); /* [cm^{-2}] */
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
* fonction PRIMARY SOURCE RADIAL DISTRIBUTION.
*
* Cette fonction donne le profil des sources de rayons cosmiques primaires dans
* le plan de notre galaxie. Ce profil radial suit la distribution des pulsars.
* return (pow(u,1.2)*exp(- 6.44*u)); VIEUX PROFIL
* return (pow(u,0.6)*exp(- 3.00*u)); NOUVEAU PROFIL
*
*/
double f_PSRD(double u)
{
  double r,r0,rs,fr1,fr2,z0,a,b,resultat;

  r = u * R_GAL;
/*
	resultat = pow((r/8.5),2.0) * exp((-3.53)*(r - (8.5))/(8.5));// :p_salati:20060413 DAVID_RICHARD
*/
/*
	r0       = 1.528;
	a        = 2.35;
	resultat = pow(r,a) * exp(-r/r0);// L04
*/
/*
	r0       = 2.0;
	resultat = exp(-r/r0);// MS10
*/
/*
  DISTRIBUTION OF YUSIFOV AND KUCUK. 
*/
/*
	r0       = 1.25;//[kpc]
	a        = 4.0;
	z0       = E_DISC;//[kpc]
	resultat = pow(r,a) * exp(-r/r0);
*/
/*
  DISTRIBUTION OF YUSIFOV AND KUCUK MODIFIEE PAR GUILEHM ET TIMUR.
*/
/*
	r0       = 0.55;//[kpc]
	rs       = 8.50;//[kpc]
	a        = 1.64;
	b        = 4.01;
	z0       = E_DISC;//[kpc]
	fr1      = (r+r0)/(rs+r0);
	fr2      = (r-rs)/(rs+r0);
	resultat = pow(fr1,a) * exp((-1.0)*b*fr2);
*/

  resultat = pow(r,0.6) * exp(- 3.0*r/20.0);
	return resultat;
}

/********************************************************************************************/
/********************************************************************************************/
