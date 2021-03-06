#include "SOLAR_MOD.h"

/********************************************************************************************/
/********************************************************************************************/
/* Cette fonction a besoin en entree du numero atomique A et de la charge Z du
noyau considere. Ensuite, il convient de lui fournir la valeur de la modulation
solaire PHI au moment des faits. ATTENTION, PHI EST EXPRIME EN [GV] !
L'ENERGIE TOTALE PAR NUCLEON DANS LE MILIEU INTERSTELLAIRE EST NOTEE EnIS.
Le flux interstellaire est note flux_IS. Cette fonction calcule alors
l'ENERGIE TOTALE PAR NUCLEON EnTOA au niveau de la Terre (Top Of the Atmosphere)
ainsi que le flux correspondant flux_TOA.
*/
void FFA_IS_to_TOA(double A,double Z,double PHI,
double EnIS,double flux_IS,double *EnTOA,double *flux_TOA)
{
  double pnIS,pnTOA;
  double pnC,EnC,En_min,En_trans;
  
  pnC = Z * RIGIDITY_MS_C / A; /* Dans nos unites, la charge de l'electron vaut 1. */
  EnC = sqrt(pow(MASSE_PROTON,2) + pow(pnC,2));
  
  En_trans = EnC + Z*PHI/A;
  En_min   = En_trans + pnC*log(MASSE_PROTON/(EnC+pnC));
  
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
  else if (EnIS <= En_trans)
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
/* Cette fonction a besoin en entree du numero atomique A et de la charge Z du
noyau considere. Ensuite, il convient de lui fournir la valeur de la modulation
solaire PHI au moment des faits. ATTENTION, PHI EST EXPRIME EN [GV] !
L'ENERGIE TOTALE PAR NUCLEON AU SOMMET DE L'ATMOSPHERE EST NOTE EnTOA.
Le flux terrestre est note flux_TOA. Cette fonction calcule alors
l'ENERGIE TOTALE PAR NUCLEON EnIS dans le milieu interstellaire ainsi que
le flux correspondant flux_IS.
*/
void FFA_TOA_to_IS(double A,double Z,double PHI,
double EnTOA,double flux_TOA,double *EnIS,double *flux_IS)
{
  double pnIS,pnTOA;
  double pnC,EnC;
  
  pnC = Z * RIGIDITY_MS_C / A; /* Dans nos unites, la charge de l'electron vaut 1. */
  EnC = sqrt(pow(MASSE_PROTON,2) + pow(pnC,2));
  
  if (EnTOA <= MASSE_PROTON)
  {
    printf(
    " TON NOYAU NE PEUT PAS EXISTER AINSI !\n"
    " SON ENERGIE PAR NUCLEON DOIT ETRE SUPERIEURE A %.5e [GEV]\n",
    MASSE_PROTON);
    return;
  }
  else if (EnTOA <= EnC)
  {
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
    *EnIS = EnC + pnC*log((EnTOA+pnTOA)/(EnC+pnC)) + Z*PHI/A;
    pnIS  = sqrt(pow(*EnIS,2) - pow(MASSE_PROTON,2));
    *flux_IS = pow((pnIS/pnTOA),2) * flux_TOA;
    return;
  }
  else
  {
    *EnIS = EnTOA + Z*PHI/A;
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
    pnIS  = sqrt(pow(*EnIS,2) - pow(MASSE_PROTON,2));
    *flux_IS = pow((pnIS/pnTOA),2) * flux_TOA;
    return;
  }
}
/********************************************************************************************/
/********************************************************************************************/
