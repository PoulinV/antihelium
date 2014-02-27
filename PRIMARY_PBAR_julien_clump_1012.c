#include "PRIMARY_PBAR_julien_clump_1012.h"

/********************************************************************************************/
/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Pbar->BESSEL_PBAR_PRI_Epbar_i
* en fonction de l'énergie CINETIQUE des antiprotons T_pbar et des coefficients de BESSEL i.
*
*/
void calculation_BESSEL_PBAR_PRIMARY_Epbar_i(long n_vert, long n_rad, double alpha_i[NDIM+1],
                                             struct Structure_Pbar* pt_Pbar,
                                             struct Structure_Propagation* pt_Propagation,
                                             struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
  long i_pbar,i;
  double T_pbar,E_pbar,impulsion_pbar,v_pbar,K_pbar;
  double Si,Abar_i;

  long i_vert,i_rad;
  double x_vert,z_vert,dx_vert,weight_SIMPSON_vert;
  double x_rad, r_rad, dx_rad, weight_SIMPSON_rad;
  static double q_pbar_primary_i_z[2001];
/*
* Protection concernant le fait que n_vert doit être inférieur ou égal à 1000
* mais pas supérieur !
*/
  if (n_vert > 1000)
  {
    printf(" ATTENTION ! n_vert > 1000 !");
    exit (0);
  }
/*
* On remet à zéro le tableau pt_Pbar->BESSEL_PBAR_PRI_Epbar_i.
* On pourrait d'ailleurs le remettre à zero plus loin.
*/
  for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
  {
    for (i=0;i<=NDIM;i++)
    {
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
    }
  }
/*
* On remplit maintenant le tableau pt_Pbar->BESSEL_PBAR_PRI_Epbar_i.
*/
  for (i=1;i<=NDIM;i++)
  {
/*
*   Nous commençons par le calcul de l'intégrale radiale q_pbar_primary_i_z.
*   \beq
*   q_pbar_primary_i_z = \frac{2}{J_{1}^{2}(\alpha_{i})} \times
*   \int_{0}^{1} \, x \, dx \, J_{0} \left( \alpha_{i} x \right) \,
*   \left\{ \frac{\rho_{\chi}}{\rho_0} \right\}^{2} \;\; .
*   \eeq
*
*   Cette intégrale dépend de i et de z et nous stockons les résultats
*   correspondant à l'ordre de Bessel i dans le tableau q_pbar_primary_i_z
*   qui ne dépend que de z.
*
*   Nous entamons une boucle sur z suivie par une boucle sur r.
*/
    dx_vert = 1. / (double) (2*n_vert);
    x_vert  = 0.;
    for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
    {
      z_vert = pt_Propagation->E_DIFFUS * x_vert;
/*    z_vert est la coordonnée verticale exprimée en [kpc].
*/
      q_pbar_primary_i_z[i_vert] = 0.0;
      dx_rad = 1. / (double) (2*n_rad);
      x_rad  = 0.0;
      for (i_rad=0;i_rad<=(2*n_rad);i_rad++)
      {
        r_rad = R_GAL * x_rad;
/*      r_rad est le rayon galactocentrique exprimé en [kpc].
*/
        if (i_rad==0 || i_rad==(2*n_rad)) {weight_SIMPSON_rad = 1./3.;}
        else {weight_SIMPSON_rad = (1. + (double)(i_rad % 2)) * 2. / 3.;}
	
		q_pbar_primary_i_z[i_vert] += (x_rad * dx_rad * weight_SIMPSON_rad) * besselj0(alpha_i[i]*x_rad) * pow(rapport_rho_chi_sur_rho_0_Einasto(r_rad,z_vert),2.0); //[NO UNIT].
			
        x_rad += dx_rad;
      }
      q_pbar_primary_i_z[i_vert] *= 2. / pow(besselj1(alpha_i[i]),2.0);
/*    [NO UNIT].
*/
      x_vert += dx_vert;
    }
/*
*   Nous continuons ensuite par une boucle sur l'énergie des antiprotons
*   englobant la boucle sur la variable verticale permettant ainsi de
*   calculer pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i].
*/
    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      T_pbar = T_PBAR_MIN *
      pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
      E_pbar = T_pbar + MASSE_PROTON;
      impulsion_pbar = sqrt(pow(E_pbar,2) - pow(MASSE_PROTON,2));
      v_pbar         = CELERITY_LIGHT * impulsion_pbar / E_pbar;
      K_pbar         = K_space_diffusion(E_pbar,MASSE_PROTON,1.0,pt_Propagation);
	  
/*    Si est exprimé en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_pbar,2));


/*    Abar_i est exprimé en [cm s^{-1}].  */
      Abar_i  = pt_Propagation->VENT_GALACTIQUE;
	  Abar_i += 2.0*E_DISC*CM_PAR_KPC *
      ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) * v_pbar * (DENSITE_H_DISC + pow(4.,(2./3.))*1.0*DENSITE_HE_DISC));
      Abar_i += K_pbar * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);

	  pt_Pbar->TABLE_Abar_i[i_pbar][i] = Abar_i;


/*    Il convient maintenant d'intégrer sur la variable verticale x_vert
*     ainsi que sur la variable radiale x_rad.
*/
      dx_vert = 1. / (double) (2*n_vert);
      x_vert  = 0.;
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
      for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
      {
        z_vert = pt_Propagation->E_DIFFUS * x_vert;
/*      z_vert est la coordonnée verticale exprimée en [kpc].
*/
        if (i_vert==0 || i_vert==(2*n_vert)) {weight_SIMPSON_vert = 1./3.;}
        else {weight_SIMPSON_vert = (1. + (double)(i_vert % 2)) * 2. / 3.;}
      

//	Calcul avec les spectres de Gaelle :
		
        pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] += (dx_vert * weight_SIMPSON_vert) *
        q_pbar_primary_i_z[i_vert] *
        pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *
        exp(- pt_Propagation->VENT_GALACTIQUE*z_vert*CM_PAR_KPC / (2.*K_pbar));
        sinh((Si/2.)*(pt_Propagation->E_DIFFUS-z_vert)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS); //[antiprotons cm^{-3} s^{-1} GeV^{-1}].
						
		
		x_vert += dx_vert;
      }
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] *= 2. * pt_Propagation->E_DIFFUS*CM_PAR_KPC / Abar_i;
	  
	  
/*    S'exprime maintenant en unités de [antiprotons cm^{-3} GeV^{-1}].
*/
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
 * Nous calculons ici le rapport $\left\{ \frac{\rho_{\chi}}{\rho_0} \right\}$
 * où $\rho_{\chi}$ désigne la masse volumique des neutralinos exprimée
 * en [GeV cm^{-3}] et dépendante des coordonnées cylindriques rr et z.
 * La masse volumique de référence est $\rho_0$ = 1 [GeV cm^{-3}].
 * Le résultat final est sans dimension et dépend de DEUX variables :
 * - Le rayon galactocentrique rr est exprimé en [kpc].
 * - La distance verticale z est également exprimée en [kpc].
 *
 */
double rapport_rho_chi_sur_rho_0(double rr,double z)
{
	double alpha,beta,gamma,core;
	double r,x,rho_chi;
	double enhancement = 1.0;
/*
	double RC_SMBH      = 0.5;
	double RC_SMBH      = 0.001;
	double RC_SMBH      = 0.1; VALEUR UTILISEE PAR PIERRE BRUN DANS SON MONTE_CARLO
*/
	double RC_SMBH      = 0.5;
	double ETA_CANONIC  = 1.0;
	double ETA_KRAVTSOV = (3.0 / (3.0 - (2.0*0.4)));
	double ETA_NFW      = 3.0;
	double ETA_MOORE    = (3.0 / (3.0 - (2.0*1.3)));
/*
	CANONICAL ISOTHERMAL CORE PROFILE
	alpha = 2.0;
	beta  = 2.0;
	gamma = 0.0;
	//core  = 4.0;
	//core  = 0.5;
	//core  = 2.8;
	core  = 5.0;

	PROFILE DE KRAVTSOV
	alpha = 2.0;
	beta  = 3.0;
	gamma = 0.4;
	core  = 10.0;

	PROFILE DE NFW
	alpha = 1.0;
	beta  = 3.0;
	gamma = 1.0;
	core  = 25.0;

	PROFILE DE NFW UTILISE PAR PIERRE BRUN DANS SON MONTE_CARLO
	alpha = 1.0;
	beta  = 3.0;
	gamma = 1.0;
	core  = 20.0;

	PROFILE DE MOORE
	alpha = 1.5;
	beta  = 3.0;
	gamma = 1.5;
	core  = 30.0;

	PROFILE DE MOORE NOUVEAU
	alpha = 1.5;
	beta  = 3.0;
	gamma = 1.3;
	core  = 30.0;
*/
	alpha = 1.0;
	beta  = 3.0;
	gamma = 1.0;
  //core  = 25.0; /* [kpc] */
	core  = 20.0; /* [kpc] */

	double ETA_TIMUR_prime = 8.0 * gamma * (pow(PI,(2.0)) - 9.0 + 6.0*gamma) / 9.0 / (3.0 - 2.0*gamma);
	double ETA_TIMUR       = ETA_TIMUR_prime + (2.0*gamma);

	r = sqrt(rr*rr + z*z);
	if (r <= RC_SMBH)
	{
		enhancement   = ETA_TIMUR + ETA_TIMUR_prime;
		x = r / RC_SMBH;
		if (1.e-3 <= x)
		{
			//enhancement = sin(PI*x) / (PI*x);
			//enhancement = pow(enhancement,2.0);
			enhancement = (ETA_TIMUR * sin(PI*x) / (PI*x)) + (ETA_TIMUR_prime * sin(2.0*PI*x) / (2.0*PI*x));
		}
	  //enhancement  *= (2.*PI*PI/3.) * (ETA_CANONIC - 1.);
		enhancement  += 1.0;
		enhancement   = sqrt(enhancement);
	  //enhancement   = 1.0;
		r             = RC_SMBH;
	}
	rho_chi  = (1. + pow((R_EARTH/core),alpha)) / (1. + pow((r/core),alpha));
	rho_chi  = RHO_CHI_SOLAR * /* rho_{CDM \odot} est exprimé en [GeV cm^{-3}] */
	pow((R_EARTH/r),gamma) * pow(rho_chi,((beta - gamma)/alpha));
	rho_chi /= RHO_CHI_0;

	return (rho_chi*enhancement);
}

/********************************************************************************************/
/********************************************************************************************/


// Nous calculons ici le rapport $\left\{ \frac{\rho_{\chi}}{\rho_0} \right\}$
// ou $\rho_{\chi}$ designe la masse volumique des neutralinos exprimee
// en [GeV cm^{-3}] et dependante des coordonnees cylindriques rr et z.
// La masse volumique de reference est $\rho_0$ = 1 [GeV cm^{-3}].
// Le resultat final est sans dimension et depend de DEUX variables :
// - Le rayon galactocentrique rr est exprime en [kpc].
// - La distance verticale z est egalement exprimee en [kpc].

double rapport_rho_chi_sur_rho_0_Einasto(double rr,double z)
{
	double r,rho_khi;

	r        = sqrt(rr*rr + z*z);
	rho_khi  = rhos_Ein * exp(-2.0/alpha_Ein * (pow(r/rs_Ein, alpha_Ein) - 1.0));
	rho_khi /= RHO_CHI_0;
	
	return rho_khi; //[NO UNIT]
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/

void DNPBAR_ON_DTPBAR_gaelle_read_file(struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	long i_pbar,i_channel;

    FILE *gaelle_file;
	FILE *gaelle_file_test;

    gaelle_file      = fopen(gaelle_file_name,"r");
	gaelle_file_test = fopen("./gaelle_files/gaelle_file_test.dat","w");

//  We clean the tables by putting the contents to 0.
	for (i_pbar=0;i_pbar<=N_x_pbar_scan;i_pbar++)
    {
		pt_Primary_Source_Term->X_PBAR[i_pbar] = 0.0;
		for (i_channel=0;i_channel<=number_channels;i_channel++)
		{
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[i_channel][i_pbar] = 0.0;
		}				
  	}


//  We read Gaelle's file on the antiproton multiplicity for various channels.
	for (i_pbar=0;i_pbar<N_x_pbar_scan;i_pbar++)
	{
		fscanf(gaelle_file,
		" %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&pt_Primary_Source_Term->X_PBAR                      [i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 1][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 2][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 3][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 4][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 5][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 6][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 7][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 8][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 9][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[10][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[11][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[12][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[13][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[14][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[15][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[16][i_pbar]);
	}


//  This is a test ! We print what has been read in Gaelle's file.
	for (i_pbar=0;i_pbar<N_x_pbar_scan;i_pbar++)
    {
		fprintf(gaelle_file_test,
		" %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e \n",
		pt_Primary_Source_Term->X_PBAR                      [i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 1][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 2][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 3][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 4][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 5][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 6][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 7][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 8][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[ 9][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[10][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[11][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[12][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[13][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[14][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[15][i_pbar],
		pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[16][i_pbar]);
	}
	fclose(gaelle_file);
	fclose(gaelle_file_test);
}

/********************************************************************************************/
/********************************************************************************************/

double dNpbar_on_dEpbar_primary_calculation(double mass_chi, int channel, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	int    i_pbar, i_scan_pbar;
	double xi_x_pbar, T_pbar, x_pbar;

	if ((channel <= 0) || (channel > number_channels))
	{
		printf("\n ERREUR ! \n Fonction : 'dNpbar_on_dEpbar_primary_calculation'  \n channel is out of its scope ! \n");
		exit (0);
	}

    for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
      pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar] = 0.0;
  	}

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		x_pbar = T_pbar/mass_chi;

		if (x_pbar < x_pbar_scan_min)
		{
			x_pbar = x_pbar_scan_min;
		}
		if ((x_pbar_scan_min <= x_pbar) && (x_pbar < x_pbar_scan_max))
		{
			xi_x_pbar = ((double)N_x_pbar_scan) * log(x_pbar/((double)x_pbar_scan_min)) / log(((double)x_pbar_scan_max)/((double)x_pbar_scan_min));
			if (xi_x_pbar >= ((double)N_x_pbar_scan))
			{
				printf("\n ERREUR ! \n Fonction : 'dNpbar_on_dEpbar_primary_calculation'  \n xi_x_pbar >= N_x_pbar_scan ! \n");
				exit (0);
			}
			i_scan_pbar = xi_x_pbar;

			pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  =  pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[channel][i_scan_pbar];
			pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar] += (pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[channel][i_scan_pbar+1] - pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL[channel][i_scan_pbar]) * (xi_x_pbar - (double)i_scan_pbar);
		}
		if (x_pbar_scan_max <= x_pbar)
		{
			pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  = 0.0;
		}	
	}
}

/********************************************************************************************/
/********************************************************************************************/

//	Attention! Cette fonction doit etre precedee des fonctions 'DNPBAR_ON_DTPBAR_gaelle_read_file' et 'dNpbar_on_dEpbar_primary_calculation'

void primary_source_calculation (double mass_chi, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	int i_pbar;

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
    {
		pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] = 0.0;
	}
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{			
		pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar]  = 0.5 * sigma_v_annihilation * pow((RHO_CHI_0/mass_chi),2.0);
		pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *= pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]; //[#pbar cm{-3} s{-1} GeV{-1}]
	}
}

/********************************************************************************************/
/********************************************************************************************/
