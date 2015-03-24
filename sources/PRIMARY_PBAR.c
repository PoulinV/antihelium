#include "PRIMARY_PBAR.h"

/********************************************************************************************/
/********************************************************************************************/
/*
* Dans ce module, nous calculons les coefficients du tableau pt_Pbar->BESSEL_PBAR_PRI_Epbar_i
* en fonction de l'energie CINETIQUE des antiprotons T_pbar et des coefficients de BESSEL i.
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
* Protection concernant le fait que n_vert doit etre inferieur ou egal a 1000
* mais pas superieur !
*/
  if (n_vert > 1000)
  {
    printf(" ATTENTION ! n_vert > 1000 !");
    exit (0);
  }
/*
* On remet a zero le tableau pt_Pbar->BESSEL_PBAR_PRI_Epbar_i.
* On pourrait d'ailleurs le remettre a zero plus loin.
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
*   Nous commencons par le calcul de l'integrale radiale q_pbar_primary_i_z.
*   \beq
*   q_pbar_primary_i_z = \frac{2}{J_{1}^{2}(\alpha_{i})} \times
*   \int_{0}^{1} \, x \, dx \, J_{0} \left( \alpha_{i} x \right) \,
*   \left\{ \frac{\rho_{\chi}}{\rho_0} \right\}^{2} \;\; .
*   \eeq
*
*   Cette integrale depend de i et de z et nous stockons les resultats
*   correspondant a l'ordre de Bessel i dans le tableau q_pbar_primary_i_z
*   qui ne depend que de z.
*
*   Nous entamons une boucle sur z suivie par une boucle sur r.
*/
    dx_vert = 1. / (double) (2*n_vert);
    x_vert  = 0.;
    for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
    {
      z_vert = pt_Propagation->E_DIFFUS * x_vert;
/*    z_vert est la coordonnee verticale exprimee en [kpc].
*/
      q_pbar_primary_i_z[i_vert] = 0.0;
      dx_rad = 1. / (double) (2*n_rad);
      x_rad  = 0.0;
      for (i_rad=0;i_rad<=(2*n_rad);i_rad++)
      {
        r_rad = R_GAL * x_rad;
/*      r_rad est le rayon galactocentrique exprime en [kpc].
*/
        if (i_rad==0 || i_rad==(2*n_rad)) {weight_SIMPSON_rad = 1./3.;}
        else {weight_SIMPSON_rad = (1. + (double)(i_rad % 2)) * 2. / 3.;}
		
		#if defined (WIMP_annihilation)
			q_pbar_primary_i_z[i_vert] += (x_rad * dx_rad * weight_SIMPSON_rad) * besselj0(alpha_i[i]*x_rad) * pow(rapport_rho_chi_sur_rho_0 (r_rad,z_vert),2.0); //[NO UNIT].
		#elif defined (WIMP_decay)	
        	q_pbar_primary_i_z[i_vert] += (x_rad * dx_rad * weight_SIMPSON_rad) * besselj0(alpha_i[i]*x_rad) * rapport_rho_chi_sur_rho_0 (r_rad,z_vert); //[NO UNIT].
		#else
			printf("Error! \n Function : 'calculation_BESSEL_PBAR_PRIMARY_Epbar_i' \n You have to specify in COMMON.h WIMP_annihilation or WIMP_decay \n");
			exit (0);
		#endif	
			
		x_rad += dx_rad;
      }
      q_pbar_primary_i_z[i_vert] *= 2. / pow(besselj1(alpha_i[i]),2.0);
/*    [NO UNIT].
*/
      x_vert += dx_vert;
    }
/*
*   Nous continuons ensuite par une boucle sur l'energie des antiprotons
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
	  
/*    Si est exprime en [kpc^{-1}].  */
      Si =
      sqrt(pow(2.0*alpha_i[i]/R_GAL,2) + pow(pt_Propagation->VENT_GALACTIQUE*CM_PAR_KPC/K_pbar,2));


/*    Abar_i est exprime en [cm s^{-1}].  */
      Abar_i  = pt_Propagation->VENT_GALACTIQUE;
	  Abar_i += 2.0*E_DISC*CM_PAR_KPC *
      ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) * v_pbar * (DENSITE_H_DISC + pow(4.,(2./3.))*1.0*DENSITE_HE_DISC));
      Abar_i += K_pbar * Si / CM_PAR_KPC / tanh(Si*pt_Propagation->E_DIFFUS/2.);

	  pt_Pbar->TABLE_Abar_i[i_pbar][i] = Abar_i;


/*    Il convient maintenant d'integrer sur la variable verticale x_vert
*     ainsi que sur la variable radiale x_rad.
*/
      dx_vert = 1. / (double) (2*n_vert);
      x_vert  = 0.;
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
      for (i_vert=0;i_vert<=(2*n_vert);i_vert++)
      {
        z_vert = pt_Propagation->E_DIFFUS * x_vert;
/*      z_vert est la coordonnee verticale exprimee en [kpc].
*/
        if (i_vert==0 || i_vert==(2*n_vert)) {weight_SIMPSON_vert = 1./3.;}
        else {weight_SIMPSON_vert = (1. + (double)(i_vert % 2)) * 2. / 3.;}
      

//	Calcul avec les spectres de Gaelle :
		
        pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] += (dx_vert * weight_SIMPSON_vert) *
        q_pbar_primary_i_z[i_vert] *
        pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *
        exp(- pt_Propagation->VENT_GALACTIQUE*z_vert*CM_PAR_KPC / (2.*K_pbar)) *
        sinh((Si/2.)*(pt_Propagation->E_DIFFUS-z_vert)) / sinh((Si/2.)*pt_Propagation->E_DIFFUS); //[antiprotons cm^{-3} s^{-1} GeV^{-1}].
						
		
		x_vert += dx_vert;
      }
      pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] *= 2. * pt_Propagation->E_DIFFUS*CM_PAR_KPC / Abar_i;
	  
	  
/*    S'exprime maintenant en unites de [antiprotons cm^{-3} GeV^{-1}].
*/
    }
  }
  return;
}

/********************************************************************************************/
/********************************************************************************************/
/*
 * Nous calculons ici le rapport $\left\{ \frac{\rho_{\chi}}{\rho_0} \right\}$
 * ou $\rho_{\chi}$ designe la masse volumique des neutralinos exprimee
 * en [GeV cm^{-3}] et dependante des coordonnees cylindriques rr et z.
 * La masse volumique de reference est $\rho_0$ = 1 [GeV cm^{-3}].
 * Le resultat final est sans dimension et depend de DEUX variables :
 * - Le rayon galactocentrique rr est exprime en [kpc].
 * - La distance verticale z est egalement exprimee en [kpc].
 *
 */
double rapport_rho_chi_sur_rho_0_old(double rr,double z)
{
	double alpha,beta,gamma,core;
	double r,x,rho_chi;
	double enhancement = 1.0;
/*
	double RC_SMBH      = 0.5;
	double RC_SMBH      = 0.001;
	double RC_SMBH      = 0.1; VALEUR UTILISEE PAR PIERRE BRUN DANS SON MONTE_CARLO
*/
	//double RC_SMBH      = 0.5;
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
	alpha = 2.0;
	beta  = 2.0;
	gamma = 0.0;
  //core  = 25.0; /* [kpc] */
  //core  = 20.0; /* [kpc] */
	core  = 4.38; /* [kpc] */

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
	rho_chi  = RHO_CHI_SOLAR * /* rho_{CDM \odot} est exprime en [GeV cm^{-3}] */
	pow((R_EARTH/r),gamma) * pow(rho_chi,((beta - gamma)/alpha));
	rho_chi /= RHO_CHI_0;

	return (rho_chi*enhancement);
}
/********************************************************************************************/
/********************************************************************************************/
/*
 * Nous calculons ici le rapport $\left\{ \frac{\rho_{\chi}}{\rho_0} \right\}$
 * ou $\rho_{\chi}$ designe la masse volumique des neutralinos exprimee
 * en [GeV cm^{-3}] et dependante des coordonnees cylindriques rr et z.
 * La masse volumique de reference est $\rho_0$ = 1 [GeV cm^{-3}].
 * Le resultat final est sans dimension et depend de DEUX variables :
 * - Le rayon galactocentrique rr est exprime en [kpc].
 * - La distance verticale z est egalement exprimee en [kpc].
 *
 */
double rapport_rho_chi_sur_rho_0(double rr,double z)
{

  	double alpha,beta,gamma,core;
  	double r,x,rho_chi;
  	double enhancement, a1, a2;
	
	r = sqrt(rr*rr + z*z);

	// NFW
	#if defined (NFW)  

  		alpha = 1.0;						// [NO UNIT]
  		beta  = 3.0;						// [NO UNIT]
  		gamma = 1.0;						// [NO UNIT]
  		core  = 24.42;						// [kpc] 
	
		enhancement = 1.0;

		a2 = 8.0 * gamma * (pow(PI,(2.0)) - 9.0 + 6.0 * gamma) / (9.0 * (3.0 - 2.0*gamma));
		a1 = a2 + (2.0*gamma);

  		if (r <= RC_SMBH)
  		{
    		x = r / RC_SMBH;

			if (x <= 1.0e-3)
			{
				enhancement = a1 + a2;
			}
    		else
    		{
				enhancement = (a1 * sin(PI*x) / (PI*x)) + (a2 * sin(2.0*PI*x) / (2.0*PI*x));
    		}
		
    		enhancement += 1.0;
    		enhancement  = sqrt(enhancement);

    		r = RC_SMBH;
  		}
	

		// Generic Profile 

	  	rho_chi = (1. + pow((R_EARTH/core),alpha)) / (1. + pow((r/core),alpha));
	  	rho_chi = RHO_CHI_SOLAR * pow((R_EARTH/r),gamma) * pow(rho_chi,((beta - gamma)/alpha));		// [GeV cm^{-3}]
		rho_chi *= enhancement;
		rho_chi /= RHO_CHI_0;
		
	  	return rho_chi;
	
	// MOORE
	#elif defined(moore)

	  	alpha = 1.0;						// [NO UNIT]
	  	beta  = 3.0;						// [NO UNIT]
	  	gamma = 1.16;						// [NO UNIT]
	  	core  = 30.28;						// [kpc] 
	
		enhancement = 1.0;

		a2 = 8.0 * gamma * (pow(PI,(2.0)) - 9.0 + 6.0 * gamma) / (9.0 * (3.0 - 2.0*gamma));
		a1 = a2 + (2.0*gamma);

	  	if (r <= RC_SMBH)
	  	{
	    	x = r / RC_SMBH;

			if (x <= 1.0e-3)
			{
				enhancement = a1 + a2;
			}
	    	else
	    	{
				enhancement = (a1 * sin(PI*x) / (PI*x)) + (a2 * sin(2.0*PI*x) / (2.0*PI*x));
	    	}
		
	    	enhancement += 1.0;
	    	enhancement  = sqrt(enhancement);

	    	r = RC_SMBH;
	  	}


		// Generic Profile

	  	rho_chi = (1. + pow((R_EARTH/core),alpha)) / (1. + pow((r/core),alpha));
	  	rho_chi = RHO_CHI_SOLAR * pow((R_EARTH/r),gamma) * pow(rho_chi,((beta - gamma)/alpha));			// [GeV cm^{-3}]
		rho_chi *= enhancement;
		rho_chi /= RHO_CHI_0;
	  	
		return rho_chi;
	

	// EINASTO
	#elif defined(einasto)
	
		alpha = 0.17;							// [NO UNIT]
		core = 28.44;							// [kpc]
	
		rho_chi = pow(R_EARTH/core, alpha) - pow(r/core, alpha);
		rho_chi = RHO_CHI_SOLAR * exp((2.0/alpha) * rho_chi);
		rho_chi /= RHO_CHI_0;
	
		return rho_chi;

	
	
	// EINASTO B
	#elif defined(einastoB)

		alpha = 0.11;							// [NO UNIT]
		core = 35.24;							// [kpc]
	
		rho_chi = pow(R_EARTH/core, alpha) - pow(r/core, alpha);
		rho_chi = RHO_CHI_SOLAR * exp((2.0/alpha) * rho_chi);
		rho_chi /= RHO_CHI_0;
	
		return rho_chi;
	

	//	ISOTHERMAL
	#elif defined(isothermal)
	
	  	alpha = 2.0;							// [NO UNIT]
	  	beta  = 2.0;							// [NO UNIT]
		gamma = 0.0;
	  	core  = 4.38;							// [kpc] 


		// Generic Profile 

	  	rho_chi = (1. + pow((R_EARTH/core),alpha)) / (1. + pow((r/core),alpha));
	  	rho_chi = RHO_CHI_SOLAR * pow((R_EARTH/r),gamma) * pow(rho_chi,((beta - gamma)/alpha));
		rho_chi /= RHO_CHI_0;
	
		return rho_chi;
	
	
	//	BURKERT
	#elif defined(burkert)
	
		core = 12.67;
	
		rho_chi = (1.0 + r/core) * (1.0 + pow(r/core, 2.0));
		rho_chi = RHO_CHI_SOLAR * (1.0 + R_EARTH/core) * (1.0 + pow(R_EARTH/core, 2.0)) / rho_chi;
		rho_chi /= RHO_CHI_0;
	
		return rho_chi;

	#else
		
		printf("\n ERROR function 'rapport_rho_chi_sur_rho_0' \n No DM profile defined \n\n ");
		exit(0);

	#endif 	
}

/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/

void DNPBAR_ON_DTPBAR_gaelle_read_file(double mass_chi, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	long i_pbar,i_channel;
	long i_mass;
	char* name_gaelle_file_m_inf;
	char* name_gaelle_file_m_sup;
	double mass_min, mass_max;
	
	name_gaelle_file_m_inf = NULL;
	name_gaelle_file_m_sup = NULL;
	
    FILE *gaelle_file_m_inf;
	FILE *gaelle_file_m_inf_test;
    FILE *gaelle_file_m_sup;
	FILE *gaelle_file_m_sup_test;
	
	mass_min = pt_Primary_Source_Term->GAELLE_MASSES[1];
	mass_max = pt_Primary_Source_Term->GAELLE_MASSES[N_gaelle_masses];
	
	if((mass_chi < mass_min) || (mass_chi > mass_max))
	{
		printf("\n ERREUR ! \n Fonction : 'DNPBAR_ON_DTPBAR_gaelle_read_file'  \n mass_chi is out of range ! \n");
		exit (0);
	}
	if ((mass_chi >= mass_min) && (mass_chi <= mass_max))
	{
		i_mass = 1;
		
		while ((mass_chi >= pt_Primary_Source_Term->GAELLE_MASSES[i_mass]) && (i_mass < N_gaelle_masses))
		{
			pt_Primary_Source_Term->mass_inf = pt_Primary_Source_Term->GAELLE_MASSES[i_mass];
			name_gaelle_file_m_inf = pt_Primary_Source_Term->GAELLE_FILES_NAME[i_mass];
			
			pt_Primary_Source_Term->mass_sup = pt_Primary_Source_Term->GAELLE_MASSES[i_mass + 1];
			name_gaelle_file_m_sup = pt_Primary_Source_Term->GAELLE_FILES_NAME[i_mass + 1];
			
			//printf("i_mass = %d \n", i_mass);
			
			i_mass++;
		}	
	}
		
	gaelle_file_m_inf = fopen(name_gaelle_file_m_inf,"r");
	gaelle_file_m_sup = fopen(name_gaelle_file_m_sup,"r");
	
	gaelle_file_m_inf_test = fopen("../sources/FluxSourceAntip/gaelle_file_m_inf_test.dat","w");
	gaelle_file_m_sup_test = fopen("../sources/FluxSourceAntip/gaelle_file_m_sup_test.dat","w");
	
	
//  We clean the tables by putting the contents to 0.
	for (i_pbar=0;i_pbar<=N_x_pbar_scan;i_pbar++)
    {
		pt_Primary_Source_Term->X_PBAR_M_INF[i_pbar] = 0.0;
		pt_Primary_Source_Term->X_PBAR_M_SUP[i_pbar] = 0.0;
		for (i_channel=0;i_channel<=number_channels;i_channel++)
		{
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[i_channel][i_pbar] = 0.0;
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[i_channel][i_pbar] = 0.0;
		}				
  	}
			
	
//  We read Gaelle's file on the antiproton multiplicity for various channels.
	for (i_pbar=0;i_pbar<N_x_pbar_scan;i_pbar++)
	{
		fscanf(gaelle_file_m_inf,
		" %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&pt_Primary_Source_Term->X_PBAR_M_INF                      [i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 1][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 2][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 3][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 4][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 5][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 6][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 7][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 8][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 9][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[10][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[11][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[12][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[13][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[14][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[15][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[16][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[17][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[18][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[19][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[20][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[21][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[22][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[23][i_pbar]);
	
		fscanf(gaelle_file_m_sup,
		" %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&pt_Primary_Source_Term->X_PBAR_M_SUP                     [i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 1][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 2][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 3][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 4][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 5][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 6][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 7][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 8][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 9][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[10][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[11][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[12][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[13][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[14][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[15][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[16][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[17][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[18][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[19][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[20][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[21][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[22][i_pbar],
		&pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[23][i_pbar]);
						
		}	
		
	
//  This is a test ! We print what has been read in Gaelle's file.
		for (i_pbar=0;i_pbar<N_x_pbar_scan;i_pbar++)
	    {
			fprintf(gaelle_file_m_inf_test,
			" %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e \n",
			pt_Primary_Source_Term->X_PBAR_M_INF                      [i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 1][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 2][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 3][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 4][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 5][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 6][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 7][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 8][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[ 9][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[10][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[11][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[12][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[13][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[14][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[15][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[16][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[17][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[18][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[19][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[20][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[21][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[22][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[23][i_pbar]);
			
			fprintf(gaelle_file_m_sup_test,
			" %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e\t %.15e \n",
			pt_Primary_Source_Term->X_PBAR_M_SUP                      [i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 1][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 2][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 3][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 4][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 5][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 6][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 7][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 8][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[ 9][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[10][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[11][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[12][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[13][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[14][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[15][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[16][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[17][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[18][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[19][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[20][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[21][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[22][i_pbar],
			pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[23][i_pbar]);
		}
	

	fclose(gaelle_file_m_inf);
	fclose(gaelle_file_m_inf_test);
	fclose(gaelle_file_m_sup);
	fclose(gaelle_file_m_sup_test);
				

}

/********************************************************************************************/
/********************************************************************************************/

double dNpbar_on_dEpbar_primary_calculation(double mass_chi, int channel, struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	int    i_pbar, i_scan_pbar;
	double xi_x_pbar, T_pbar, x_pbar;
	double dNpbar_on_dEpbar_xi_m, dNpbar_on_dEpbar_xi_plus_un_m;

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
			if (xi_x_pbar > ((double)N_x_pbar_scan))
			{
				printf("\n ERREUR ! \n Fonction : 'dNpbar_on_dEpbar_primary_calculation'  \n xi_x_pbar >= N_x_pbar_scan ! \n");
				exit (0);
			}
			else if (xi_x_pbar == ((double)N_x_pbar_scan))
			{
				pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  = pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][N_x_pbar_scan];
				pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar] += (pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[channel][N_x_pbar_scan] - pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][N_x_pbar_scan]) * (mass_chi - pt_Primary_Source_Term->mass_inf) / (pt_Primary_Source_Term->mass_sup - pt_Primary_Source_Term->mass_inf);
			}
			else
			{
				i_scan_pbar = xi_x_pbar;
				
				dNpbar_on_dEpbar_xi_m  = pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][i_scan_pbar];
				dNpbar_on_dEpbar_xi_m += (pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[channel][i_scan_pbar] - pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][i_scan_pbar]) * (mass_chi - pt_Primary_Source_Term->mass_inf) / (pt_Primary_Source_Term->mass_sup - pt_Primary_Source_Term->mass_inf);
				
				
				dNpbar_on_dEpbar_xi_plus_un_m  = pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][i_scan_pbar+1];
				dNpbar_on_dEpbar_xi_plus_un_m += (pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_SUP[channel][i_scan_pbar+1] - pt_Primary_Source_Term->DNPBAR_ON_DTPBAR_CHANNEL_M_INF[channel][i_scan_pbar+1]) * (mass_chi - pt_Primary_Source_Term->mass_inf) / (pt_Primary_Source_Term->mass_sup - pt_Primary_Source_Term->mass_inf);
				

				pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  =  dNpbar_on_dEpbar_xi_m;
				pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar] += (dNpbar_on_dEpbar_xi_plus_un_m - dNpbar_on_dEpbar_xi_m) * (xi_x_pbar - (double)i_scan_pbar);
			
			}
		}
		if (x_pbar_scan_max <= x_pbar)
		{
			pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]  = 0.0;
		}	
	}
	return (1.0);
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
	#if defined (WIMP_annihilation)
		for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
		{			
			pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar]  = 0.5 * sigma_v_annihilation * pow((RHO_CHI_0/mass_chi),2.0);
			pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *= pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]; //[#pbar cm{-3} s{-1} GeV{-1}]
		}
	#elif defined (WIMP_decay)
		for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
		{			
			pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar]  =  decay_rate *(RHO_CHI_0/mass_chi);
			pt_Primary_Source_Term->PRIMARY_SOURCE_TERM[i_pbar] *= pt_Primary_Source_Term->DNPBAR_ON_DEPBAR[i_pbar]; //[#pbar cm{-3} s{-1} GeV{-1}]
		}
	#else
		printf("Error! \n Function : 'primary_source_calculation'\n You have to specify in COMMON.h WIMP_annihilation or WIMP_decay \n");
		exit (0);
	#endif
}
	
/********************************************************************************************/
/********************************************************************************************/

//	Nous remplissons les tableaux GAELLE_MASSES[N_gaelle_masses+1] et GAELLE_FILES_NAME[N_gaelle_masses+1].

void DM_preliminary(struct Structure_Primary_Source_Term* pt_Primary_Source_Term)
{
	int i;
	
	for(i=0;i<=N_gaelle_masses;i++)
	{
		pt_Primary_Source_Term->GAELLE_MASSES[i] = 0.0;
		pt_Primary_Source_Term->GAELLE_FILES_NAME[i] = NULL;
	}

	pt_Primary_Source_Term->GAELLE_MASSES[1] = 5.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[1] = "../sources/FluxSourceAntip/mDM=5GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[2] = 6.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[2] = "../sources/FluxSourceAntip/mDM=6GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[3] = 8.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[3] = "../sources/FluxSourceAntip/mDM=8GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[4] = 10.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[4] = "../sources/FluxSourceAntip/mDM=10GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[5] = 15.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[5] = "../sources/FluxSourceAntip/mDM=15GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[6] = 20.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[6] = "../sources/FluxSourceAntip/mDM=20GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[7] = 25.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[7] = "../sources/FluxSourceAntip/mDM=25GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[8] = 30.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[8] = "../sources/FluxSourceAntip/mDM=30GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[9] = 40.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[9] = "../sources/FluxSourceAntip/mDM=40GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[10] = 50.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[10] = "../sources/FluxSourceAntip/mDM=50GeV.dat";
	
	pt_Primary_Source_Term->GAELLE_MASSES[11] = 60.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[11] = "../sources/FluxSourceAntip/mDM=60GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[12] = 70.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[12] = "../sources/FluxSourceAntip/mDM=70GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[13] = 80.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[13] = "../sources/FluxSourceAntip/mDM=80GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[14] = 90.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[14] = "../sources/FluxSourceAntip/mDM=90GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[15] = 100.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[15] = "../sources/FluxSourceAntip/mDM=100GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[16] = 110.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[16] = "../sources/FluxSourceAntip/mDM=110GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[17] = 120.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[17] = "../sources/FluxSourceAntip/mDM=120GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[18] = 130.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[18] = "../sources/FluxSourceAntip/mDM=130GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[19] = 140.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[19] = "../sources/FluxSourceAntip/mDM=140GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[20] = 150.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[20] = "../sources/FluxSourceAntip/mDM=150GeV.dat";
	
	pt_Primary_Source_Term->GAELLE_MASSES[21] = 160.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[21] = "../sources/FluxSourceAntip/mDM=160GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[22] = 180.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[22] = "../sources/FluxSourceAntip/mDM=180GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[23] = 200.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[23] = "../sources/FluxSourceAntip/mDM=200GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[24] = 220.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[24] = "../sources/FluxSourceAntip/mDM=220GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[25] = 240.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[25] = "../sources/FluxSourceAntip/mDM=240GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[26] = 260.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[26] = "../sources/FluxSourceAntip/mDM=260GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[27] = 280.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[27] = "../sources/FluxSourceAntip/mDM=280GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[28] = 300.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[28] = "../sources/FluxSourceAntip/mDM=300GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[29] = 330.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[29] = "../sources/FluxSourceAntip/mDM=330GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[30] = 360.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[30] = "../sources/FluxSourceAntip/mDM=360GeV.dat";
	
	pt_Primary_Source_Term->GAELLE_MASSES[31] = 400.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[31] = "../sources/FluxSourceAntip/mDM=400GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[32] = 450.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[32] = "../sources/FluxSourceAntip/mDM=450GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[33] = 500.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[33] = "../sources/FluxSourceAntip/mDM=500GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[34] = 550.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[34] = "../sources/FluxSourceAntip/mDM=550GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[35] = 600.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[35] = "../sources/FluxSourceAntip/mDM=600GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[36] = 650.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[36] = "../sources/FluxSourceAntip/mDM=650GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[37] = 700.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[37] = "../sources/FluxSourceAntip/mDM=700GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[38] = 750.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[38] = "../sources/FluxSourceAntip/mDM=750GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[39] = 800.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[39] = "../sources/FluxSourceAntip/mDM=800GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[40] = 900.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[40] = "../sources/FluxSourceAntip/mDM=900GeV.dat";
	
	pt_Primary_Source_Term->GAELLE_MASSES[41] = 1000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[41] = "../sources/FluxSourceAntip/mDM=1000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[42] = 1100.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[42] = "../sources/FluxSourceAntip/mDM=1100GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[43] = 1200.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[43] = "../sources/FluxSourceAntip/mDM=1200GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[44] = 1300.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[44] = "../sources/FluxSourceAntip/mDM=1300GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[45] = 1500.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[45] = "../sources/FluxSourceAntip/mDM=1500GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[46] = 1700.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[46] = "../sources/FluxSourceAntip/mDM=1700GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[47] = 2000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[47] = "../sources/FluxSourceAntip/mDM=2000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[48] = 2500.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[48] = "../sources/FluxSourceAntip/mDM=2500GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[49] = 3000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[49] = "../sources/FluxSourceAntip/mDM=3000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[50] = 4000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[50] = "../sources/FluxSourceAntip/mDM=4000GeV.dat";
	
	pt_Primary_Source_Term->GAELLE_MASSES[51] = 5000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[51] = "../sources/FluxSourceAntip/mDM=5000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[52] = 6000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[52] = "../sources/FluxSourceAntip/mDM=6000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[53] = 7000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[53] = "../sources/FluxSourceAntip/mDM=7000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[54] = 8000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[54] = "../sources/FluxSourceAntip/mDM=8000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[55] = 9000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[55] = "../sources/FluxSourceAntip/mDM=9000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[56] = 10000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[56] = "../sources/FluxSourceAntip/mDM=10000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[57] = 12000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[57] = "../sources/FluxSourceAntip/mDM=12000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[58] = 15000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[58] = "../sources/FluxSourceAntip/mDM=15000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[59] = 20000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[59] = "../sources/FluxSourceAntip/mDM=20000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[60] = 30000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[60] = "../sources/FluxSourceAntip/mDM=30000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[61] = 50000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[61] = "../sources/FluxSourceAntip/mDM=50000GeV.dat";
	pt_Primary_Source_Term->GAELLE_MASSES[62] = 100000.0;
	pt_Primary_Source_Term->GAELLE_FILES_NAME[62] = "../sources/FluxSourceAntip/mDM=100000GeV.dat";
	
}

/********************************************************************************************/
/********************************************************************************************/

void primary_spectra_BCGS_2014(struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation, struct Structure_Primary_Source_Term* pt_Primary_Source_Term, double alpha_i[NDIM+1])
{
	char file_name[1024], directory_name[1024];
	
	int channel, i_channel;
	int i_mass_chi;
	
	long i_data,i_iteration,i_pbar,i;
	
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double mass_chi;
	double FLUX_PBAR_MIN,    FLUX_PBAR_MED,    FLUX_PBAR_MAX;

	static double PBAR_SPECTRUM_MIN       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MED       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MAX       [DIM_TAB_PBAR+1];

	FILE* results;
	
	results = NULL;
	
	DM_preliminary(pt_Primary_Source_Term);

	//mass_chi = mass_chi_choice;
	channel  = channel_choice;
	
	results = fopen(pbar_IS_spectrum_file_name,"w");
	
	//fprintf(results, " \n\n channel = %d \n\n", channel);

// Boucle sur les canaux d'annihilation	
	for(i_channel=1;i_channel<=number_channels;i_channel++)
	{
		channel = i_channel;
		
		printf(" \n\n channel = %d \n\n", channel);
		
		sprintf(directory_name, "./results/mass_scan_wc/channel_%d", channel);
		mkdir(directory_name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

// Boucle sur la masse des WIMPs	
		for(i_mass_chi=1;i_mass_chi<=N_gaelle_masses;i_mass_chi++)	
		{
			#if defined (WIMP_annihilation)
				mass_chi = pt_Primary_Source_Term->GAELLE_MASSES[i_mass_chi];
			#elif defined (WIMP_decay)
				mass_chi = 2.0 * pt_Primary_Source_Term->GAELLE_MASSES[i_mass_chi];
			#else
				printf("Error! \n Function : 'main' \n You have to specify in COMMON.h WIMP_annihilation or WIMP_decay \n");
				exit (0);
			#endif	
		
			printf("\n mass_chi = %.2e GeV \n", mass_chi);
		
			sprintf(file_name, "./results/mass_scan_wc/channel_%d/mDM=%gGeV.txt", channel, mass_chi);
			results = fopen(file_name,"w");
		
			#if defined (WIMP_annihilation)
				DNPBAR_ON_DTPBAR_gaelle_read_file   (mass_chi, pt_Primary_Source_Term);
				dNpbar_on_dEpbar_primary_calculation(mass_chi, channel, pt_Primary_Source_Term);
				primary_source_calculation          (mass_chi, pt_Primary_Source_Term);
			#elif defined (WIMP_decay)
				DNPBAR_ON_DTPBAR_gaelle_read_file   (mass_chi/2.0, pt_Primary_Source_Term);
				dNpbar_on_dEpbar_primary_calculation(mass_chi/2.0, channel, pt_Primary_Source_Term);
				primary_source_calculation          (mass_chi, pt_Primary_Source_Term);
			#else
				printf("Error! \n Function : 'main' \n You have to specify in COMMON.h WIMP_annihilation or WIMP_decay \n");
				exit (0);
			#endif	
				
				
///////////////////////////////////////////////////////////////////////////////////////////
//	ON PEUT Y ALLER !
///////////////////////////////////////////////////////////////////////////////////////////

			for (i_data=1;i_data<=3;i_data++)
			//for (i_data=2;i_data<=2;i_data++)
			{
//			Nous definissons a ce niveau les parametres que FIORENZA, DAVID et RICHARD
//			-- hereafter called FDR -- ont determines.

				if (i_data == 1)       // CAS MAX 
				{
					pt_Propagation->PUISSANCE_COEFF_DIFF = 0.46;
					pt_Propagation->DIFFUSION_0_GV  = 0.0765 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
					pt_Propagation->E_DIFFUS        = 15.0;                                       		// [kpc]
					pt_Propagation->VENT_GALACTIQUE = (5.0   * 1.0e5);                            		// [cm s^{-1}]
					pt_Propagation->V_ALFEN         = (117.6 * 1.0e5);                            		// [cm s^{-1}]
				}
				else if (i_data == 2) // CAS MED 
				{
					pt_Propagation->PUISSANCE_COEFF_DIFF = 0.70;
					pt_Propagation->DIFFUSION_0_GV  = 0.0112 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
					pt_Propagation->E_DIFFUS        = 4.0;                                        		// [kpc]
					pt_Propagation->VENT_GALACTIQUE = (12.0  * 1.0e5);                            		// [cm s^{-1}]
					pt_Propagation->V_ALFEN         = (52.9  * 1.0e5);                            		// [cm s^{-1}]
				}
				else if (i_data == 3) // CAS MIN 
				{
					pt_Propagation->PUISSANCE_COEFF_DIFF = 0.85;
					pt_Propagation->DIFFUSION_0_GV  = 0.0016 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
					pt_Propagation->E_DIFFUS        = 1.0;                                        		// [kpc]
					pt_Propagation->VENT_GALACTIQUE = (13.5  * 1.0e5);                            		// [cm s^{-1}]
					pt_Propagation->V_ALFEN         = (22.4  * 1.0e5);                            		// [cm s^{-1}]
				}

//				Nous imprimons les coefficients de diffusion_propagation choisis dans le calcul.
/*		
				printf(" CAS NUMERO      = %ld \n",i_data);
				printf(" DELTA           = %.5e [NO UNIT]\n",pt_Propagation->PUISSANCE_COEFF_DIFF);
				printf(" DIFFUSION_0_GV  = %.5e [cm^{2} s^{-1}]\n",pt_Propagation->DIFFUSION_0_GV);
				printf(" E_DIFFUS        = %.5e [kpc]\n",pt_Propagation->E_DIFFUS);
				printf(" VENT_GALACTIQUE = %.5e [cm s^{-1}]\n",pt_Propagation->VENT_GALACTIQUE);
				printf(" V_ALFEN         = %.5e [cm s^{-1}]\n\n",pt_Propagation->V_ALFEN);
*/		
//				On remet a zero les tableaux pt_Pbar->BESSEL_PBAR_SEC_Epbar_i et pt_Pbar->BESSEL_PBAR_TER_Epbar_i.

				for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
				{
					for (i=0;i<=NDIM;i++)
					{
						pt_Pbar->BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
						pt_Pbar->BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] = 0.0;
						pt_Pbar->BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = 0.0;
					}
				}

//				CALCUL DE LA CONTRIBUTION PRIMAIRE PROVENANT DE L'ANNIHILATION DES NEUTRALINOS.

				calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,500, alpha_i, pt_Pbar, pt_Propagation, pt_Primary_Source_Term);
				//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);

//				CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE
//				PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.

				//calculation_BESSEL_PROTON_Ep_i(alpha_i, &Proton, &Propagation);
				//calculation_BESSEL_HELIUM_Ep_i(alpha_i, &Helium, &Propagation);
				//calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation);
				calculation_BESSEL_PBAR_SUM_123_Epbar_i(pt_Pbar);

//				CALCUL DU SPECTRE FINAL DES ANTIPROTONS.

				//goto TEST;
				for (i_iteration=1;i_iteration<=5;i_iteration++)
				{
					calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, pt_Pbar);
					calculation_BESSEL_PBAR_SUM_123_Epbar_i(pt_Pbar);
				}

				//goto TEST;
				for (i_iteration=1;i_iteration<=5;i_iteration++)
				{
					calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, pt_Pbar);
					calculation_BESSEL_PBAR_TOT_direct_inversion_A(pt_Pbar, pt_Propagation);
					//calculation_BESSEL_PBAR_TOT_direct_inversion_B(&Pbar, &Propagation);
					//calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(&Pbar, &Propagation);
					//calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(15., 1200., &Pbar, &Propagation);
				}

//			We compute now the antiproton spectrum and store for each KINETIC
//			ENERGY the lowest -- PBAR_SPECTRUM_MIN -- the medium
//			-- PBAR_SPECTRUM_MED -- and the largest -- PBAR_SPECTRUM_MAX --
//			values which we meet.

TEST:

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

					if      (i_data == 1){PBAR_SPECTRUM_MAX[i_pbar] = flux_antiproton_IS;}
					else if (i_data == 2){PBAR_SPECTRUM_MED[i_pbar] = flux_antiproton_IS;}
					else if (i_data == 3){PBAR_SPECTRUM_MIN[i_pbar] = flux_antiproton_IS;}
				}
			}
		
		
////////////////////////////////////////////////////////////////////////////////////////////
//			On imprime le resultat

		
			for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
			{
				T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
				E_pbar_IS = T_pbar_IS + MASSE_PROTON;

				FLUX_PBAR_MIN = PBAR_SPECTRUM_MIN[i_pbar];
				FLUX_PBAR_MED = PBAR_SPECTRUM_MED[i_pbar];
				FLUX_PBAR_MAX = PBAR_SPECTRUM_MAX[i_pbar];

				fprintf(results, " %.10e\t %.10e\t %.10e\t %.10e\t \n", T_pbar_IS, (1.0e04*FLUX_PBAR_MIN), (1.0e04*FLUX_PBAR_MED), (1.0e04*FLUX_PBAR_MAX));	
			}
			fclose(results);
		}
	}
}

/********************************************************************************************/
/********************************************************************************************/




