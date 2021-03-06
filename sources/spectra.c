#include "spectra.h"

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

void print_total_pbar_spectra_MIN_MED_MAX(struct Structure_Nuclei* pt_Proton, struct Structure_Nuclei* pt_Helium, struct Structure_Pbar* pt_Pbar, struct Structure_Cross_Section* pt_Cross_Section, struct Structure_Propagation* pt_Propagation, struct Structure_Primary_Source_Term* pt_Primary_Source_Term, double alpha_i[NDIM+1])
{
	long   i_data,i_iteration,i_pbar,i;

	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA,flux_proton_TOA;

	double FLUX_PBAR_MIN, FLUX_PBAR_MED, FLUX_PBAR_MAX;
	double FLUX_PBAR_TOA_MIN,FLUX_PBAR_TOA_MED,FLUX_PBAR_TOA_MAX;

	int    channel, i_channel;
	double mass_chi;
	int i_mass_chi;

	static double PBAR_SPECTRUM_MIN       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MED       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MAX       [DIM_TAB_PBAR+1];
	static double RESULTS_T_PBAR_TOA           [DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MIN[DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MED[DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MAX[DIM_TAB_PBAR+1];

	FILE* results;


	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
	    PBAR_SPECTRUM_MIN[i_pbar] = 0.0;
	    PBAR_SPECTRUM_MED[i_pbar] = 0.0;
	    PBAR_SPECTRUM_MAX[i_pbar] = 0.0;
	}



///////////////////////////////////////////////////////////////////////////////////////////
//	On charge les valeurs des multiplicites
//	$g \left( T_{\pbar} \right) = \frac{d N_{\pbar}}{d T_{\pbar}}$ et l'on calcule

//	\beq
// 	result \; = \; \frac{1}{2} \, < \sigma v >_{channel f} \,
// 	\left\{ \frac{\rho_{0}}{m_{\chi}} \right\}^{2} \, \frac{d N_{\pbar}}{d T_{\pbar}}_{channel f} \;\; .
// 	\eeq

// 	que l'on stocke dans le tableau PRIMARY_SOURCE_TERM[DIM_TAB_PBAR+1].

	DM_preliminary(pt_Primary_Source_Term);

	mass_chi = mass_chi_choice;
	channel  = channel_choice;



	printf("\n mass_chi = %.2e GeV \n", mass_chi);

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
//				On remet a zero les tableaux Pbar.BESSEL_PBAR_SEC_Epbar_i et Pbar.BESSEL_PBAR_TER_Epbar_i.

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
		//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, pt_Pbar, pt_Propagation, pt_Primary_Source_Term);

//				CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE
//				PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.

		calculation_BESSEL_PROTON_Ep_i(alpha_i, pt_Proton, pt_Propagation);
		calculation_BESSEL_HELIUM_Ep_i(alpha_i, pt_Helium, pt_Propagation);
		calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, pt_Proton, pt_Helium, pt_Pbar, pt_Cross_Section, pt_Propagation);
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
			//calculation_BESSEL_PBAR_TOT_direct_inversion_B(pt_Pbar, pt_Propagation);
			//calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(pt_Pbar, pt_Propagation);
			//calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(15., 1200., pt_Pbar, pt_Propagation);
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
			//flux_antiproton_IS = GENERIC_FLUX(R_EARTH,0.,E_pbar_IS,MASSE_PROTON,1.,alpha_i,pt_Pbar->BESSEL_PBARi, pt_Propagation);

			if      (i_data == 1){PBAR_SPECTRUM_MAX[i_pbar] = flux_antiproton_IS;}
			else if (i_data == 2){PBAR_SPECTRUM_MED[i_pbar] = flux_antiproton_IS;}
			else if (i_data == 3){PBAR_SPECTRUM_MIN[i_pbar] = flux_antiproton_IS;}
		}
	}


////////////////////////////////////////////////////////////////////////////////////////////
//			On imprime le resultat

	results = fopen(pbar_IS_spectra_MIN_MED_MAX_file_name,"w");

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


////////////////////////////////////////////////////////////////////////////////////////////
//		Nous modulons maintenant les spectres PBAR obtenus.


	pt_Propagation->PHI_FISK = phi_fisk;

	results = fopen(pbar_TOA_spectra_MIN_MED_MAX_file_name,"w");


	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

		FFA_IS_to_TOA(1.,1.,pt_Propagation->PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MIN[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MIN);

		if (E_pbar_TOA <= MASSE_PROTON)
		{
			RESULTS_T_PBAR_TOA[i_pbar]       = 0.0;
			RESULTS_SPECTRUM_TOA_MAX[i_pbar] = 0.0;
			RESULTS_SPECTRUM_TOA_MED[i_pbar] = 0.0;
			RESULTS_SPECTRUM_TOA_MIN[i_pbar] = 0.0;
			//continue;
		}

		FFA_IS_to_TOA(1.,1.,pt_Propagation->PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MED[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MED);
		FFA_IS_to_TOA(1.,1.,pt_Propagation->PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MAX[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MAX);

//		Nous les imprimons !

		T_pbar_TOA = E_pbar_TOA - MASSE_PROTON;
		fprintf(results, " %.10e\t %.10e\t %.10e\t %.10e\t \n", T_pbar_TOA, (1.0e04*FLUX_PBAR_TOA_MIN), (1.0e04*FLUX_PBAR_TOA_MED), (1.0e04*FLUX_PBAR_TOA_MAX));

//		Nous les stockons en memoire dans les tableaux RESULTS_T_PBAR_TOA[DIM_TAB_PBAR+1] et RESULTS_SPECTRUM_TOA_MIN_MED_MAX[DIM_TAB_PBAR+1];

		RESULTS_T_PBAR_TOA[i_pbar] = T_pbar_TOA;
		RESULTS_SPECTRUM_TOA_MAX[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MAX);
		RESULTS_SPECTRUM_TOA_MED[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MED);
		RESULTS_SPECTRUM_TOA_MIN[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MIN);
	}
	fclose(results);
}

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

//	On imprime le resultat

void print_PBAR_IS_SPECTRUM(double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], double A_nuclei)
{
	long i_pbar;
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double flux_pbar;

	FILE* results;
	// results = fopen(pbar_IS_spectrum_file_name,"w");
	if(A_nuclei == 1){
		results = fopen(pbar_IS_spectrum_file_name,"w");
		printf("Printing in file: %s\n", pbar_IS_spectrum_file_name);
	}
	else if (A_nuclei == 2){
		results = fopen(Debar_IS_spectrum_file_name,"w");
		printf("Printing in file: %s\n", Debar_IS_spectrum_file_name);
	}
	else if (A_nuclei == 3){
		results = fopen(He3bar_IS_spectrum_file_name,"w");
		printf("Printing in file: %s\n", He3bar_IS_spectrum_file_name);
	}
	else if (A_nuclei == 4){
		results = fopen(He4bar_IS_spectrum_file_name,"w");
		printf("Printing in file: %s\n", He4bar_IS_spectrum_file_name);
	}
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

		flux_pbar = PBAR_IS_SPECTRUM[i_pbar];

		fprintf(results, " %.10e\t %.10e\t \n", T_pbar_IS, (1.0e04*flux_pbar));											// [#pbar m^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

//	On imprime le resultat

void print_PBAR_TOA_SPECTRUM(double PBAR_TOA_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_TOA[DIM_TAB_PBAR+1], double A_nuclei)
{
	long i_pbar;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA;
	double flux_pbar_TOA;

	FILE* results;
	if(A_nuclei == 1){
		results = fopen(pbar_TOA_spectrum_file_name,"w");
		printf("Printing in file: %s\n", pbar_TOA_spectrum_file_name);
	}
	else if (A_nuclei == 2){
		results = fopen(Debar_TOA_spectrum_file_name,"w");
		printf("Printing in file: %s\n", Debar_TOA_spectrum_file_name);
	}
	else if (A_nuclei == 3){
		results = fopen(He3bar_TOA_spectrum_file_name,"w");
		printf("Printing in file: %s\n", He3bar_TOA_spectrum_file_name);
	}
			else if (A_nuclei == 4){
		results = fopen(He4bar_TOA_spectrum_file_name,"w");
		printf("Printing in file: %s\n", He4bar_TOA_spectrum_file_name);
	}
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_TOA = T_PBAR_TOA[i_pbar];
		flux_pbar_TOA = PBAR_TOA_SPECTRUM[i_pbar];

		fprintf(results, " %.10e\t %.10e\t \n", T_pbar_TOA, (1.0e04*flux_pbar_TOA));									// [#pbar m^{-2} sr^{-1} s^{-1} GeV/n^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

void print_PROTON_IS_SPECTRUM(double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1])
{
	long i_proton;
	double T_proton ,E_proton ,flux_proton_IS;
	double flux_proton;

	FILE* results;
	results = fopen(proton_IS_spectrum_file_name,"w");

	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
		T_proton = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
		E_proton = T_proton + MASSE_PROTON;

		flux_proton_IS = PROTON_IS_SPECTRUM[i_proton];

		fprintf(results, " %.10e\t %.10e\t \n", T_proton, (1.0e04*flux_proton_IS));											// [#pbar m^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

void print_PROTON_TOA_SPECTRUM(double PROTON_TOA_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], double T_PROTON_TOA[DIM_TAB_PROTON_SPECTRUM+1])
{
	long i_proton;
	double T_proton_TOA ,E_proton_TOA ,flux_proton_TOA;

	FILE* results;
	results = fopen(proton_TOA_spectrum_file_name,"w");

	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
		T_proton_TOA = T_PROTON_TOA[i_proton];

		flux_proton_TOA = PROTON_TOA_SPECTRUM[i_proton];

		fprintf(results, " %.10e\t %.10e\t \n", T_proton_TOA, (1.0e04*flux_proton_TOA));											// [#pbar m^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

void print_PROTON_SPECTRUM_exp(void)
{
	long i_proton;
	double T_proton ,E_proton ,flux_proton_exp;
	double flux_proton;

	FILE* results;
	results = fopen(proton_exp_spectrum_file_name,"w");

	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
		T_proton = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
		E_proton = T_proton + MASSE_PROTON;

		flux_proton_exp = flux_proton_EXP(E_proton);

		fprintf(results, " %.10e\t %.10e\t \n", T_proton, (1.0e04*flux_proton_exp));											// [#proton m^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

void print_HELIUM_SPECTRUM_exp(void)
{
	long i_proton;
	double T_proton ,E_proton ,flux_proton_exp;
	double flux_proton;

	FILE* results;
	results = fopen(helium_exp_spectrum_file_name,"w");

	for (i_proton=0;i_proton<=DIM_TAB_PROTON_SPECTRUM;i_proton++)
	{
		T_proton = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
		E_proton = T_proton + MASSE_PROTON;

		flux_proton_exp = flux_helium_EXP(E_proton);

		fprintf(results, " %.10e\t %.10e\t \n", T_proton, (1.0e04*flux_proton_exp));											// [#pbar m^{-3} sr^{-1} s^{-1} GeV^{-1}]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

//	On imprime le resultat

void print_PBAR_OVER_P_IS_SPECTRUM(double PBAR_OVER_P_IS_SPECTRUM[DIM_TAB_PBAR+1])
{
	long i_pbar;
	double T_IS ,E_IS , pbar_over_p_IS;

	FILE* results;
	results = fopen(pbar_over_p_IS_spectrum_file_name,"w");

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_IS = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_IS = T_IS + MASSE_PROTON;

		pbar_over_p_IS = PBAR_OVER_P_IS_SPECTRUM[i_pbar];																					// [NO UNIT]

		fprintf(results, " %.10e\t %.10e\t \n", T_IS, pbar_over_p_IS);
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

//	On imprime le resultat

void print_PBAR_OVER_P_TOA_SPECTRUM(double PBAR_OVER_P_TOA_SPECTRUM[DIM_TAB_PBAR+1], double T_PBAR_OVER_P_TOA[DIM_TAB_PBAR+1])
{
	long i_pbar;
	double T_pbar_over_p_TOA,pbar_over_p_TOA;

	FILE* results;
	results = fopen(pbar_over_p_TOA_spectrum_file_name,"w");

	//for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	for (i_pbar=DIM_TAB_PBAR;i_pbar>=0;i_pbar--)
	{
		T_pbar_over_p_TOA = T_PBAR_OVER_P_TOA[i_pbar];
		pbar_over_p_TOA = PBAR_OVER_P_TOA_SPECTRUM[i_pbar];

		fprintf(results, " %.10e\t %.10e\t \n",T_pbar_over_p_TOA,pbar_over_p_TOA);																					// [NO UNIT]
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

//	On imprime le resultat

void print_PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY(double PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2])
{
	long i_pbar;
	double T_IS ,E_IS , pbar_over_p_IS_INF, pbar_over_p_IS_SUP;

	FILE* results;
	results = fopen(pbar_over_p_IS_uncertainty_spectrum_file_name,"w");

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_IS = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_IS = T_IS + MASSE_PROTON;

		pbar_over_p_IS_INF = PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY[i_pbar][0];																					// [NO UNIT]
		pbar_over_p_IS_SUP = PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY[i_pbar][1];																					// [NO UNIT]

		fprintf(results, " %.10e\t %.10e\t %.10e\t \n", T_IS, pbar_over_p_IS_INF, pbar_over_p_IS_SUP);
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/


//	On imprime le resultat

void print_PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY(double PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2], double T_PBAR_OVER_P_TOA[DIM_TAB_PBAR+1])
{
	long i_pbar;
	double T_TOA ,E_TOA , pbar_over_p_TOA_INF, pbar_over_p_TOA_SUP;

	FILE* results;
	results = fopen(pbar_over_p_TOA_uncertainty_spectrum_file_name,"w");

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_TOA = T_PBAR_OVER_P_TOA[i_pbar];

		pbar_over_p_TOA_INF = PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY[i_pbar][0];																					// [NO UNIT]
		pbar_over_p_TOA_SUP = PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY[i_pbar][1];																					// [NO UNIT]

		fprintf(results, " %.10e\t %.10e\t %.10e\t \n", T_TOA, pbar_over_p_TOA_INF, pbar_over_p_TOA_SUP);
	}

	fclose(results);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/


void print_secondary_source_term(double PRIMARY_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], struct Structure_Cross_Section* pt_Cross_Section, double TARGET_DENSITY,  double A_nuclei){

		double T_proton ,E_proton,T_proton_previous ,E_proton_previous,T_pbar,flux_proton_IS, flux_helium_IS,CS_H_ON_H,CS_H_ON_He,CS_He_ON_H,CS_He_ON_He, resultat;

		FILE* results;
		results = fopen(secondary_source_term_file_name,"w");
		double weight_SIMSPON[DIM_TAB_PROTON+1];
		double E,dE,E_min,E_max,dlog_E_proton;
		int i_proton, i_pbar;



		dlog_E_proton = pow((E_PROTON_MAX/E_PROTON_MIN),(1./(double)DIM_TAB_PROTON)) - 1.;
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		resultat = 0;
		for (i_proton=0;i_proton<=DIM_TAB_PROTON;i_proton++)
		{
						E_proton = E_PROTON_MIN *
						pow((E_PROTON_MAX/E_PROTON_MIN),((double)i_proton/(double)DIM_TAB_PROTON));


						if (i_proton==0 || i_proton==DIM_TAB_PROTON) {weight_SIMSPON[i_proton] = 1./3.;}
						else {weight_SIMSPON[i_proton] = (1. + (double)(i_proton % 2)) * 2. / 3.;}
						// T_proton = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),((double)i_proton/(double)DIM_TAB_PROTON_SPECTRUM));
						// T_proton_previous = T_PROTON_SPECTRUM_MIN * pow((T_PROTON_SPECTRUM_MAX/T_PROTON_SPECTRUM_MIN),(((double) i_proton-1)/(double)DIM_TAB_PROTON_SPECTRUM));
						// E_proton = T_proton + MASSE_PROTON;
						// E_proton_previous = T_proton_previous + MASSE_PROTON;
						// dE = E_proton - E_proton_previous;
						// printf("dE %e E_proton %E\n", dE,E_proton_previous);
						flux_proton_IS = flux_proton_EXP(E_proton);
						flux_helium_IS = flux_helium_EXP(E_proton);
						CS_H_ON_H = pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton];
						CS_H_ON_He = pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]*pow(4,2.2/3);
						CS_He_ON_H = pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]*pow(4,2.2/3);
						CS_He_ON_He = pt_Cross_Section->DSPBAR_SUR_DEPBAR_H_ON_H[i_pbar][i_proton]*pow(16,2.2/3);
						resultat += 4*PI*(flux_proton_IS*DENSITE_H_DISC*CS_H_ON_H+flux_proton_IS*DENSITE_HE_DISC*CS_H_ON_He+flux_helium_IS*CS_He_ON_H*DENSITE_H_DISC+flux_helium_IS*CS_He_ON_He*DENSITE_HE_DISC)
							*weight_SIMSPON[i_proton] * dlog_E_proton;
						//* E_proton;
				}
				T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR)); //T/n
				// fprintf(stderr, " %.10e\t %.10e\t \n", T_pbar, (1.0e06*resultat));											// [#pbar m^{-3} s^{-1} GeV^{-1}]
				fprintf(results, " %.10e\t %.10e\t \n", T_pbar, (1.0e06*resultat));											// [#pbar m^{-3} s^{-1} GeV^{-1}]
		}
		printf("printing in file: %s\n", secondary_source_term_file_name);

}
