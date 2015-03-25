#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include <sys/types.h> 
#include <sys/stat.h> 

#include "COMMON.h"
#include "./sources/STRUCTURES.h"

#include "./sources/BESSEL_PRELIM.h"
#include "./sources/CROSS_SECTIONS.h"
#include "./sources/DIFFUSION_PROPAGATION.h"
#include "./sources/SOLAR_MOD.h"
#include "./sources/PROTON.h"
#include "./sources/HELIUM.h"
#include "./sources/PRIMARY_PBAR.h"
#include "./sources/ANTI_PROTON.h"
#include "./sources/spectra.h"


/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

// DECLARATIONS DES FONCTIONS TEMPORAIRES

void PBAR_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PBAR+1]);
void print_propagation_parameters(struct Structure_Propagation* pt_Propagation);
void PBAR_BESSEL_TABLES_123_initialization(struct Structure_Pbar* pt_Pbar);
void PBAR_IS_SPECTRUM_calculation(double PBAR_SPECTRUM[DIM_TAB_PBAR+1], struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1]);

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

FILE *probleme;

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

int main(void)
{
	time_t t1,t2;
	double duree;
	t1 = time(NULL);
	printf("\n");


//	DECLARATION DES VARIABLES
	/////////////////////////	
	
	struct Structure_Nuclei Proton;
	struct Structure_Nuclei Helium;
	struct Structure_Pbar Pbar;
	struct Structure_Cross_Section Cross_Section;
	struct Structure_Propagation Propagation;
	struct Structure_Primary_Source_Term Primary_Source_Term;	
	
	long   i_data,i_iteration,i_pbar,i;
	
	double alpha_i[NDIM+1];
	double PBAR_SPECTRUM[DIM_TAB_PBAR+1];
	double T_TOA[DIM_TAB_PBAR+1];
	double PBAR_SPECTRUM_TOA[DIM_TAB_PBAR+1];
	
	
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA,flux_proton_TOA;
	double flux_pbar, flux_pbar_TOA;


//	INITALISATION DES VARIABLES
	///////////////////////////
	
	FILE* results;
	
	results = NULL;

	probleme = fopen("PB_SUMMARY","w");
	results = fopen(pbar_IS_spectrum_file_name,"w");
	
	
// 	CALCULS PRELIMINAIRES 
	/////////////////////
	
	MIN_MED_MAX_loading(&Propagation);
	print_propagation_parameters(&Propagation);
		
	bessel_preliminary_write_file(alpha_i, &Proton, &Helium);
	bessel_preliminary_read_file (alpha_i, &Proton, &Helium);

	CLEANING_ALL_THE_DSPBAR_SUR_DEPBAR   (&Cross_Section);
	//DSPBAR_SUR_DEPBAR_H_ON_H_write_file  (&Cross_Section);
	//DSPBAR_SUR_DEPBAR_H_ON_HE_write_file (&Cross_Section);
	//DSPBAR_SUR_DEPBAR_HE_ON_H_write_file (&Cross_Section);
	//DSPBAR_SUR_DEPBAR_HE_ON_HE_write_file(&Cross_Section);
	DSPBAR_SUR_DEPBAR_H_ON_H_read_file   (&Cross_Section);
	DSPBAR_SUR_DEPBAR_H_ON_HE_read_file  (&Cross_Section);
	DSPBAR_SUR_DEPBAR_HE_ON_H_read_file  (&Cross_Section);
	DSPBAR_SUR_DEPBAR_HE_ON_HE_read_file (&Cross_Section);
	
	DM_source_term_calculation(&Primary_Source_Term);
	
	PBAR_BESSEL_TABLES_123_initialization(&Pbar);
	
	
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	

	PBAR_SPECTRUM_initialization(PBAR_SPECTRUM);
	PBAR_SPECTRUM_initialization(PBAR_SPECTRUM_TOA);
	
	
	
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	ON PEUT Y ALLER !
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	

			

	

//	CALCUL DE LA CONTRIBUTION PRIMAIRE PROVENANT DE L'ANNIHILATION DES NEUTRALINOS.
	calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,500, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);
	//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);

//	CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.
	calculation_BESSEL_PROTON_Ep_i(alpha_i, &Proton, &Propagation);
	calculation_BESSEL_HELIUM_Ep_i(alpha_i, &Helium, &Propagation);
	calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation);
	calculation_BESSEL_PBAR_SUM_123_Epbar_i(&Pbar);

//	CALCUL DU SPECTRE FINAL DES ANTIPROTONS.

	//goto TEST;
	for (i_iteration=1;i_iteration<=5;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, &Pbar);
		calculation_BESSEL_PBAR_SUM_123_Epbar_i(&Pbar);
	}

	//goto TEST;
	for (i_iteration=1;i_iteration<=5;i_iteration++)
	{
		calculation_BESSEL_PBAR_TERTIARY_Epbar_i(alpha_i, &Pbar);
		calculation_BESSEL_PBAR_TOT_direct_inversion_A(&Pbar, &Propagation);
		//calculation_BESSEL_PBAR_TOT_direct_inversion_B(&Pbar, &Propagation);
		//calculation_BESSEL_PBAR_TOT_direct_inversion_GJ_NR(&Pbar, &Propagation);
		//calculation_BESSEL_PBAR_TOT_diffusion_soluce_A(15., 1200., &Pbar, &Propagation);
	}

//			We compute now the antiproton spectrum and store for each KINETIC
//			ENERGY the lowest -- PBAR_SPECTRUM_MIN -- the medium
//			-- PBAR_SPECTRUM_MED -- and the largest -- PBAR_SPECTRUM_MAX --
//			values which we meet.




TEST:
/*
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;
		for (i=1;i<=NDIM;i++)
		{
			Pbar.BESSEL_PBARi[i] = Pbar.BESSEL_PBAR_TOT_Epbar_i[i_pbar][i];
		}
		flux_antiproton_IS = GENERIC_FLUX_04(R_EARTH,0.,E_pbar_IS,MASSE_PROTON,1.,alpha_i,Pbar.BESSEL_PBARi, &Propagation);
		//flux_antiproton_IS = GENERIC_FLUX(R_EARTH,0.,E_pbar_IS,MASSE_PROTON,1.,alpha_i,Pbar.BESSEL_PBARi, &Propagation);

		
		PBAR_SPECTRUM[i_pbar] = flux_antiproton_IS;
		
	}
*/
	PBAR_IS_SPECTRUM_calculation(PBAR_SPECTRUM, &Pbar, &Propagation, alpha_i);



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//			On imprime le resultat


	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

		flux_pbar = PBAR_SPECTRUM[i_pbar];

		fprintf(results, " %.10e\t %.10e\t \n", T_pbar_IS, (1.0e04*flux_pbar));	
	}
	
	fclose(results);

	//goto LA_FIN;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Modulation des spectres IS ===> TOA et impression.

	Propagation.PHI_FISK = fisk_potential;
	
	results = fopen(pbar_TOA_spectrum_file_name,"w");

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

//		Nous modulons maintenant les spectres PBAR obtenus.

		if (E_pbar_TOA <= MASSE_PROTON)
		{
			T_TOA[i_pbar]             = 0.0;
			PBAR_SPECTRUM_TOA[i_pbar] = 0.0;
			//continue;
		}

		FFA_IS_to_TOA(1.,1.,Propagation.PHI_FISK,E_pbar_IS,PBAR_SPECTRUM[i_pbar],&E_pbar_TOA,&flux_pbar_TOA);

//		Nous les imprimons !

		T_pbar_TOA = E_pbar_TOA - MASSE_PROTON;
		fprintf(results, " %.10e\t %.10e\t \n", T_pbar_TOA, (1.0e04*flux_pbar_TOA));

//		Nous les stockons en memoire dans les tableaux RESULTS_T_TOA[DIM_TAB_PBAR+1] et RESULTS_SPECTRUM_TOA_MIN_MED_MAX[DIM_TAB_PBAR+1];
	
		T_TOA[i_pbar] = T_pbar_TOA;
		PBAR_SPECTRUM_TOA[i_pbar] = (1.0e04*flux_pbar_TOA);
		
	}
	fclose(results);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	//primary_spectra_BCGS_2014(&Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	//print_total_pbar_spectra_MIN_MED_MAX(&Proton, &Helium, &Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


LA_FIN :

	t2 = time(NULL);
	duree = difftime(t2,t1);
	printf("\n\n");
	printf(" duree = %f \n\n",duree);
	fprintf(probleme,"\n\n");
	fprintf(probleme," duree = %f \n\n",duree);
	fclose(probleme);
	return EXIT_SUCCESS;
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

// FONCTIONS TEMPORAIRES

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/


void PBAR_SPECTRUM_initialization(double SPECTRUM[DIM_TAB_PBAR+1])
{
	long i_pbar;
	
	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
	    SPECTRUM[i_pbar] = 0.0;
	}
	
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

//	Nous imprimons les coefficients de diffusion_propagation choisis dans le calcul.

void print_propagation_parameters(struct Structure_Propagation* pt_Propagation)
{
	printf(" DELTA           = %.5e [NO UNIT]\n",pt_Propagation->PUISSANCE_COEFF_DIFF);
	printf(" DIFFUSION_0_GV  = %.5e [cm^{2} s^{-1}]\n",pt_Propagation->DIFFUSION_0_GV);
	printf(" E_DIFFUS        = %.5e [kpc]\n",pt_Propagation->E_DIFFUS);
	printf(" VENT_GALACTIQUE = %.5e [cm s^{-1}]\n",pt_Propagation->VENT_GALACTIQUE);
	printf(" V_ALFEN         = %.5e [cm s^{-1}]\n\n",pt_Propagation->V_ALFEN);
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

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

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/

void PBAR_IS_SPECTRUM_calculation(double PBAR_SPECTRUM[DIM_TAB_PBAR+1], struct Structure_Pbar* pt_Pbar, struct Structure_Propagation* pt_Propagation, double alpha_i[NDIM+1])
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

		
		PBAR_SPECTRUM[i_pbar] = flux_antiproton_IS;	
	}
}

/****************************************************************************************************************************************************************************************/
/****************************************************************************************************************************************************************************************/





























