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
	double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1], PBAR_TOA_SPECTRUM[DIM_TAB_PBAR+1], T_PBAR_TOA[DIM_TAB_PBAR+1];
	double PROTON_IS_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], PROTON_TOA_SPECTRUM[DIM_TAB_PROTON_SPECTRUM+1], T_PROTON_TOA[DIM_TAB_PROTON_SPECTRUM+1];
	double PBAR_OVER_P_IS_SPECTRUM[DIM_TAB_PBAR+1], PBAR_OVER_P_TOA_SPECTRUM[DIM_TAB_PBAR+1], T_PBAR_OVER_P_TOA[DIM_TAB_PBAR+1];
	double PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2], PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY[DIM_TAB_PBAR+1][2];

	double EnTOA, pnTOA, flux_TOA, TnTOA;

//	INITALISATION DES VARIABLES
	///////////////////////////
	
	probleme = fopen("PB_SUMMARY","w");
	
	
// 	CALCULS PRELIMINAIRES 
	/////////////////////
	
	TABLE_PROPAGATION_loading(&Propagation);
	MIN_MED_MAX_loading(&Propagation);
	print_propagation_parameters(&Propagation);
	//propagation_parameters_loading(&Propagation, 1);
	//print_propagation_parameters(&Propagation);
	
		
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
	
	PROTON_SPECTRUM_initialization(PROTON_IS_SPECTRUM);
	PROTON_SPECTRUM_initialization(PROTON_TOA_SPECTRUM);
	PBAR_BESSEL_TABLES_123_initialization(&Pbar);
	PBAR_SPECTRUM_initialization(PBAR_IS_SPECTRUM);
	PBAR_SPECTRUM_initialization(PBAR_TOA_SPECTRUM);
	PBAR_SPECTRUM_initialization(PBAR_OVER_P_IS_SPECTRUM);
	PBAR_SPECTRUM_initialization(PBAR_OVER_P_TOA_SPECTRUM);


//	CALCUL DU FLUX DE PROTONS
//////////////////////////////////////////
	
	PROTON_IS_SPECTRUM_calculation(PROTON_IS_SPECTRUM, &Proton, &Propagation, alpha_i);
	PROTON_TOA_SPECTRUM_calculation(PROTON_IS_SPECTRUM, PROTON_TOA_SPECTRUM, T_PROTON_TOA, &Propagation);
	
	
//	CALCUL DES Pi DES ANTIPROTONS
/////////////////////////////////

//	CALCUL DE LA CONTRIBUTION PRIMAIRE PROVENANT DE L'ANNIHILATION DES NEUTRALINOS.
	calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,500, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);
	//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);

//	CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.
	calculation_BESSEL_PROTON_Ep_i(alpha_i, &Proton, &Propagation);
	calculation_BESSEL_HELIUM_Ep_i(alpha_i, &Helium, &Propagation);
	//calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation);
	calculation_BESSEL_PBAR_SUM_123_Epbar_i(&Pbar);
	
	
//	CALCUL DES Pi EN PRENANT EN COMPTE LES PERTES D'ENERGIE ET LA REACCELERATION DIFFUSIVE.
	tertiary_component_effect_calculation(&Pbar, alpha_i);
	ELDR_effect_calculation(&Propagation, &Pbar, alpha_i);


//	CALCUL DU FLUX D'ANTIPROTONS
////////////////////////////////

	PBAR_IS_SPECTRUM_calculation(PBAR_IS_SPECTRUM, &Pbar, &Propagation, alpha_i);	
	//PBAR_TOA_SPECTRUM_calculation(PBAR_IS_SPECTRUM, PBAR_TOA_SPECTRUM, T_PBAR_TOA, &Propagation);

//	CALCUL DU RAPPORT Pbar/P
////////////////////////////
	
	PBAR_OVER_P_IS_SPECTRUM_calculation(PBAR_OVER_P_IS_SPECTRUM, &Proton, &Pbar, &Propagation, alpha_i);
	PBAR_OVER_P_TOA_SPECTRUM_calculation(PBAR_OVER_P_TOA_SPECTRUM, T_PBAR_OVER_P_TOA, &Proton, &Pbar, &Propagation, alpha_i);
	
	//PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY_calculation(PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation, alpha_i);
	//PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY_calculation_1(PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY, PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY, T_PBAR_OVER_P_TOA, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation, alpha_i);
		

//	AFFICHAGE DES SPECTRES
//////////////////////////
	
	print_PROTON_IS_SPECTRUM(PROTON_IS_SPECTRUM);
	print_PROTON_TOA_SPECTRUM(PROTON_TOA_SPECTRUM, T_PROTON_TOA);
	print_PROTON_SPECTRUM_exp();
	print_HELIUM_SPECTRUM_exp();
	
	
	
	print_PBAR_IS_SPECTRUM(PBAR_IS_SPECTRUM);
	//print_PBAR_TOA_SPECTRUM(PBAR_TOA_SPECTRUM, T_PBAR_TOA);
	//print_total_pbar_spectra_MIN_MED_MAX(&Proton, &Helium, &Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	
	//print_PBAR_OVER_P_IS_SPECTRUM(PBAR_OVER_P_IS_SPECTRUM);
	print_PBAR_OVER_P_TOA_SPECTRUM(PBAR_OVER_P_TOA_SPECTRUM, T_PBAR_OVER_P_TOA);
	
	//print_PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY(PBAR_OVER_P_IS_SPECTRUM_UNCERTAINTY);
	//print_PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY(PBAR_OVER_P_TOA_SPECTRUM_UNCERTAINTY, T_PBAR_OVER_P_TOA);
		
	

//	CALCUL DES SPECTRES D'ANTIPROTONS PRIMAIRES POUR LE PPPC DE M.Cirelli
	/////////////////////////////////////////////////////////////////////
	
	
	primary_spectra_BCGS_2014(&Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	
	
	
//	TEST
	
	TnTOA = 100.0;
	EnTOA = TnTOA + MASSE_PROTON;
    pnTOA = sqrt(pow(EnTOA,2) - pow(MASSE_PROTON,2));
	
	flux_TOA = fit_proton_flux_AMS02_yoann(pnTOA);
	flux_TOA *= EnTOA / pnTOA;
	
	
	
	printf("TnTOA = %.5e \t flux_TOA = %.5e  \n", EnTOA, flux_TOA);
		
	
	
	

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



