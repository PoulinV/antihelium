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
	double PBAR_IS_SPECTRUM[DIM_TAB_PBAR+1];
	double T_TOA[DIM_TAB_PBAR+1];
	double PBAR_TOA_SPECTRUM[DIM_TAB_PBAR+1];
	
	
	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA,flux_proton_TOA;
	double flux_pbar, flux_pbar_TOA;


//	INITALISATION DES VARIABLES
	///////////////////////////
	
	probleme = fopen("PB_SUMMARY","w");
	
	
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
	PBAR_SPECTRUM_initialization(PBAR_IS_SPECTRUM);
	PBAR_SPECTRUM_initialization(PBAR_TOA_SPECTRUM);

//	CALCUL DES Pi DES PROTONS ET D'HELIUM
/////////////////////////////////////////
	
	calculation_BESSEL_PROTON_Ep_i(alpha_i, &Proton, &Propagation);
	calculation_BESSEL_HELIUM_Ep_i(alpha_i, &Helium, &Propagation);
	

//	CALCUL DES FLUX DE PROTONS ET D'HELIUM
//////////////////////////////////////////
	

	
//	CALCUL DES Pi DES ANTIPROTONS
/////////////////////////////////

//	CALCUL DE LA CONTRIBUTION PRIMAIRE PROVENANT DE L'ANNIHILATION DES NEUTRALINOS.
	calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,500, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);
	//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);

//	CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.
	calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation);
	
//	CALCUL DES Pi EN PRENANT EN COMPTE LES PERTES D'ENERGIE ET LA REACCELERATION DIFFUSIVE.
	//tertiary_component_effect_calculation(&Pbar, alpha_i);
	ELDR_effect_calculation(&Propagation, &Pbar, alpha_i);


//	CALCUL DU FLUX D'ANTIPROTONS
////////////////////////////////

	PBAR_IS_SPECTRUM_calculation(PBAR_IS_SPECTRUM, &Pbar, &Propagation, alpha_i);	
	PBAR_TOA_SPECTRUM_calculation(PBAR_IS_SPECTRUM, PBAR_TOA_SPECTRUM, T_TOA, &Propagation);


//	AFFICHAGE DU FLUX D'ANTIPROTONS
///////////////////////////////////
	
	print_PBAR_IS_SPECTRUM(PBAR_IS_SPECTRUM);
	print_PBAR_TOA_SPECTRUM(PBAR_TOA_SPECTRUM, T_TOA);
	//print_total_pbar_spectra_MIN_MED_MAX(&Proton, &Helium, &Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	

//	CALCUL DES SPECTRES D'ANTIPROTONS PRIMAIRES POUR LE PPPC DE M.Cirelli
	/////////////////////////////////////////////////////////////////////
	
	//primary_spectra_BCGS_2014(&Pbar, &Cross_Section, &Propagation, &Primary_Source_Term, alpha_i);
	

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


	








