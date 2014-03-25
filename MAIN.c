#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <locale.h>

#include "COMMON.h"
#include "STRUCTURES.h"

#include "BESSEL_PRELIM.h"
#include "CROSS_SECTIONS.h"
#include "DIFFUSION_PROPAGATION.h"
#include "SOLAR_MOD.h"
#include "PROTON.h"
#include "HELIUM.h"
#include "PRIMARY_PBAR.h"
#include "ANTI_PROTON.h"

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

FILE *probleme;

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/

int main(void)
{
	time_t t1,t2;
	double duree;

	long   i_data,i_iteration,i_pbar,i;

	double T_pbar_IS ,E_pbar_IS ,flux_antiproton_IS ,flux_proton_IS;
	double T_pbar_TOA,E_pbar_TOA,flux_antiproton_TOA,flux_proton_TOA;

	double FLUX_PBAR_MIN,    FLUX_PBAR_MED,    FLUX_PBAR_MAX;
	double FLUX_PBAR_TOA_MIN,FLUX_PBAR_TOA_MED,FLUX_PBAR_TOA_MAX;

	static double PBAR_SPECTRUM_MIN       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MED       [DIM_TAB_PBAR+1];
	static double PBAR_SPECTRUM_MAX       [DIM_TAB_PBAR+1];
	static double RESULTS_T_TOA           [DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MIN[DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MED[DIM_TAB_PBAR+1];
	static double RESULTS_SPECTRUM_TOA_MAX[DIM_TAB_PBAR+1];

	double alpha_i[NDIM+1];

	struct Structure_Nuclei              Proton;
	struct Structure_Nuclei              Helium;
	struct Structure_Pbar                Pbar;
	struct Structure_Cross_Section       Cross_Section;
	struct Structure_Propagation         Propagation;
	struct Structure_Primary_Source_Term Primary_Source_Term;

	int    channel;
	double mass_chi;

	FILE *results;

	t1 = time(NULL);
	probleme = fopen("PB_SUMMARY","w");

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

	mass_chi = mass_chi_choice;
	channel  = channel_choice;
	
	printf("\n"
	"mass_chi = %.2e GeV\n"
	"channel  = %d\n\n", mass_chi, channel);
	
	gaelle_preliminary					(&Primary_Source_Term);
	DNPBAR_ON_DTPBAR_gaelle_read_file   (mass_chi, &Primary_Source_Term);
	dNpbar_on_dEpbar_primary_calculation(mass_chi, channel, &Primary_Source_Term);
	primary_source_calculation          (mass_chi, &Primary_Source_Term);
	
///////////////////////////////////////////////////////////////////////////////////////////
//	ON PEUT Y ALLER !
///////////////////////////////////////////////////////////////////////////////////////////

	//for (i_data=1;i_data<=3;i_data++)
	for (i_data=2;i_data<=2;i_data++)
	{
//		Nous definissons a ce niveau les parametres que FIORENZA, DAVID et RICHARD
//		-- hereafter called FDR -- ont determines.

		if (i_data == 1)       // CAS MAX 
		{
			Propagation.PUISSANCE_COEFF_DIFF = 0.46;
			Propagation.DIFFUSION_0_GV  = 0.0765 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
			Propagation.E_DIFFUS        = 15.0;                                       		// [kpc]
			Propagation.VENT_GALACTIQUE = (5.0   * 1.0e5);                            		// [cm s^{-1}]
			Propagation.V_ALFEN         = (117.6 * 1.0e5);                            		// [cm s^{-1}]
		}
		else if (i_data == 2) // CAS MED 
		{
			Propagation.PUISSANCE_COEFF_DIFF = 0.70;
			Propagation.DIFFUSION_0_GV  = 0.0112 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
			Propagation.E_DIFFUS        = 4.0;                                        		// [kpc]
			Propagation.VENT_GALACTIQUE = (12.0  * 1.0e5);                            		// [cm s^{-1}]
			Propagation.V_ALFEN         = (52.9  * 1.0e5);                            		// [cm s^{-1}]
		}
		else if (i_data == 3) // CAS MIN 
		{
			Propagation.PUISSANCE_COEFF_DIFF = 0.85;
			Propagation.DIFFUSION_0_GV  = 0.0016 * pow(CM_PAR_KPC,2.) / SEC_PAR_MGYR; 		// [cm^{2} s^{-1}]
			Propagation.E_DIFFUS        = 1.0;                                        		// [kpc]
			Propagation.VENT_GALACTIQUE = (13.5  * 1.0e5);                            		// [cm s^{-1}]
			Propagation.V_ALFEN         = (22.4  * 1.0e5);                            		// [cm s^{-1}]
		}

//		Nous imprimons les coefficients de diffusion_propagation choisis dans le calcul.
		
		printf(" CAS NUMERO      = %ld \n",i_data);
		printf(" DELTA           = %.5e [NO UNIT]\n",Propagation.PUISSANCE_COEFF_DIFF);
		printf(" DIFFUSION_0_GV  = %.5e [cm^{2} s^{-1}]\n",Propagation.DIFFUSION_0_GV);
		printf(" E_DIFFUS        = %.5e [kpc]\n",Propagation.E_DIFFUS);
		printf(" VENT_GALACTIQUE = %.5e [cm s^{-1}]\n",Propagation.VENT_GALACTIQUE);
		printf(" V_ALFEN         = %.5e [cm s^{-1}]\n\n",Propagation.V_ALFEN);
		
//		On remet a zero les tableaux Pbar.BESSEL_PBAR_SEC_Epbar_i et Pbar.BESSEL_PBAR_TER_Epbar_i.

		for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
		{
			for (i=0;i<=NDIM;i++)
			{
				Pbar.BESSEL_PBAR_PRI_Epbar_i[i_pbar][i] = 0.0;
				Pbar.BESSEL_PBAR_SEC_Epbar_i[i_pbar][i] = 0.0;
				Pbar.BESSEL_PBAR_TER_Epbar_i[i_pbar][i] = 0.0;
			}
		}

//		CALCUL DE LA CONTRIBUTION PRIMAIRE PROVENANT DE L'ANNIHILATION DES NEUTRALINOS.

		calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,500, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);
		//calculation_BESSEL_PBAR_PRIMARY_Epbar_i(100,1000, alpha_i, &Pbar, &Propagation, &Primary_Source_Term);

//		CALCUL DE LA CONTRIBUTION SECONDAIRE PROVENANT DE LA SPALLATION DU GAZ INTERSTELLAIRE
//		PAR LES PROTONS ET LES HELIONS DU RAYONNEMENT COSMIQUE.

		calculation_BESSEL_PROTON_Ep_i(alpha_i, &Proton, &Propagation);
		calculation_BESSEL_HELIUM_Ep_i(alpha_i, &Helium, &Propagation);
		calculation_BESSEL_PBAR_SECONDARY_Epbar_i(alpha_i, &Proton, &Helium, &Pbar, &Cross_Section, &Propagation);
		calculation_BESSEL_PBAR_SUM_123_Epbar_i(&Pbar);

//		CALCUL DU SPECTRE FINAL DES ANTIPROTONS.

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

//		We compute now the antiproton spectrum and store for each KINETIC
//		ENERGY the lowest -- PBAR_SPECTRUM_MIN -- the medium
//		-- PBAR_SPECTRUM_MED -- and the largest -- PBAR_SPECTRUM_MAX --
//		values which we meet.

TEST:

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

			if      (i_data == 1){PBAR_SPECTRUM_MAX[i_pbar] = flux_antiproton_IS;}
			else if (i_data == 2){PBAR_SPECTRUM_MED[i_pbar] = flux_antiproton_IS;}
			else if (i_data == 3){PBAR_SPECTRUM_MIN[i_pbar] = flux_antiproton_IS;}
		}
	}

////////////////////////////////////////////////////////////////////////////////////////////
//	On imprime le resultat

	results = fopen(antiproton_spectrum_file_name,"w");
		
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
	goto LA_FIN;

////////////////////////////////////////////////////////////////////////////////////////////	
//	Modulation des spectres IS ===> TOA et impression.

	Propagation.PHI_FISK = fisk_potential;
	
	results = fopen(antiproton_spectrum_file_name,"w");

	for (i_pbar=0;i_pbar<=DIM_TAB_PBAR;i_pbar++)
	{
		T_pbar_IS = T_PBAR_MIN *
		pow((T_PBAR_MAX/T_PBAR_MIN),((double)i_pbar/(double)DIM_TAB_PBAR));
		E_pbar_IS = T_pbar_IS + MASSE_PROTON;

//		Nous modulons maintenant les spectres PBAR obtenus.

		FFA_IS_to_TOA(1.,1.,Propagation.PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MIN[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MIN);

		if (E_pbar_TOA <= MASSE_PROTON)
		{
			RESULTS_T_TOA[i_pbar]            = 0.0;
			RESULTS_SPECTRUM_TOA_MAX[i_pbar] = 0.0;
			RESULTS_SPECTRUM_TOA_MED[i_pbar] = 0.0;
			RESULTS_SPECTRUM_TOA_MIN[i_pbar] = 0.0;
			continue;
		}

		FFA_IS_to_TOA(1.,1.,Propagation.PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MED[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MED);
		FFA_IS_to_TOA(1.,1.,Propagation.PHI_FISK,E_pbar_IS,PBAR_SPECTRUM_MAX[i_pbar],&E_pbar_TOA,&FLUX_PBAR_TOA_MAX);

//		Nous les imprimons !

		T_pbar_TOA = E_pbar_TOA - MASSE_PROTON;
		fprintf(results, " %.10e\t %.10e\t %.10e\t %.10e\t \n", T_pbar_TOA, (1.0e04*FLUX_PBAR_TOA_MIN), (1.0e04*FLUX_PBAR_TOA_MED), (1.0e04*FLUX_PBAR_TOA_MAX));

//		Nous les stockons en memoire dans les tableaux RESULTS_T_TOA[DIM_TAB_PBAR+1] et RESULTS_SPECTRUM_TOA_MIN_MED_MAX[DIM_TAB_PBAR+1];
	
		RESULTS_T_TOA[i_pbar] = T_pbar_TOA;
		RESULTS_SPECTRUM_TOA_MAX[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MAX);
		RESULTS_SPECTRUM_TOA_MED[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MED);
		RESULTS_SPECTRUM_TOA_MIN[i_pbar] = (1.0e04*FLUX_PBAR_TOA_MIN);
	}
	fclose(results);

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

/**************************************************************************************************************************************************************************************************/
/**************************************************************************************************************************************************************************************************/
