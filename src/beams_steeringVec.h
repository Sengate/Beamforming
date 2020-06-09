#inclued "steering_vector.h"
#include "convert_hadec_to_enu.h"
#include "convert_theta_phi_to_enu.h"
#include "convert_enu_to_azalt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include <math.h>


/**
 *@brief Structure containing TAB beams positions (ras and decs) and observor's latitude and local sidereal time(lst).
 **/
typedef struct{
    float nbeams;
    float *ras;
    float *decs;
}beams_info;



/***-------------------------------------------------------------------------------------------------------*/
void make3Dmatrix_complex(unsigned dim1, unsigned dim2, unsigned dim3, complex (**aptr)[dim2][dim3]) {
    *aptr = malloc(dim1 * sizeof(**aptr));
    
    if (!aptr) {printf("Memo not allocated");}
    else
    {printf("sure");}
}

float *ha;
float *beams_azimuths;
float *beams_altitudes;



//Allocate memories
void beams_phases_mem(Data_parameters *dat, beams_info *TABs){
    
    complex (*EW_beams_stv)[dat->nchan][dat->nants];
    complex (*NS_beams_stv)[dat->nchan][dat->nants];
    
    make3Dmatrix_complex(TABs->nbeams, dat->nchan, dat->nants, &EW_beams_stv);
    make3Dmatrix_complex(TABs->nbeams, dat->nchan, dat->nants, &NS_beams_stv);
    
    ha = (float*)malloc(sizeof(float) * TABs->nbeams);
    beams_azimuths = (float*)malloc(sizeof(float) * TABs->nbeams);
    beams_altitudes = (float*)malloc(sizeof(float) * TABs->nbeams);
    
}

/*void beams_phases_destroy(){
 free(ha);
 free(beams_altitudes);
 free(beams_altitudes);
 free(EW_beams_stv);
 free(NS_beams_stv);
 }
 */





void get_TAB_phases(Data_parameters *dat, beams_info *beams_param, float lst, float obs_latitude){
    
    //Allocate memories
    float *ha = (float*)malloc(sizeof(float) * beams_param->nbeams);
    float *beams_azimuths = (float*)malloc(sizeof(float) * beams_param->nbeams);
    float *beams_altitudes = (float*)malloc(sizeof(float) * beams_param->nbeams);
    
    for (int ibeam=0; ibeam<beams_param->nbeams; ibeam++){
        ha[ibeam] = lst - beams_param->ras[ibeam];
        
        hadec equ = {ha[ibeam], beams_param->decs[ibeam]};
        //convert hadec to enu coordinates.
        equ2enu *hadec2enu = hadec_2_enu(equ, obs_latitude);
        
        //convert enu to altaz.
        altaz *enu2altaz = hor2Altaz(hadec2enu);
        
        beams_azimuths[ibeam] = enu2altaz->az;
        beams_altitudes[ibeam] = enu2altaz->alt;
        
        for (unsigned int ifreq =0; ifreq < dat->nchan; ifreq++){
            
            for (unsigned int iant = 0; iant < dat->nants; iant++){
                
                EW_beams_stv[ibeam][ifreq][iant].im = EWphases[ifreq][iant].im
                EW_beams_stv[ibeam][ifreq][iant].re = EWphases[ifreq][iant].re
                
                NS_beams_stv[ibeam][ifreq][iant].im = NSphases[ifreq][iant].im
                NS_beams_stv[ibeam][ifreq][iant].re = NSphases[ifreq][iant].re
                
            }
        }
    }
    
    free(ha);
    free(beams_azimuths);
    free(beams_altitudes);
}


    

