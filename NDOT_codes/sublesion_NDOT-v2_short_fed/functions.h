#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "class.h"
#include "parameters.h"
// #include "parameters.cpp"

#include <string>
#include <cmath>
using namespace std ;

void reaction_plaque(Ion* iobj, Reaction_ion_bio* robj,   double conc[][N+NT-1], int ni,
 int nrr, int x, string inh_type[], double K_inh[], int n_inh);


void reaction_enamel(Ion* iobj, double conc_enamel[][N+NT-1],  
 int ni, int x, double E_por[] , double rd_HAP[],double rr_HAP[], double rd_ca[], double rr_ca[], double DSe[]);

void diffusion(Ion* iobj, double conc[][N+NT-1], double reaction[], int ni, int x, double E_por[]);

void electro_diffusion(Ion* iobj, double conc[][N+NT-1], double reaction[], 
 double phi[N],  int ni, int x, double E_por[]);

void equillibriate_ions_film( Ion* iobj, Reaction_ion_equi* robj, double conc[], int ni,int nr );


void equillibriate_ions(Ion* iobj, Reaction_ion_equi* robj, double conc[][N+NT-1], int ni,int nr, int x);

void boundary_update(Ion* iobj, double conc_film[], double conc[][N+NT-1],  
int ni, double E_por[]);

void phi_update(Ion* iobj, double conc[][N+NT-1],double phi[N], double new_phi[N], int ni);

int id_ion(Ion* iobj, string ion_name, int ni);

void exchange_bulk_film(double conc_bulk[], double conc_film[], double reaction[],int size);

void exchange_film_plaque(Ion* iobj, double conc[][N+NT-1], double conc_film[], 
double reaction[], int size);

void film_conc_update(double conc_film[], double rate1[],double rate2[],int ni);

void reaction(Ion* iobj, Reaction_ion_bio* robj, double conc[][N+NT-1],
int ni, int nrr, int x,  string inh_type[],  double K_inh[], int n_inh,double E_por[],double rd_HAP[],double rr_HAP[],
 double rd_ca[], double rr_ca[], double DSe[]);

double glu_pulse(double current_time, double glu_conc);







#endif