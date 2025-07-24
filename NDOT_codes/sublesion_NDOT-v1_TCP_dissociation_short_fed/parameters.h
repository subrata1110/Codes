#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <string>
#include <cmath>
using namespace std ;


// Define the pararameters neccesarry for the Simulation
// Define the pararameters neccesarry for the Simulation

const int N= 100; // number of spatial grid points in Plaque
const int NT=200; // number of spatial grid points in teeth
const int NF=3; // number of spatial points in film saliva
const double L= 250.0*pow(10,-6.0) ; // 100 um the size of the Plaque in m
const double LT=100.0*pow(10,-6.0) ; // 100 um the size of the teeth in m
const double LF=100.0*pow(10,-6.0) ;// 100 um the size of the film in m

const double dx=L/(N-1); // the grid length in Plaque in m
const double dxT=LT/(NT-1); // the grid length in the enamel in m
const double dxF=LF/(NF-1); //the grid length in the film in m
const double simu_time=5400*5; // total simulation time in sec
const int simu_step=pow(10,5); // equillibration time step for ions
const double dt=0.0001,dte=0.00001; // the time steps in simulation in sec
const double t_hf=0.5*60; // the residence time in seconds in film
const double t_hs=2.0*60; // the residence time in seconds in bulk 

const double a_vf=pow(10,4); // the Area/ Volume of film in m^-1
const double e=1.60217663*pow(10,-19.0); // Charge of electron in columbs
const double kb=1.380649*pow(10,-23.0); // boltzman constant in m^2 kg s^-2 K^-1
const double T=310; // Temperature in kelvin
const double R=8.3145; // Gas constant in J mol^-1 K^-1
const double F= 96485.0; //Faraday constant in Columb/mol
const double K_HAP_enamel=5.5*pow(10,-28.0); // enamel solubility constant mol^9/m3^9
const double nd=0.3; // Deminerlisation reaction order
const double nr=1.25; // Reminerlisation reaction order
const double md=2.8; // Deminerlisation reaction order
const double kd=0.42*pow(10,-3.0) /60.0; // Deminerlisation rate constant in mol^0.7m^-1.1s^-1
const double kr=7.44*1.0e-3; // Reminerlisation rate constant in mol^0.25m^0.75s^-1
const double epsilon=80*8.854*pow(10,-12.0) ; // Absolute dielectric constant of water in F/m or C/V.m
const int glu_status=1; // Glucose pulse is on or off
const int freq_ion_equilli=1; // frequency for ions equillibration
const int n_fixed_ions=0; // number of fixed ions in plaque
const int nr_fixed=0; // number of reaction which are for fixed charges
const double Ke_h2o=pow(10,-8.0); // Water dissociation constant
const double Ke_pho2n=pow(10.0,-9.32); // dissociation constant for pho--
// Max rd_HAP 1.2*10^-6 mol/m^2/s

const double xf= 0.0; // Degree of substitution of Fluoride ion
const double K_HAP_enamel_F= 7.41*pow(10.0,-33.0); //solubility constant for FHAP surface mol^9/m^27
const double t_rest=100.0, t_step=10.0, t_feed=120.0 ;
const double t_gap=5270.0;
const double t_cycle=t_step+t_feed+t_gap;
const int ncycle=15;
const double glu_peak=560.0;
const double t_total=t_cycle*ncycle;
const double tcp_sp_area=0.378e10; //specific area of the TCP /m 
 const double K_sp_tcp=3.0e-15; //solubiity product if TCP
//const double K_sp_tcp=2.51e-15; //solubiity product if TCP
// const double K_sp_tcp=exp(10.18); //solubiity product if TCP

const double r_rod=2.5*pow(10.0,-8.0); //radius of the HAP cystal in m
const double M_HAP=0.502 ; // HAP molar weight kg/mol
const double rho_HAP=3160.0; //density of HAP kg/m^3 
// Class for ions in the system
const string output_filename="all_ion_tim_tcp_1wt.csv";
const string output_remin_demin="Remin_demin_tcp_1wt.csv";


#endif
