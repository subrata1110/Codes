#include "functions.h"
#include "class.h"
#include "parameters.h"
#include <fstream>
#include<bits/stdc++.h>
#include<sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <chrono>
using namespace std;



void reaction_plaque(Ion* iobj, Reaction_ion_bio* robj,   double conc[][N+NT-1], int ni,
 int nrr, int x, string inh_type[], double K_inh[], int n_inh){
    double rate;
    int id_i;
    double reac_rate[ni];
    for(int i=0;i<ni;i++){
        // cout<<i<<endl;
        reac_rate[i]=0;
    }
    
    for (int i=0; i< nrr; i++){ // reaction id's
    
        rate=robj[i].q_reac*robj[i].bac_conc;
        for(int j=0; j<robj[i].nreactant;j++){ // ion ids for reactanats
            int ion_id= robj[i].reactants_id[j];
            rate *= conc[ion_id][x]/(conc[ion_id][x]+robj[i].Ks_reactants[j]);
            // cout<< conc[ion_id][x]/(conc[ion_id][x]+robj[i].Ks_reactants[j])<<endl;
        }

        for( int j=0; j<robj[i].n_Inh;j++){ /// ion id for inhibition
            string namebac=robj[i].bac_name;
            string nameion= robj[i].inh_ions[j];

            for (int k=0; k< ni;k++ ){        
                if (iobj[k].name==nameion){                    
                    id_i=k;
                    break;
                }
            }


            

            for (int k=0;k<n_inh;k++){
                if(inh_type[k].find(namebac) != string::npos &&inh_type[k].find(nameion) != string::npos ){

                    // cout<< id_i<<" "<<iobj[id_i].name<< endl;
                    rate *= K_inh[k]/(K_inh[k]+conc[id_i][x]);
                    // cout<<K_inh[k]/(K_inh[k]+conc[id_i][x])<<endl;;
                    break;
                }
            }

        }
        
        for (int j=0; j<robj[i].nreactant;j++){
           int ionid=robj[i].reactants_id[j];          
           reac_rate[ionid] += robj[i].stoi_co_reactants[j]* rate;
        //    cout<< reac_rate[ionid]<<endl;
        }

        for (int j=0; j<robj[i].nproduct;j++){
           int ionid=robj[i].products_id[j];          
           reac_rate[ionid] += robj[i].stoi_co_products[j]* rate;
        }
        
    }

    for (int i=0 ; i<ni;i++){
            conc[i][x]=conc[i][x]+dt*reac_rate[i];
        }
    

}




void reaction_enamel(Ion* iobj, double conc_enamel[][N+NT-1],  
 int ni, int x, double E_por[] , double rd_HAP[],double rr_HAP[], double rd_ca[], double rr_ca[], double DSe[], double fap_fraction){

double  IP, DS; // Ilie 2012
double conc[ni];

    int id_ca,id_h,id_pho, id, id_lach,id_pho1n, id_f1n,id_sod;
     // Ca2++, Pho2--, H+ ion's id has to be calculated

    for (int i=0; i< ni;i++ ){
        // cout << iobj[i].name<< endl;

    if(iobj[i].name== "f1n"){
        id_f1n=i;
        break;
    }}
    for (int i=0;i<ni;i++){
    if (iobj[i].name=="ca2p"){        
            id_ca=i;
        break;
    }
    }

    // cout<< id_ca << " "<< ni<< endl;

    for (int i=0; i< ni;i++ ){
    if (iobj[i].name=="pho2n"){
            id_pho=i;
        break;
    }
    }

    for (int i=0; i< ni;i++ ){
    if (iobj[i].name=="h1p"){
            id_h=i;
        break;
    }
    } 

    for (int i=0; i< ni;i++ ){
    if (iobj[i].name=="lach"){
            id_lach=i;
        break;
    }
    } 

    for (int i=0; i< ni;i++ ){
    if (iobj[i].name=="pho1n"){
            id_pho1n=i;
        break;
    }
    } 

    for (int i=0; i< ni;i++ ){
    if (iobj[i].name=="sod1p"){
            id_sod=i;
        break;
    }
    } 





    double IP_HAP=pow(conc_enamel[id_ca][x],5.0)*pow(conc_enamel[id_pho][x],3.0)*
    pow(conc_enamel[id_h][x],-(4.0-xf))*pow(Ke_pho2n,3)*pow(10,-8.0)*pow(conc_enamel[id_f1n][x],xf);
    double DS_HAP=pow(IP_HAP/K_HAP_enamel,1/9.0);
    DSe[x]=DS_HAP;

    string acid[]={"pho1n","lach"};
    double sum=0.0;
    int len=sizeof(acid)/sizeof(acid[0]);
    // cout<< "len= "<<len << endl;
    for(int i=0;i<len;i++){
        for(int j=0; j<ni;j++){
            if(iobj[j].name==acid[i]){
                id=j;
                break;
            }
        }
        
        sum += conc_enamel[id][x];
        
    }

    double IP_FAP=pow(conc_enamel[id_ca][x],5.0)*pow(conc_enamel[id_pho][x],3.0)*
    pow(conc_enamel[id_h][x],-3.0)*pow(Ke_pho2n,3)*conc_enamel[id_f1n][x];
    double DS_FAP=pow(IP_FAP/K_HAP_enamel_F,1/9.0);

    double sigma=1-DS_FAP;


//     cout<<" "<<endl;
// cout << " DS=   "<<DS_HAP<<" "<<conc_plaq[id_ca][N-1] <<endl;
// cout<<" "<<endl;
    double rd_flu=0.0;
    if (DS_HAP<=1)
    {
        if(DS_FAP<=1){
            rd_ca[x]=kd*pow((1-DS_HAP),md)*pow((sum),nd)*(1-fap_fraction)+fap_fraction*kd_F*sigma*sigma;
            rd_flu=fap_fraction*kd_F*sigma*sigma;

        }else{
            rd_ca[x]=kd*pow((1-DS_HAP),md)*pow((sum),nd)*(1-fl_frac);            

        }
         
        // cout<<rd_HAP<<endl;
        // cout<<conc_plaq[id_h][N-1]<<" "<<conc_plaq[id_ca][N-1]<< " "<<conc_plaq[id_pho][N-1]<<endl;
    }else{
         rd_ca[x]=0;
    }  


    double a_d= 2.0/r_rod*(1.0-E_por[x]); // this the specific crystal surface for demin process

    rd_HAP[x]=0.2*rd_ca[x]*a_d; 
    rd_flu=rd_flu*0.2*a_d;
    

/// Similar process for Reminerilasation 
if (DS_HAP>=1 & E_por[x]>0.01)
    {
         rr_ca[x]=kr*pow((DS_HAP-1.0),nr);
        // cout<<rd_HAP<<endl;
        // cout<<conc_plaq[id_h][N-1]<<" "<<conc_plaq[id_ca][N-1]<< " "<<conc_plaq[id_pho][N-1]<<endl;
    }else{
         rr_ca[x]=0;
    }  
    double a_r;

    if(E_por[x]<=0.01){
        a_r=0.0;
        // cout<<E_por<<endl;
    }else{
        a_r=a_d*(1-1/(1+exp(45.0*(E_por[x]-0.15)))); // specific crystal surface for the Remin
    }

    // a_r=a_d*(1-1/(1+exp(45.0*(E_por-0.15)))); // specific crystal surface for the Remin
     

    rr_HAP[x]=rr_ca[x]*a_r;

    E_por[x] += dt*M_HAP/rho_HAP*(rd_HAP[x]-rr_HAP[x]); //porosity change equation
    if(E_por[x]>1.0){
        E_por[x]=1.0;
    }
    // cout<< dt*M_HAP/rho_HAP*(rd_HAP-rr_HAP)<<endl;

    // if (dt*M_HAP/rho_HAP*(rd_HAP-rr_HAP)<0){
    //     cout<< E_por<<endl;
    // }




//  For Hydrohen ion
    conc[id_h]=conc_enamel[id_h][x] -(4.0-xf)*(rd_HAP[x]-rr_HAP[x])*dt;

    if(conc[id_h]<0  ){
        conc[id_h]=0.0001;
        
    }

     if(conc[id_h]>0.01  ){
        conc[id_h]=0.01;        
    }


    // for Fluoride ion

    conc[id_f1n]=conc_enamel[id_f1n][x] +rd_flu*dt;
    if(conc[id_f1n] <0) {
        cerr<<"F";
        exit(EXIT_FAILURE);
        conc[id_f1n]=0.0;

    }
    


    //For Ca ion
    conc[id_ca]=conc_enamel[id_ca][x] +5.0*(rd_HAP[x]-rr_HAP[x])*dt;
    if(conc[id_ca] <0){
        // cerr<<pow((DS_HAP-1.0),nr)*kr*a_r<<endl;
        //  cerr<<"Ca "<<conc_plaq[id_ca][N-2]<<" "<<dx/iobj[id_ca].Diffu*5.0*(rd_ca-rr_ca)<<endl;;
        //  cerr<< rd_ca<<" "<< rr_ca<<" " << dx/iobj[id_ca].Diffu*5.0<<" "<< endl;
        // exit(EXIT_FAILURE);
        conc[id_ca]=0.0;
        }
        // cout<<conc[id_ca]<< " "<<conc_plaq[id_ca][N-2] <<endl;

    // For HPO4--  ion
    conc[id_pho]=conc_enamel[id_pho][x] +3*(rd_HAP[x]-rr_HAP[x])*dt;
    if(conc[id_pho] <0){
        //  cerr<<"pho2--";
        // exit(EXIT_FAILURE);
        conc[id_pho]=0.0;

    } 


    


    for ( int i=0; i<ni;i++){
        if(i != id_ca && i != id_h && i!= id_pho && i!= id_f1n){
            conc[i]=conc_enamel[i][x];
            // cout<< "hello "<<iobj[i].name<< endl;
        }
    }
    // cout << " rd= "<< rd_HAP<<endl;

    for (int i=0; i<ni;i++){
        conc_enamel[i][x]=conc[i];
    }

 }

void diffusion(Ion* iobj, double conc[][N+NT-1], double reaction[], int ni, int x, double E_por[]){
for (int i=0;i<ni;i++){ 
    double alpha;
    if(x>0 & x<N-1){
        alpha=0.25*iobj[i].Diffu/pow(dx,2.0);
        reaction[i]=alpha*(conc[i][x-1]+conc[i][x+1]-2*conc[i][x]);

    }else{
        if(iobj[i].name== "glu"){
            reaction[i]=0.0;
        }else{
            alpha=E_por[x]*iobj[i].Diffu/pow(dxT,2.0);
            double diffu1=(E_por[x+1]+E_por[x])/2.0*iobj[i].Diffu;
            double diffu2=(E_por[x]+E_por[x-1])/2.0*iobj[i].Diffu;
            double alpha1=diffu1/pow(dxT,2.0);
            double alpha2=diffu2/pow(dxT,2.0);
            reaction[i]=alpha1*(conc[i][x+1]-conc[i][x])-alpha2*(conc[i][x]-conc[i][x-1]);
        }
        
    }

    // reaction[i]=alpha*(conc[i][x-1]+conc[i][x+1]-2*conc[i][x]);
    // if(abs(reaction[i])>1.0){
    //         cerr<< "error difffion"<<endl;
    //         cerr<<reaction[i]<<" "<<iobj[i].Diffu/pow(dx,2.0)<<  endl;
    //         cerr<< conc[i][x-1]+conc[i][x+1]-2*conc[i][x]<<endl;
    //         cerr<< conc[i][x-1]<<" "<<conc[i][x]<<" "<<conc[i][x+1]<< " "<<endl;
    //         cerr<< iobj[i].name <<" " <<x<<endl;
    //         exit(EXIT_FAILURE);
    //     }  
   } 
}


void electro_diffusion(Ion* iobj, double conc[][N+NT-1], double reaction[], 
 double phi[],  int ni, int x, double E_por[]){


    for (int i=0;i<ni;i++){  
        double alpha;
         if(x>0 & x<N-1){
            alpha=0.25*iobj[i].Diffu/pow(dx,2.0)*iobj[i].charge*F/R/T;
        }else{
            if(iobj[i].name=="glu"){
                alpha=0;
            }else {
                alpha=E_por[x]*iobj[i].Diffu/pow(dxT,2.0)*iobj[i].charge*F/R/T;

            }
            
            
        }
        
        reaction[i]=alpha/4*(conc[i][x+1]-conc[i][x-1])*
        (phi[x+1]-phi[x-1])+
        alpha*conc[i][x]*(phi[x+1]+phi[x-1]-2*phi[x]); 
        if(abs(reaction[i])>1.0){
            cerr<< alpha<<endl;
            cerr<<reaction[i]<<endl;
            cerr<< phi[x+1]+phi[x-1]-2*phi[x]<<endl;
            cerr<< conc[i][x+1]<<" "<<conc[i][x-1]<<" "<<x<< " "<<i<<endl;
            cerr<< conc[i][x+1]-conc[i][x-1]<<endl;
            exit(EXIT_FAILURE);
        }   
   }
 }



void equillibriate_ions_film( Ion* iobj, Reaction_ion_equi* robj, double conc[], int ni,int nr ){    
    
    int id_H, id_oh;    
    int ir1,ip1,ip2;
    double c1,c2,c3,ka,xc;

    // for(int i=0;i<ni;i++){
    //     if(iobj[i].name=="h1p"){
    //         id_H=i;            
    //         break;
    //     }
    // } 
    
   

        for (int j = 0; j < nr; j++){
            
            ir1=robj[j].reactants_id[0];
            ip1=robj[j].products_id[0];
            ip2=robj[j].products_id[1];  

            c1=conc[ir1];
            c2=conc[ip1];
            c3=conc[ip2];
            ka=robj[j].K_eq;

            xc= (-(c2+c3+ka)+sqrt((c2+c3+ka)*(c2+c3+ka)-4*(c2*c3-ka*c1)))/2.0;

            conc[ir1] -= xc;
            conc[ip1] += xc;
            conc[ip2] += xc; 

        }     

               
        

    }
   
    

void equillibriate_ions(Ion* iobj, Reaction_ion_equi* robj, double conc[][N+NT-1], int ni,int nr, int x){    
    int id_H, id_oh;    
    int ir1,ip1,ip2;
    double c1,c2,c3,ka,xc;
    int nu,nru;

    if(x>=0 & x<=N-1){
        nu=ni;
        nru=nr;
    }else{
        nu=ni-n_fixed_ions;
        nru=nr-nr_fixed;

    }

    // for(int i=0;i<nu;i++){
    //     if(iobj[i].name=="h1p"){
    //         id_H=i;            
    //         break;
    //     }
    // } 

    
   

        for (int j = 0; j < nru; j++){
            
            ir1=robj[j].reactants_id[0];
            ip1=robj[j].products_id[0];
            ip2=robj[j].products_id[1];  

            c1=conc[ir1][x];
            c2=conc[ip1][x];
            c3=conc[ip2][x];
            ka=robj[j].K_eq;

            xc= (-(c2+c3+ka)+sqrt((c2+c3+ka)*(c2+c3+ka)-4*(c2*c3-ka*c1)))/2.0;

            conc[ir1][x] -= xc;
            conc[ip1][x] += xc;
            conc[ip2][x] += xc; 

        }     
       
}


void boundary_update(Ion* iobj, double conc_film[], double conc[][N+NT-1],  
int ni, double E_por[]){

for (int i=0; i< ni-n_fixed_ions;i++ ){
    conc[i][0]=conc_film[i];              //// Bulk-thin Film

}

// for (int i=0; i< ni-n_fixed_ions ;i++ ){
//     double diffu1=iobj[i].Diffu;
//     double diffu2=iobj[i].Diffu*0.25;
//     double dx1= dxF;
//     double dx2=dx;

//     double alpha= diffu1/dx1*dx2/diffu2;               ///// Thin Film- Plaque

//     conc[i][NF-1]=(conc[i][NF]+alpha*conc[i][NF-2])/(1+alpha);            

// } 


for (int i=0; i< ni-n_fixed_ions ;i++ ){

    if(iobj[i].name== "glu"){
        conc[i][N-1]=conc[i][N-2];
    }else  {
        double diffu1=iobj[i].Diffu*(E_por[N-1]+ 0.25)/2.0;
        double diffu2=iobj[i].Diffu*(E_por[N]+E_por[N-1])/2.0;       //// Plaque- Enamel
        double dx1= dx;
        double dx2=dxT;

        double alpha= diffu1/dx1*dx2/diffu2;

        conc[i][N-1]=(conc[i][N]+alpha*conc[i][N-2])/(1+alpha); 


    }
}

for (int i=0; i< ni-n_fixed_ions ;i++ ){
    conc[i][N+NT-2]=conc[i][N+NT-3];        //// Enamel exterme boundary
}
}




void phi_update(Ion* iobj, double conc[][N+NT-1],double phi[], double new_phi[], int ni){

double charge_density;




for (int i=0; i<=N-1; i++){
    charge_density=0;

    for (int j=0; j<ni; j++){
        charge_density += e*iobj[j].charge *conc[j][i];
    }

    
    new_phi[i]=0.5*(phi[i-1]+phi[i+1]+dx*dx/epsilon*charge_density);
}

for (int i=N; i<N+NT-2; i++){
    charge_density=0;

    for (int j=0; j<ni-n_fixed_ions; j++){
        charge_density += e*iobj[j].charge *conc[j][i];
    }

    
    new_phi[i]=0.5*(phi[i-1]+phi[i+1]+dx*dx/epsilon*charge_density);
}


//update boundary phi
new_phi[0]=0; //saliva plaque boundary is nutral

new_phi[N+NT-2]=new_phi[N+NT-3]; //Plaque tooth surface is reflective boundary

}







int id_ion(Ion* iobj, string ion_name, int ni){
    int id;
    for (int i=0; i< ni; i++){
        if(iobj[i].name== ion_name){
            id=i;
            break;

        }
    }
    return id;
 }


void exchange_bulk_film(double conc_bulk[], double conc_film[], double reaction[],int size){
    for (int i=0; i< size; i++){
        reaction[i]=log(2.0)/t_hf*(conc_bulk[i]-conc_film[i]);
    }
}

void exchange_film_plaque(Ion* iobj, double conc[][N+NT-1], double conc_film[], 
double reaction[], int size){
    for (int i=0; i< size; i++){
        reaction[i]=a_vf*iobj[i].Diffu*0.25 *(conc[i][1]-conc[i][0])/dx;
    }
}

void film_conc_update(double conc_film[], double rate1[],double rate2[],int ni){
for(int j=0; j<ni;j++){
                double reaction_term=dt*(rate1[j]+rate2[j]);
                conc_film[j]=conc_film[j]+ reaction_term; 

}
}



void reaction(Ion* iobj, Reaction_ion_bio* robj, double conc[][N+NT-1],
int ni, int nrr, int x,  string inh_type[],  double K_inh[], int n_inh,double E_por[],double rd_HAP[],double rr_HAP[],
 double rd_ca[], double rr_ca[], double DSe[], double fap_fraction){

    

    if(x>= 0 & x<N-1 ){ ///x=0,  x=NF-1 and x=NF+N-2 is the boundary

        reaction_plaque(iobj,robj,conc,ni, nrr,x,inh_type,K_inh,n_inh);

    } else {

        reaction_enamel(iobj,conc, ni-n_fixed_ions,x,E_por,rd_HAP,rr_HAP,rd_ca,rr_ca, DSe,fap_fraction);

    }

    
}

double glu_pulse(double current_time, double glu_conc){
    double bulk_glucose;
     double gluc_prev=0.07;
    int opulse=1;

    if(current_time <t_rest){
            bulk_glucose=0.07;     

    }else{
        double t0=current_time-t_rest;
        int ipulse=t0/t_cycle+1;
        
        if(ipulse <= ncycle){
            if(ipulse>opulse){
            gluc_prev=glu_conc;
            opulse=ipulse;
        }
            double ti=t0-(ipulse-1)*t_cycle;
            if(ti>=0 & ti<t_step){
                bulk_glucose=gluc_prev+(glu_peak-gluc_prev)/t_step*ti;
                
            }else if( ti>=t_step & ti<(t_step+t_feed)) {
                bulk_glucose=glu_peak;
                                   

            }else {
                bulk_glucose=0.07+glu_peak*exp(-log(2.0)/t_hs*(ti-(t_step+t_feed)));
                
            }
        } else{
            double ti=current_time-t_rest-(ncycle-1)*t_cycle;
            bulk_glucose=0.07+glu_peak*exp(-log(2.0)/t_hs*(ti-(t_step+t_feed)));  
        }     


 
    }
return bulk_glucose;

 }
