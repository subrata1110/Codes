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
#include "class.h"
// #include "reaction_ion_equi.h"
// #include "reaction_ion_bio.h"
#include "parameters.h"
#include "functions.h"


using namespace std;





 

int main(){
    // Reading Ion file and create objects neccessary for further calcuation 

    string Ionfile= "ion.txt";
    string react_equi="reaction.txt";
    string react_bio= "reaction_bio.txt";
    string ionname;
    ifstream ionfile(Ionfile);

    if(!ionfile.is_open()){
        cerr << "Error opening file" << endl;
        return 1;
    }

    int n_ion=0;
    string line;
    getline(ionfile,line); // neglecting the first line
    while (getline(ionfile,line))
    {
        n_ion++;
    }

    ionfile.close();
    ionfile.open(Ionfile);

    getline(ionfile,line); // ignore the first line

    Ion* iobj= new Ion[n_ion];

    for (int i=0; i<n_ion; i++){
        Ion& ilocal=iobj[i];
        getline(ionfile,line);
        istringstream iss(line);
        iss >> ilocal.name>> ilocal.conc_intial >> ilocal.charge >> ilocal.Diffu ;
        // ilocal.Diffu *= 0.25;

        ilocal.id=i;
    }

    ionfile.close();
    cout << "ion file read"<<endl;

    //// Now read the reaction file for equillibration

    ifstream reaceq_file(react_equi);

    getline(reaceq_file,line); // igonre the first line

    int nr_eq=0;

    while(getline(reaceq_file,line)){
        nr_eq++;
    }

    reaceq_file.close();
    reaceq_file.open(react_equi);

    getline(reaceq_file,line); // igonore the first line

    Reaction_ion_equi* robj_eq= new Reaction_ion_equi[nr_eq];

    for (int i=0; i< nr_eq;i++){
        Reaction_ion_equi& rlocal=robj_eq[i];

        getline(reaceq_file,line);
        istringstream iss(line);
        iss >> rlocal.K_eq >> rlocal.nreactant>>rlocal.nproduct;
        rlocal.reactants_id=new int[rlocal.nreactant];
        rlocal.products_id=new int[rlocal.nproduct];
        
        for (int j=0; j< rlocal.nreactant;j++){
            iss >> ionname;
            for (int k=0; k<n_ion;k++){
                if(iobj[k].name==ionname){
                    rlocal.reactants_id[j]= iobj[k].id;
                    cout<< ionname<< " reactanats "<< rlocal.reactants_id[j]<< endl;
                    break;
                }
            }
        }


        for (int j=0; j< rlocal.nproduct;j++){
            iss >> ionname;
            for (int k=0; k<n_ion;k++){
                if(iobj[k].name==ionname){
                    rlocal.products_id[j]= iobj[k].id;
                    cout<< ionname<< " product "<< rlocal.products_id[j]<< endl;
                    break;
                }
            }
        }

    }

    reaceq_file.close();
    cout << "reaction file read"<<endl;


    ifstream reabio_file(react_bio);

    int nr_bio=0;

    getline(reabio_file,line); // igonre the first line

    while(getline(reabio_file,line)){
        nr_bio++;
        
    }

    reabio_file.close();
    reabio_file.open(react_bio);

    Reaction_ion_bio* rbioobj =new Reaction_ion_bio[nr_bio];
    getline(reabio_file,line); // ignore the first line


    for (int i=0; i<nr_bio; i++){
        getline(reabio_file,line);
        istringstream iss(line);
        Reaction_ion_bio& rlocal=rbioobj[i];

        iss >> rlocal.bac_name>> rlocal.bac_conc>> rlocal.nreactant >>rlocal.nproduct ;

        rlocal.reactants_id=new int[rlocal.nreactant];
        rlocal.products_id=new int[rlocal.nproduct];
        rlocal.stoi_co_reactants=new double[rlocal.nreactant];
        rlocal.stoi_co_products=new double[rlocal.nproduct];
        rlocal.Ks_reactants= new double[rlocal.nreactant];
        
        for (int j=0; j< rlocal.nreactant;j++){
            iss >> ionname;
            for (int k=0; k<n_ion;k++){
                if(iobj[k].name==ionname){
                    rlocal.reactants_id[j]= iobj[k].id;                    
                    break;
                }
            }
        }


        for (int j=0; j< rlocal.nproduct;j++){
            iss >> ionname;
            for (int k=0; k<n_ion;k++){
                if(iobj[k].name==ionname){
                    rlocal.products_id[j]= iobj[k].id;
                    break;
                }
            }
        }

        for (int j=0; j< rlocal.nreactant;j++){
            iss >> rlocal.stoi_co_reactants[j];            
        }


        for (int j=0; j< rlocal.nproduct;j++){
            iss >> rlocal.stoi_co_products[j];            
        }

        for (int j=0; j< rlocal.nreactant;j++){
            iss >> rlocal.Ks_reactants[j];            
        }



        iss>> rlocal.n_Inh;

        rlocal.inh_ions= new string[rlocal.n_Inh];


        for (int j=0; j<rlocal.n_Inh;j++){

            iss >> rlocal.inh_ions[j]; 

        }

        iss>> rlocal.q_reac; 


    }

    reabio_file.close();

    cout << "reaction bio file read"<<endl;

    ifstream file("inhibition_constant.txt");

    int n_inh_reac=0;
    getline(file,line); // igonre the first line

    while (getline(file,line))
    {
        n_inh_reac++;
    }

    file.close();
    file.open("inhibition_constant.txt");

    getline(file,line); //igonre the first line

    string inh_type[n_inh_reac];
    double K_inh[n_inh_reac];
    

    for (int i=0; i<n_inh_reac;i++){
        getline(file,line);
        istringstream iss(line);
        iss >> inh_type[i]>> K_inh[i];
        
    }
    file.close();

    cout << "inhibition file read" << n_inh_reac<<endl;


    /// intitalise the concentration of ions

    double conc_bulk_saliva[n_ion-n_fixed_ions]; // define the ion concentrations
    double conc_film[n_ion-n_fixed_ions];
    double conc[n_ion][N+NT-1],conc_temp[n_ion][N+NT-1] ; 
    
    double ratef1[n_ion], ratef2[n_ion], ratef3[n_ion];
    
    double phi[N+NT-1], phi_new[N+NT-1];    

    


    for (int i=0; i< n_ion-n_fixed_ions;i++){
        conc_bulk_saliva[i]=iobj[i].conc_intial; 
        conc_film[i]= iobj[i].conc_intial;         

    }

    for (int i=0; i<N+NT-1;i++){
         if(i >=0 & i<N){
             for (int j=0; j<n_ion;j++){
                conc[j][i]=iobj[j].conc_intial;
            }

        } else{
            for (int j=0; j<n_ion-n_fixed_ions;j++){
                if(iobj[j].name=="glu"){
                    conc[j][i]=0;
                }else{
                    conc[j][i]=iobj[j].conc_intial;
                }
                
            }
            
        }
    }
    int id_tcp=id_ion(iobj,"tcp",n_ion);
    for (int i=0; i<N+NT-1;i++){
        conc[id_tcp][i]=0.0;

    }




    for (int i=0; i<N+NT-1; i++){
        phi[i]=0;

    }

    



    //////////////
    /// the intial parameters


    int id_glu=id_ion(iobj,"glu",n_ion);    
    int id_h= id_ion(iobj,"h1p",n_ion);
    int id_ca2p=id_ion(iobj,"ca2p",n_ion);
    int id_f1n=id_ion(iobj,"f1n",n_ion);
    int id_pho2n=id_ion(iobj,"pho2n",n_ion);
    int id_k1p= id_ion(iobj,"k1p",n_ion);
    for (int i=0; i<n_ion;i++){
        cout << iobj[i].name<< " "<< iobj[i].id<< " dt*D/dx**2= "<< dt*iobj[i].Diffu/dx/dx<< endl;
        if(dt*iobj[i].Diffu/dxF/dxF>=0.5 ){
            cerr << dx << " "<< iobj[i].conc_intial*iobj[i].Diffu/1.2*1e6<< endl;
            exit(EXIT_FAILURE);
        }        

    }

    if( dx>iobj[id_h].conc_intial*iobj[id_h].Diffu/1.2*1e6){
            cerr << dx << " "<< iobj[id_h].conc_intial*iobj[id_h].Diffu/1.2*1e6<< endl;
            // exit(EXIT_FAILURE);
        }

    

    


/////////////////////////////////////////////////////
///////////SIMULATION//////////////////////////////
////////////////////////////////////////////////
    // Now start the the system for time dynamics
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
auto start_time = std::chrono::high_resolution_clock::now();

ofstream outfile(output_filename);
ofstream outfile_rd(output_remin_demin);

outfile<< "Time in min"<<","<<"position"<<",";
for(int i=0;i<n_ion;i++){
    outfile<< iobj[i].name <<",";
    }
    
    outfile <<endl;


outfile_rd << "Time in min"<<","<< "Position"<< ","<<"rd_HAP"<<","<<"rr_HAP"<<","<<"rd_Ca"
<<","<<"rr_ca"<<"," <<"Porosity"<<","<<"DS"<<"," <<endl;
   
    double current_time=0.0;
    double E_porosity[N+NT-1]; //porosity 
    for (int i=0; i<N-1; i++){
        E_porosity[i]=0.0;
    }
    for (int i=N-1; i<N+NT-1; i++){
        E_porosity[i]=0.01;
    }


    double rate_dHAP[N+NT-1],rate_rHAP[N+NT-1], rate_dca[N+NT-1],rate_rca[N+NT-1], DS_enamel[N+NT-1];

    // outfile<< "Time in min"<<",";
    //         for(int i=0;i<n_ion;i++){
    //             outfile<< iobj[i].name <<",";
    //         }
    //         outfile<< "glu_bulk"<<","<<"glu_film"<<"," << "glu_teeth"<<","<<"pH"<<",";
    //         outfile<< endl;


    cout << "Started simulation"<<endl;
    // cout<<"k in enamel x=0 = "<< conc_enamel[id_k1p][0]<<endl;
    
    for(int it=0; it<=simu_time/dt;it++){

        


        /// Glucose ion id 

        // int id_glu=id_ion(iobj,"glu",n_ion);
        // int id_h= id_ion(iobj,"h1p",n_ion);

        // double t_rest=10.0; // in sec
        // double t_step=10.0;
        // double t_feed=120.0;

        // // Glucose pulse for first 2 min is without pulse 
        // if(glu_status==1){
        //     if(current_time<t_rest){
        //         conc_bulk_saliva[id_glu]=0.07;
        //     } else if (current_time >=t_rest && current_time<(t_rest+t_step)){
        //          conc_bulk_saliva[id_glu]=0.07+560.0/t_step*(current_time-t_rest);
        //     } 
        //     else if (current_time >=(t_rest+t_step) && current_time<(t_rest+t_step+t_feed)) {
        //         conc_bulk_saliva[id_glu]=0.07+560.0;
        //     }else{
        //         conc_bulk_saliva[id_glu]=0.07+560.0*exp(-log(2)/t_hs*(current_time-(t_rest+t_step+t_feed)));
        //     }

        // }

        conc_bulk_saliva[id_glu]=glu_pulse(current_time,conc_bulk_saliva[id_glu]);

         



        


        


        /// First the reaction of the plaque


        

        double rate1[n_ion], rate2[n_ion];
    exchange_bulk_film(conc_bulk_saliva,conc_film,rate1,n_ion);
    exchange_film_plaque(iobj,conc,conc_film,rate2,n_ion);

    // cout<< rate1[id_glu]<<" "<<rate2[id_glu]<<endl;

    // double conc_film_updated[n_ion-n_fixed_ions];

    film_conc_update(conc_film, rate1,rate2,n_ion);

    equillibriate_ions_film(iobj,robj_eq, conc_film,
            n_ion-n_fixed_ions,nr_eq-nr_fixed);
    
    // cout << conc_bulk_saliva[id_glu]<< " "<< conc_film[id_glu]<<endl;

    
           double charge_density=0;

        for (int j=0; j<n_ion-n_fixed_ions; j++){
            if(iobj[j].name != "k1p"){
                charge_density=charge_density+iobj[j].charge *conc_film[j];
            }
        }

        // calclating the indices of a potassium
        int id_K1p;
        id_K1p=id_ion(iobj,"k1p",n_ion-n_fixed_ions);
        

        conc_film[id_K1p]=-charge_density/iobj[id_K1p].charge;

        
    
        // cout<<" fine " << n_fixed_ions<<" "<<n_ion<<" "<<nr_eq<<" "<< nr_bio<<" "<<" "<< n_inh_reac<< endl;

       #pragma omp parallel for shared(iobj,rbioobj,conc,E_porosity,rate_dHAP,rate_rHAP,rate_dca,rate_rca,DS_enamel) 
        for (int ix=0;ix<N+NT-1;ix++){
            // cout<<ix<<endl;
            
            reaction(iobj,rbioobj, conc,n_ion,nr_bio,ix,inh_type,K_inh,n_inh_reac,E_porosity,rate_dHAP,
            rate_rHAP,rate_dca,rate_rca,DS_enamel);
          
        }
       #pragma omp barrier
    //    cout << conc_film[id_glu]<<endl;

// cout<<" fine " << conc[id_glu][1]<< endl;


       int nxi;
       #pragma omp parallel for shared(conc, iobj,E_porosity) private(ratef1,ratef2)
            for (int ix=0;ix<N+NT-1;ix++){
            if(ix !=0 & ix != (N-1) & ix != (N+NT-2)  ){
                diffusion(iobj,conc,ratef1,n_ion,ix,E_porosity);
                electro_diffusion(iobj,conc,ratef2,phi,n_ion,ix,E_porosity);

                if(nxi >= 0 & nxi<= N-1 ){
                    nxi=n_ion;

                }else{
                   nxi=n_ion-n_fixed_ions;
                }

                

                for (int i=0;i<nxi; i++){
                    conc_temp[i][ix]=conc[i][ix]+(ratef1[i]+ratef2[i])*dt;

                }


            }
        }
       #pragma omp barrier

       for (int ix=0; ix<(N+NT-1);ix++){
         if(ix >= 0 & ix<= N-1 ){
            nxi=n_ion;

        }else{
            nxi=n_ion-n_fixed_ions;
        }

        for (int i=0;i<nxi; i++){
            conc[i][ix]=conc_temp[i][ix];

        }
       }


       boundary_update(iobj,conc_film,conc,n_ion,E_porosity);









       
    

     
// cout<< it<<" k in enamel x=0 = "<< conc_enamel[id_k1p][0]<<endl;



       #pragma omp parallel for shared(conc,robj_eq,iobj) 
        for (int ix=0;ix<N+NT-1;ix++){

            equillibriate_ions(iobj,robj_eq,conc,n_ion,nr_eq,ix);

                 
        }

       #pragma omp barrier  


       
    // Now nutralise charge by calcluating the potasium ions

         charge_density;
        #pragma omp parallel for shared(conc,iobj) private(charge_density) 
        for (int i=0; i<N+NT-1;i++){
            charge_density=0;
            for (int j=0; j<n_ion-n_fixed_ions; j++){

                if(iobj[j].name != "k1p"){
                    charge_density=charge_density+iobj[j].charge *conc[j][i];
                }
            }
            int id_K1p;
            id_K1p=id_ion(iobj,"k1p",n_ion-n_fixed_ions);
            conc[id_K1p][i]=-charge_density/iobj[id_K1p].charge;            

        }
        #pragma omp barrier  

        

       

        phi_update(iobj,conc,phi,phi_new,n_ion);  

        for (int i=0; i<(N+NT-1);i++){
            phi[i]=phi_new[i];
        }     

        

        
       
// if(it%120000   == 0 ){
//             cout<<current_time/60.0 <<" "<< conc_film[id_glu] <<endl;
// }


        if(it%1200000   == 0 ){
            cout<<current_time/60.0 <<" "<< conc_film[id_glu] <<endl;

            outfile<<current_time/60.0<<","<< -150 <<",";
            for(int i=0;i<n_ion-n_fixed_ions;i++){
                outfile<< conc_bulk_saliva[i]<<",";
            }
            outfile<<endl;

            outfile<<current_time/60.0<<","<< -100 <<",";
            for(int i=0;i<n_ion-n_fixed_ions;i++){
                outfile<< conc_bulk_saliva[i]<<",";
            }
            outfile<<endl;

            outfile<<current_time/60.0<<","<< -100 <<",";
            for(int i=0;i<n_ion-n_fixed_ions;i++){
                outfile<< conc_film[i]<<",";
            }
            outfile<<endl;

            outfile<<current_time/60.0<<","<< 0 <<",";
            for(int i=0;i<n_ion-n_fixed_ions;i++){
                outfile<< conc[i][0] <<",";
            }
            outfile<<endl;




            for (int ipos=0;ipos<(N-1);ipos++) {
                outfile<<current_time/60.0<<","<< ipos*dx*pow(10,6.0) <<",";
            for(int i=0;i<n_ion;i++){
                outfile<< conc[i][ipos] <<",";
            }
            outfile<<endl;

            } 

            
            for (int ipos=N-1;ipos<(NT+N-1);ipos++) {
                outfile<<current_time/60.0<<","<< ((N-1)*dx+(ipos-(N-1))*dxT)*pow(10,6.0) <<",";
            for(int i=0;i<n_ion-n_fixed_ions;i++){
                outfile<< conc[i][ipos] <<",";
            }
            outfile<<endl;

            } 
            

            
            for (int ipos=(N-1);ipos<(N+NT-1);ipos++) {
              
                outfile_rd<<current_time/60.0<<","<<((ipos-N+1)*dxT)*pow(10,6.0)<<"," <<rate_dHAP[ipos]<<","<<rate_rHAP[ipos]<<","<<rate_dca[ipos]<<","
                <<rate_rca[ipos]<<","<< E_porosity[ipos]<<"," << DS_enamel[ipos]<<"," <<endl;

            }

            

        }

        current_time=current_time+dt;


    }

    delete [] iobj;
    delete [] robj_eq;
    delete [] rbioobj;

    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    // Print the duration
    std::cout << "Time taken by code: " << duration.count()*1.0e-6 << " seconds" << std::endl;

 outfile.close();
 outfile_rd.close();

}






