#ifndef ION_H
#define ION_H

#include <string>
using namespace std;

class Ion {
public:
    string name;   
    double conc_intial;
    double Diffu;
    double charge;
    int id;

    Ion(){}    
    
    ~Ion(){}
};



class Reaction_ion_equi {
public:
    double K_eq;

    int nreactant;
    int nproduct;

    int* reactants_id;
    int* products_id;

    int id;

    Reaction_ion_equi(){}
    Reaction_ion_equi(int nreactant, int nproduct){
        this->nreactant = nreactant;
        this->nproduct = nproduct;
        reactants_id = new int[nreactant];
        products_id = new int[nproduct];
    }
    ~Reaction_ion_equi() {
        delete [] reactants_id;
        delete [] products_id;
    }
};

// this class will creates objects for the biological process occuring in the plaque upon reading the reaction
//file in the plaque
class Reaction_ion_bio {
public:

    string bac_name; 
    double bac_conc;   
    int nreactant;
    int nproduct;   
    int n_Inh;
    

    int* reactants_id;
    double* stoi_co_reactants;
    int* products_id;
    double* stoi_co_products;   
    string* inh_ions;
    double* Ks_reactants;

    double q_reac;
    int id;

    Reaction_ion_bio(){}
    Reaction_ion_bio(int nreactant, int nproduct, int n_Mc,int n_Inh){
        this->nreactant = nreactant;
        this->nproduct = nproduct;      
        this-> n_Inh= n_Inh;
        reactants_id = new int[nreactant];
        stoi_co_reactants=new double[nreactant];
        products_id = new int[nproduct];        
        stoi_co_products=new double[nproduct];       
        inh_ions= new string[n_Inh];
        Ks_reactants= new double[nreactant];
    }
    ~Reaction_ion_bio() {
        delete [] reactants_id;
        delete [] products_id;
        delete [] stoi_co_products;
        delete [] stoi_co_reactants;       
        delete [] inh_ions;
        delete [] Ks_reactants;
    }
};





class ion_sys {
    public:
    double conc;
    ion_sys(){}
    ~ion_sys(){}
};

#endif