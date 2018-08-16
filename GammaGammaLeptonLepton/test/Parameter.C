#include <string>
#include <limits>
#include <iostream>

#include "Parameter.h"

using namespace std;

Parameter::Parameter(string name, double *pAddr){
    //cout << "Creating parameter" << endl;

    
    pNam = name;
    pAdd = pAddr;
    cutPoints = 0;
    high = std::numeric_limits<double>::infinity();
    low = -high;
    inc = true;
    //

}

void Parameter::setRange(double llow, double hhigh, bool incl){
    low = llow;
    high = hhigh;
    inc = incl;
    return;
}

void Parameter::setnPoints(int points){
    nPoints = points;
    return;
}

bool Parameter::cut(double val){
    bool c = low <= val && val <= high;
    bool out =  (c && inc)||(!c && !inc);
    if(!out) cutPoints += 1;
    return out;
};

int Parameter::printCuts(){
   if(cutPoints > 0){ 
    cout << cutPoints <<" ("<< cutPoints*100.0/nPoints <<"\%)" << "\t points did not pass '" << low << " < " << pNam << " < "  << high << "'\t selection" << endl;
   }
    return cutPoints;
}

void Parameter::changeAddr(double *addr){
    cutPoints = 0;
    pAdd = addr;
    return;
}
