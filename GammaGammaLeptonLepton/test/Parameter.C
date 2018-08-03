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
    //

}

void Parameter::setRange(double llow, double hhigh){
    low = llow;
    high = hhigh;
    return;
}

bool Parameter::cut(double val){
    bool c = low <= val && val <= high;
    if(!c) cutPoints += 1;
    return c;
};

int Parameter::printCuts(){
   if(cutPoints > 0){ 
    cout << cutPoints << "\t points did not pass '" << low << " < " << pNam << " < "  << high << "'\t selection" << endl;
   }
    return cutPoints;
}

void Parameter::changeAddr(double *addr){
    cutPoints = 0;
    pAdd = addr;
    return;
}
