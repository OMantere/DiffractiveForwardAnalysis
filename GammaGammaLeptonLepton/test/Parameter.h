#ifndef __PARAMETER_H_INCLUDED__
#define __PARAMETER_H_INCLUDED__

#include <string>
#include <limits>
#include <iostream>


using namespace std;

class Parameter{ 
    
    public:
        Parameter(string name="def", double *pAddr=NULL);
       // Parameter(string pNam, double *pAddr);
 
        double *pAdd;
        bool cut(double val = 0);
        int printCuts();
        string pNam;
        void setRange(double llow = 0, double hhigh= 0);
        void changeAddr(double *addr=NULL);
        double low;
        double high;
    private:
        int cutPoints;
};


#endif
