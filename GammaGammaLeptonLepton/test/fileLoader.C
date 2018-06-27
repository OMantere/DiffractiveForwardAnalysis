#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>

using namespace std;

std::ifstream infile("fileList.txt");

std::string line;

void fileLoader()
{

    vector<string> filelist;
    string line;
    while(getline(infile, line)){
        filelist.push_back(line);
    }

    vector<string> names;
    vector<double> xsecs;
    vector<string> outputs;

    for(int i=0; i < filelist.size(); i++){
        vector<string> words;
        
        istringstream iss(filelist[i]);

        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(words));
        names.push_back(words[0]);
        double xsec;
        std::string::size_type sz;
        xsec = std::stof(words[1], &sz);
        xsecs.push_back(xsec);
        outputs.push_back(words[2]);


    }

            
}
