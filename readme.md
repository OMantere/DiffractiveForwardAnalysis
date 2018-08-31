#Running the code

The code is run with the command 

    root -l Parameter.C combinedHisto.C

The code then takes instructions from a .txt file (filelist2.txt) where the names of datasets, cross sections, errors, names of datasets and number of events in that file are included. 

The options for combinedHisto.C are included in the top of the file. There are 4 booleans corresponding to different plots. plot1 produces 3x2 matrix of 6 first elements in the paramList-vector. The paramList contains the parameters that are loaded in to the program for which one can make cuts. 

plot2 does not work in the current version. It is supposed to make 2d COLZ histograms.

plot3 does "scanplots" where the first data-entry is used as a background and the following entries in the filelist are used as signals that are compared to the background. The parameter used for this is param2.

plot4 plots only the first parameter in paramList. plot1 needs to be enabled for this to work.

The cuts are defined within the setRanges -function. The syntax is
    setParamRange(pList, lowerbound, upperbound, inclusive/exclusive);
The pList is just an array with the pointers to the parameters. The lowerbound and upperbound define the region that is cut. The inclusive/exclusive is a boolean where true means that the region is cut away and false means that everything else is cut away.

If something was left unclear, please contact me at    kristian.arjas95@gmail.com
