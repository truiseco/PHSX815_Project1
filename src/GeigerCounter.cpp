/*************************************************************
* @author   Triston Ruiseco
* @file     GeigerCounter.cpp
* @date     02/22/2021
* @brief    Generates an exponentially distributed data sample and exports it
            to a file according to users specifications.
*************************************************************/

// Std Includes
#include <iostream>
#include <fstream>
#include <cstring>

// Local Includes
#include "Random.h"

using std::string, std::stoi, std::stod, std::stol, std::cout, std::ofstream;

// Program-specific helper function declarations
/*
* param a: first string to be compared
* param b: second string to be compared
* return 1 if a & b are identical character-wise and in length, 0 otherwise
*/
bool strsame(string a, string b);

// Begin primary program function
int main(int argc, char** argv){

  // Command line options parsing variable declarations
  bool printhelp = 0;
  bool argexists = 0;
  long seed = 314159;
  string filename = "data.txt";
  bool output = 0;
  double rate = 1;   // Rate parameter for exponential distribution
  int num_meas = 1;  // Number of measurements to simulate

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    if(strsame(argv[i],"--help")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"-h")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"-seed")){
      argexists = 1;
      seed = stol(argv[++i]);
    }
    if(strsame(argv[i],"-rate")){
      argexists = 1;
      double r = stod(argv[++i]);
      if(r > 0.0){
	       rate = r;
      }
    }
    if(strsame(argv[i],"-measures")){
      argexists = 1;
      int arg_num = stoi(argv[++i]);
      if(arg_num > 0){
	       num_meas = arg_num;
      }
    }
    if(strsame(argv[i],"-output")){
      argexists = 1;
      filename = string(argv[++i]);
      output = 1;
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }

  /* Print the executable usage instructions if the user adds -h or --help
     or provides an undefined option */
  if(printhelp || !output){
    cout << "Usage: " << argv[0] << " [options]\n"
         << "  options:\n"
         << "   --help(-h)          print options\n"
         << "   -seed [number]      random seed to use\n"
         << "   -rate [number]      rate radioactive events (per second)\n"
         << "   -measures [number]  number of time measurements\n"
         << "   -output [filename]  name of ouptut file\n";

    return 0;
  }

  // Instantiate random object to generate random sample distribution
  Random random(seed);

  // Generate and write to file an exponentially distributed sample
  double percent = 0.0;                 // Variables to help print progress
  double denom = double(num_meas)/100;  // as generation of large samples files
                                        // can be very time consuming
  ofstream outFile;
  outFile.open(filename);
  outFile << "rate: " << rate << "\n";
  for(int m = 0; m < num_meas; ++m){
    percent = double(m)/denom;                  // Print progress
    if(int(percent) == percent){                // --------------
      cout << int(percent) << "% complete.\n";  // --------------
    }                                           // --------------
    outFile << random.Exponential(rate) << " ";
  }
  outFile.close();

  return 0;
}

// Program-specific helper function definitions
bool strsame(string a, string b){
  if(a.length()==b.length()){
    int n = a.length();
    for(int i = 0; i < n; i++){
      if(a.at(i)!=b.at(i)){
        return 0;
      }
    }
    return 1;
  }
  return 0;
}
