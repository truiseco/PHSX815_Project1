/*************************************************************
* @author   Triston Ruiseco
* @file     ExpHypoTest.cpp
* @date     02/22/2021
* @brief    Analyzes two exponentially distributed data samples and simulates
            highly configurable tests for various measurements
            per experiment then plots the signifiance of the tests vs the
            measurements per experiment used to generate them.
*************************************************************/

// Std Includes
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// ROOT Includes
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"

// Directives and declarations for namespaces and namespace members
using std::cout, std::cin, std::string, std::vector, std::stoi, std::ifstream,
      std::sort, std::erf, std::sqrt, std::stringstream, std::exp,
      std::abs, std::min, std::max;

// Program-specific helper function declarations
/*
* param a: first string to be compared
* param b: second string to be compared
* return 1 if a & b are identical character-wise and in length, 0 otherwise
*/
bool strsame(string a, string b);

/*
* param x: data
* param rate: rate parameter of exp distribution
* return probability of observing x in exp distribution
*/
double ExpPDF(double x, double rate);

/*
* param arr0: vector of some distribution for some hypothesis 0
* param arr1: vector of some distribution for some hypothesis 1
* return significance level of test such that alpha = beta
*/
double FindABSame(const vector<double> &arr0, const vector<double> &arr1);

/*
* return normPDF(x,0,1)
*/
double sigma(int x);

/*
* return first index of arr such that arr[index] is of lesser value than y
*/
int FirstIndexLess(const vector<double>& arr, double y);

void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title);

// Begin primary program function
int main(int argc, char** argv){

  // Command line options parsing variable declarations
  bool printhelp = 0;
  bool argexists = 0;
  bool haveInput[2] = {0,0};
  string InputFile[2] = {"\0","\0"};
  int Nexp = 0;
  int mpe = 0;
  int step = 1;

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    argexists = 0;
    if(strsame(argv[i],"--help") || strsame(argv[i],"-h")){
      argexists = 1;
      printhelp = 1;
    }
    if(strsame(argv[i],"-mpe")){
      argexists = 1;
      mpe = stoi(argv[++i]);
      if(mpe <= 0){
        printhelp = 1;
        cout << "-mpe must be a positive integer\n";
      }
    }
    if(strsame(argv[i],"-Nexp")){
      argexists = 1;
      Nexp = stoi(argv[++i]);
      if(Nexp <= 0){
        printhelp = 1;
        cout << "-Nexp must be a positive integer\n";
      }
    }
    if(strsame(argv[i],"-H0")){
      argexists = 1;
      InputFile[0] = string(argv[++i]);
      haveInput[0] = 1;
    }
    if(strsame(argv[i],"-H1")){
      argexists = 1;
      InputFile[1] = string(argv[++i]);
      haveInput[1] = 1;
    }
    if(strsame(argv[i],"-step")){
      argexists = 1;
      step = stoi(argv[++i]);
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }

  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide input, or provides an undefined option */
  if(printhelp || !haveInput[0] || !haveInput[1]){
    cout << "Usage: " << argv[0] << " [options] -Nexp [integer] -mpe [integer] -H0 [filename] -H1 [filename]\n"
         << "  descriptions:\n"
         << "   -Nexp [integer]       number of experiments per test\n"
         << "   -mpe [integer]        max measurements/experiment\n"
         << "   -H0 [filename]        input data for hypothesis 0\n"
         << "   -H1 [filename]        input data for hypothesis 1\n"
         << "  options:\n"
         << "   --help(-h)            print options\n"
         << "   -step                 measurements/experiment increment";

    return 0;
  }

  // Experimental data storage objects
  double rate[2] = {0.0,0.0};
  vector<double> times[2];
  string testr = "\0";
  double temptime = 0.;

  // Read in experimental data
  for(int h = 0; h < 2; ++h){
    ifstream inFile(InputFile[h]);

    // Check if file can be opened
    if(!inFile.is_open()){
      cout << "Failed to open " << InputFile[h] << "\n";
      return 0;
    }

    // Check if file includes necessary valid rate parameter
    inFile >> testr;
    if(!strsame(testr,"rate:")){
      cout << "Input file " << InputFile[h] << " formatted improperly\n";
      inFile.close();
      return 0;
    } else{
      inFile >> rate[h];
      if(rate <= 0){
        cout << "Input File " << InputFile[h] << " contains invalid rate parameter\n";
        inFile.close();
        return 0;
      }
    }

    // Store all experimental data
    cout << "Reading data set " << h << "...\n";
    for(int i = 0; i < Nexp*mpe && inFile >> temptime; ++i){
      times[h].push_back(temptime);
    }
    inFile.close();

    // Check that enough data was found to perform requested calculations
    if(times[h].size() < (Nexp*mpe)){
      cout << InputFile[h] << " contains too few measurements to complete analysis\n";
      return 0;
    }
  }

  // Data analysis storage objects
  vector<double> LLR[2];    // log-likelihood of experiment for each h
  vector<double> CR;        // significance level of each generated test
  vector<double> MPE;       // keeps index-wise track of  measures/exp for CR
  double tempLLR = 0.0;


  /* Construct vector of log likelihood ratios for each experiment,
     for each hypothesis, for each value of measurements/experiment, */
  cout << "Analyzing data...\n";
  for(int M = 1; M < mpe; M+=step){    // loop over all measurements/experiment
    for(int h = 0; h < 2; ++h){        // loop over all hypotheses
      LLR[h].clear();
      for(int e = 0; e < Nexp; ++e){   // loop over all experiments
        tempLLR = 0;
        for(int m = 0; m < M; ++m){    // loop over all measurements
          tempLLR += log(ExpPDF(times[h][(M*e)+m], rate[1])); // M*e used to prevent reusing data
          tempLLR -= log(ExpPDF(times[h][(M*e)+m], rate[0])); // across experiments within the same test
        }
        LLR[h].push_back(tempLLR);
      }
    }
    cout << "Step " << ((M-1)/step)+1  << " of " << mpe/step << " complete.\n";
    // Sort test distributions
    sort(LLR[0].begin(),LLR[0].end());
    sort(LLR[1].begin(),LLR[1].end());

    // Find and record alpha
    CR.push_back(FindABSame(LLR[0], LLR[1]));
    MPE.push_back(double(M)/1000.0);
  }

  // Plot and save results
  stringstream title;
  title << Nexp << " experiments per test with rates " << rate[0] << ", " <<rate[1] << " events / second";
  PlotConfidence(MPE, CR, title.str());

  return 0;
}



// Program-specific helper function definitions
double sigma(int x){
  return erf(double(x)/sqrt(2));
}

int FirstIndexLess(const vector<double>& arr, double y){
  int n = arr.size();
  for(int i = 0; i < n; ++i){
    if(arr[i] < y){
      return i;
    }
  }
  return n;
}

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

double ExpPDF(double x, double rate){
  return (rate*exp(-rate*x));
}

double FindABSame(const vector<double>& arr0, const vector<double>& arr1){
  int n0 = arr0.size();
  int n1 = arr1.size();
  vector<int> x[2];

  int b = 0;
  for(int a = n0 - 1; a > -1; --a){
    b = 0;
    for(int i = 0; i < n1; ++i){
      if(arr1[i+1] < arr0[a]){
        ++b;
      }
    }
    x[0].push_back(a);
    x[1].push_back(abs(b + a - n0));
  }
  int k = x[1].size();
  int mindex = 0;
  for(int i = 0; i < k; ++i){
    if(x[1][i] < x[1][mindex]){
      mindex = i;
    }
  }
  return (1.0 - (double(x[0][mindex])/double(n0)));
}

void PlotConfidence(vector<double>& arr0, vector<double>& arr1, const string& title){
  TCanvas* canvas = (TCanvas*) new TCanvas("canvas", "Canvas_Title",200,10,1000,800);
  canvas->SetGrid();


  double lm = 0.15;
  double rm = 0.05;
  double bm = 0.1;
  double tm = 0.1;
  canvas->SetLeftMargin(lm);
  canvas->SetRightMargin(rm);
  canvas->SetBottomMargin(bm);
  canvas->SetTopMargin(tm);
  canvas->SetLogy();
  canvas->Draw();
  canvas->Update();

  int N = arr0.size();
  TGraph* graph = new TGraph(N, &(arr0[0]), &(arr1[0]));
  graph->SetLineColor(kAzure+1);
  graph->SetLineWidth(2);
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(0);
  graph->SetTitle(title.c_str());
  graph->GetXaxis()->SetTitle("1000s of Measurements/Experiment");
  graph->GetYaxis()->SetTitle("Test Significance #alpha (= #beta)");
  graph->Draw();
  canvas->Update();

  TLatex text;
  text.SetTextFont(42);
  text.SetTextSize(0.03);
  text.SetTextAlign(33);
  text.SetTextAngle(90);
  text.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);

  int n = arr1.size();
  int i = 1;
  int k = FirstIndexLess(arr1, 1 - sigma(i));

  while(k < n && i < 8){
    double xndc = (1.-rm-lm)*((arr0[k]-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+lm;
    line->SetLineColor(kRed+3-i);
    line->DrawLineNDC(xndc,bm,xndc,1.-tm);
    text.DrawLatex(xndc+0.005, 1-tm-0.01, Form("#alpha = %d #sigma", i ));
    ++i;
    k = FirstIndexLess(arr1, 1 - sigma(i));
  }

  canvas->SaveAs("confidence.png");
  canvas->SetLogy(0);
  canvas->SaveAs("confidenceLinear.png");
}
