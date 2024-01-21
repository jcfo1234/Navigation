//----------------------------------------------------------------------
// Distributions.h
// Probability distributions class
//----------------------------------------------------------------------

#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <map>
#include <math.h>
#include "..\..\DataTypeDefinitions.h"
using namespace std;

class Distributions
{
public:
   // Default Constructor
   Distributions();
   // Default destructor
   ~Distributions();
   // Get Normal distrubution abscise value
   DOUBLE GetNormalAbscise(DOUBLE dSignificance_);
   // Get the Chi square distribution abscise value
   DOUBLE GetChiSquareAbscise(ULONG ulDegofFreedom_, DOUBLE dSignificance_);
private:
   // Set Normal Distribution Table
   void SetNormalDistribution();
   // Set Chi square
   void SetChiSquare();
   // Set Chi square distribution 1 degree of freedom table
   void SetChiSquare1();
   // Set Chi square distribution 2 degrees of freedom table
   void SetChiSquare2();
   // Set Chi square distribution 3 degrees of freedom table
   void SetChiSquare3();
   // Set Chi square distribution 4 degrees of freedom table
   void SetChiSquare4();
   // Set Chi square distribution 5 degrees of freedom table
   void SetChiSquare5();
   // Set Chi square distribution 6 degrees of freedom table
   void SetChiSquare6();
   // Set Chi square distribution 7 degrees of freedom table
   void SetChiSquare7();
   // Set Chi square distribution 8 degrees of freedom table
   void SetChiSquare8();
   // Set Chi square distribution 9 degrees of freedom table
   void SetChiSquare9();
   // Set Chi square distribution 10 degrees of freedom table
   void SetChiSquare10();
   // Set Chi square distribution 11 degrees of freedom table
   void SetChiSquare11();
   // Set Chi square distribution 12 degrees of freedom table
   void SetChiSquare12();
   // Set Chi square distribution 13 degrees of freedom table
   void SetChiSquare13();
   // Set Chi square distribution 14 degrees of freedom table
   void SetChiSquare14();
   // Set Chi square distribution 15 degrees of freedom table
   void SetChiSquare15();
   // Set Chi square distribution 16 degrees of freedom table
   void SetChiSquare16();
   // Set Chi square distribution 17 degrees of freedom table
   void SetChiSquare17();
   // Set Chi square distribution 18 degrees of freedom table
   void SetChiSquare18();
   // Set Chi square distribution 19 degrees of freedom table
   void SetChiSquare19();
   // Set Chi square distribution 20 degrees of freedom table
   void SetChiSquare20();
   // Set Chi square distribution 21 degrees of freedom table
   void SetChiSquare21();
   // Set Chi square distribution 22 degrees of freedom table
   void SetChiSquare22();
   // Set Chi square distribution 23 degrees of freedom table
   void SetChiSquare23();
   // Set Chi square distribution 24 degrees of freedom table
   void SetChiSquare24();
   // Set Chi square distribution 25 degrees of freedom table
   void SetChiSquare25();
   // Set Chi square distribution 26 degrees of freedom table
   void SetChiSquare26();
   // Set Chi square distribution 27 degrees of freedom table
   void SetChiSquare27();
   // Set Chi square distribution 28 degrees of freedom table
   void SetChiSquare28();
   // Set Chi square distribution 29 degrees of freedom table
   void SetChiSquare29();
   // Set Chi square distribution 30 degrees of freedom table
   void SetChiSquare30();
   // Set Chi square distribution 40 degrees of freedom table
   void SetChiSquare40();
   // Set Chi square distribution 50 degrees of freedom table
   void SetChiSquare50();
   // Set Chi square distribution 60 degrees of freedom table
   void SetChiSquare60();
   // Set Chi square distribution 70 degrees of freedom table
   void SetChiSquare70();
   // Set Chi square distribution 80 degrees of freedom table
   void SetChiSquare80();
   // Set Chi square distribution 90 degrees of freedom table
   void SetChiSquare90();
   // Set Chi square distribution 100 degrees of freedom table
   void SetChiSquare100();

   // Normal probability distribution
   map<LONG, DOUBLE> mpNormalDistribution;
   // Chi square probability distribution
   map<ULONG, map<ULONG, DOUBLE>*> mpChiSquareDistribution;
   // Degrees of freedom Chi square distributions
   map<ULONG, DOUBLE> mpChiSquare1;
   map<ULONG, DOUBLE> mpChiSquare2;
   map<ULONG, DOUBLE> mpChiSquare3;
   map<ULONG, DOUBLE> mpChiSquare4;
   map<ULONG, DOUBLE> mpChiSquare5;
   map<ULONG, DOUBLE> mpChiSquare6;
   map<ULONG, DOUBLE> mpChiSquare7;
   map<ULONG, DOUBLE> mpChiSquare8;
   map<ULONG, DOUBLE> mpChiSquare9;
   map<ULONG, DOUBLE> mpChiSquare10;
   map<ULONG, DOUBLE> mpChiSquare11;
   map<ULONG, DOUBLE> mpChiSquare12;
   map<ULONG, DOUBLE> mpChiSquare13;
   map<ULONG, DOUBLE> mpChiSquare14;
   map<ULONG, DOUBLE> mpChiSquare15;
   map<ULONG, DOUBLE> mpChiSquare16;
   map<ULONG, DOUBLE> mpChiSquare17;
   map<ULONG, DOUBLE> mpChiSquare18;
   map<ULONG, DOUBLE> mpChiSquare19;
   map<ULONG, DOUBLE> mpChiSquare20;
   map<ULONG, DOUBLE> mpChiSquare21;
   map<ULONG, DOUBLE> mpChiSquare22;
   map<ULONG, DOUBLE> mpChiSquare23;
   map<ULONG, DOUBLE> mpChiSquare24;
   map<ULONG, DOUBLE> mpChiSquare25;
   map<ULONG, DOUBLE> mpChiSquare26;
   map<ULONG, DOUBLE> mpChiSquare27;
   map<ULONG, DOUBLE> mpChiSquare28;
   map<ULONG, DOUBLE> mpChiSquare29;
   map<ULONG, DOUBLE> mpChiSquare30;
   map<ULONG, DOUBLE> mpChiSquare40;
   map<ULONG, DOUBLE> mpChiSquare50;
   map<ULONG, DOUBLE> mpChiSquare60;
   map<ULONG, DOUBLE> mpChiSquare70;
   map<ULONG, DOUBLE> mpChiSquare80;
   map<ULONG, DOUBLE> mpChiSquare90;
   map<ULONG, DOUBLE> mpChiSquare100;
};

#endif