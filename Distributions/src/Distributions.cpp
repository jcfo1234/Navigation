//----------------------------------------------------------------------
// Distributions.cpp
// Source file of Distributions.h
// Source code of probability distribution tables
//----------------------------------------------------------------------
#include "..\inc\Distributions.h"

//------------------------------------------------------------------------------------------------
Distributions::Distributions()
{
   SetChiSquare1();
   SetChiSquare2();
   SetChiSquare3();
   SetChiSquare4();
   SetChiSquare5();
   SetChiSquare6();
   SetChiSquare7();
   SetChiSquare8();
   SetChiSquare9();
   SetChiSquare10();
   SetChiSquare11();
   SetChiSquare12();
   SetChiSquare13();
   SetChiSquare14();
   SetChiSquare15();
   SetChiSquare16();
   SetChiSquare17();
   SetChiSquare18();
   SetChiSquare19();
   SetChiSquare20();
   SetChiSquare21();
   SetChiSquare22();
   SetChiSquare23();
   SetChiSquare24();
   SetChiSquare25();
   SetChiSquare26();
   SetChiSquare27();
   SetChiSquare28();
   SetChiSquare29();
   SetChiSquare30();
   SetChiSquare40();
   SetChiSquare50();
   SetChiSquare60();
   SetChiSquare70();
   SetChiSquare80();
   SetChiSquare90();
   SetChiSquare100();
   SetChiSquare();
   SetNormalDistribution();
}

//------------------------------------------------------------------------------------------------
Distributions::~Distributions()
{
   mpChiSquareDistribution.clear();
   mpChiSquare1.clear();
   mpChiSquare2.clear();
   mpChiSquare3.clear();
   mpChiSquare4.clear();
   mpChiSquare5.clear();
   mpChiSquare6.clear();
   mpChiSquare7.clear();
   mpChiSquare8.clear();
   mpChiSquare9.clear();
   mpChiSquare10.clear();
   mpChiSquare11.clear();
   mpChiSquare12.clear();
   mpChiSquare13.clear();
   mpChiSquare14.clear();
   mpChiSquare15.clear();
   mpChiSquare16.clear();
   mpChiSquare17.clear();
   mpChiSquare18.clear();
   mpChiSquare19.clear();
   mpChiSquare20.clear();
   mpChiSquare21.clear();
   mpChiSquare22.clear();
   mpChiSquare23.clear();
   mpChiSquare24.clear();
   mpChiSquare25.clear();
   mpChiSquare26.clear();
   mpChiSquare27.clear();
   mpChiSquare28.clear();
   mpChiSquare29.clear();
   mpChiSquare30.clear();
   mpChiSquare40.clear();
   mpChiSquare50.clear();
   mpChiSquare60.clear();
   mpChiSquare70.clear();
   mpChiSquare80.clear();
   mpChiSquare90.clear();
   mpChiSquare100.clear();
   mpNormalDistribution.clear();
}

//------------------------------------------------------------------------------------------------
DOUBLE Distributions::GetNormalAbscise(DOUBLE dSignificance_)
{
   DOUBLE dAbscise = 0;
   // Search for significance level
   for (map<LONG, DOUBLE>::iterator itNormalIter = mpNormalDistribution.begin(); itNormalIter != mpNormalDistribution.end(); itNormalIter++)
   {
      // Found significance level
      if (abs(itNormalIter->second - dSignificance_) < 0.0001)
      {
         dAbscise = (DOUBLE)itNormalIter->first / 100.0;
         break;
      }
      // Overshoot of significance level
      else if (itNormalIter->second > dSignificance_)
      {
         map<LONG, DOUBLE>::iterator itTemp = itNormalIter--;
         dAbscise = ((DOUBLE)itNormalIter->first + (DOUBLE)itTemp->first) / 200;
         break;
      }
   }
   return dAbscise;
}

//------------------------------------------------------------------------------------------------
DOUBLE Distributions::GetChiSquareAbscise(ULONG ulDegofFreedom_, DOUBLE dSignificance_)
{
   DOUBLE dAbscise = 0;
   DOUBLE dValue = 0;
   for (map<ULONG, map<ULONG, DOUBLE>*>::iterator itChiSquare = mpChiSquareDistribution.begin(); itChiSquare != mpChiSquareDistribution.end(); itChiSquare++)
   {
      // Search first for the degrees of freedom
      if (itChiSquare->first == ulDegofFreedom_)
      {
         // Do a search on significance level once the degrees of freedom are found
         for (map<ULONG, DOUBLE>::iterator itChiSquare2 = itChiSquare->second->begin(); itChiSquare2 != itChiSquare->second->end(); itChiSquare2++)
         {
            dValue = (DOUBLE)itChiSquare2->first / 10000.0;
            if (abs(dValue - dSignificance_) < 0.0001)
            {
               dAbscise = itChiSquare2->second;
               break;
            }
         }
         break;
      }
   }

   return dAbscise;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare()
{
   mpChiSquareDistribution[1] = &mpChiSquare1;
   mpChiSquareDistribution[2] = &mpChiSquare2;
   mpChiSquareDistribution[3] = &mpChiSquare3;
   mpChiSquareDistribution[4] = &mpChiSquare4;
   mpChiSquareDistribution[5] = &mpChiSquare5;
   mpChiSquareDistribution[6] = &mpChiSquare6;
   mpChiSquareDistribution[7] = &mpChiSquare7;
   mpChiSquareDistribution[8] = &mpChiSquare8;
   mpChiSquareDistribution[9] = &mpChiSquare9;
   mpChiSquareDistribution[10] = &mpChiSquare10;
   mpChiSquareDistribution[11] = &mpChiSquare11;
   mpChiSquareDistribution[12] = &mpChiSquare12;
   mpChiSquareDistribution[13] = &mpChiSquare13;
   mpChiSquareDistribution[14] = &mpChiSquare14;
   mpChiSquareDistribution[15] = &mpChiSquare15;
   mpChiSquareDistribution[16] = &mpChiSquare16;
   mpChiSquareDistribution[17] = &mpChiSquare17;
   mpChiSquareDistribution[18] = &mpChiSquare18;
   mpChiSquareDistribution[19] = &mpChiSquare19;
   mpChiSquareDistribution[20] = &mpChiSquare20;
   mpChiSquareDistribution[21] = &mpChiSquare21;
   mpChiSquareDistribution[22] = &mpChiSquare22;
   mpChiSquareDistribution[23] = &mpChiSquare23;
   mpChiSquareDistribution[24] = &mpChiSquare24;
   mpChiSquareDistribution[25] = &mpChiSquare25;
   mpChiSquareDistribution[26] = &mpChiSquare26;
   mpChiSquareDistribution[27] = &mpChiSquare27;
   mpChiSquareDistribution[28] = &mpChiSquare28;
   mpChiSquareDistribution[29] = &mpChiSquare29;
   mpChiSquareDistribution[30] = &mpChiSquare30;
   mpChiSquareDistribution[40] = &mpChiSquare40;
   mpChiSquareDistribution[50] = &mpChiSquare50;
   mpChiSquareDistribution[60] = &mpChiSquare60;
   mpChiSquareDistribution[70] = &mpChiSquare70;
   mpChiSquareDistribution[80] = &mpChiSquare80;
   mpChiSquareDistribution[90] = &mpChiSquare90;
   mpChiSquareDistribution[100] = &mpChiSquare100;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare1()
{
   mpChiSquare1[5] = 12.116;
   mpChiSquare1[10] = 10.828;
   mpChiSquare1[50] = 7.879;
   mpChiSquare1[100] = 6.635;
   mpChiSquare1[250] = 5.024;
   mpChiSquare1[500] = 3.841;
   mpChiSquare1[750] = 3.170;
   mpChiSquare1[1000] = 2.706;
   mpChiSquare1[2000] = 1.642;
   mpChiSquare1[9000] = 0.016;
   mpChiSquare1[9500] = 0.004;
   mpChiSquare1[9750] = 0.001;
   mpChiSquare1[9990] = 0;
   mpChiSquare1[9995] = 0;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare2()
{
   mpChiSquare2[5] = 15.202;
   mpChiSquare2[10] = 13.816;
   mpChiSquare2[50] = 10.597;
   mpChiSquare2[100] = 9.210;
   mpChiSquare2[250] = 7.378;
   mpChiSquare2[500] = 5.991;
   mpChiSquare2[750] = 5.181;
   mpChiSquare2[1000] = 4.605;
   mpChiSquare2[2000] = 3.219;
   mpChiSquare2[9000] = 0.211;
   mpChiSquare2[9500] = 0.103;
   mpChiSquare2[9750] = 0.051;
   mpChiSquare2[9990] = 0.020;
   mpChiSquare2[9995] = 0.010;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare3()
{
   mpChiSquare3[5] = 17.731;
   mpChiSquare3[10] = 16.266;
   mpChiSquare3[50] = 12.838;
   mpChiSquare3[100] = 11.345;
   mpChiSquare3[250] = 9.348;
   mpChiSquare3[500] = 7.815;
   mpChiSquare3[750] = 6.905;
   mpChiSquare3[1000] = 6.251;
   mpChiSquare3[2000] = 4.642;
   mpChiSquare3[9000] = 0.584;
   mpChiSquare3[9500] = 0.352;
   mpChiSquare3[9750] = 0.216;
   mpChiSquare3[9990] = 0.115;
   mpChiSquare3[9995] = 0.072;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare4()
{
   mpChiSquare4[5] = 19.998;
   mpChiSquare4[10] = 18.467;
   mpChiSquare4[50] = 14.860;
   mpChiSquare4[100] = 13.277;
   mpChiSquare4[250] = 11.143;
   mpChiSquare4[500] = 9.488;
   mpChiSquare4[750] = 8.496;
   mpChiSquare4[1000] = 7.779;
   mpChiSquare4[2000] = 5.989;
   mpChiSquare4[9000] = 1.064;
   mpChiSquare4[9500] = 0.711;
   mpChiSquare4[9750] = 0.484;
   mpChiSquare4[9990] = 0.297;
   mpChiSquare4[9995] = 0.207;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare5()
{
   mpChiSquare5[5] = 22.106;
   mpChiSquare5[10] = 20.516;
   mpChiSquare5[50] = 16.750;
   mpChiSquare5[100] = 15.086;
   mpChiSquare5[250] = 12.833;
   mpChiSquare5[500] = 11.070;
   mpChiSquare5[750] = 10.008;
   mpChiSquare5[1000] = 9.236;
   mpChiSquare5[2000] = 7.289;
   mpChiSquare5[9000] = 1.610;
   mpChiSquare5[9500] = 1.145;
   mpChiSquare5[9750] = 0.831;
   mpChiSquare5[9990] = 0.554;
   mpChiSquare5[9995] = 0.412;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare6()
{
   mpChiSquare6[5] = 24.104;
   mpChiSquare6[10] = 22.458;
   mpChiSquare6[50] = 18.548;
   mpChiSquare6[100] = 16.812;
   mpChiSquare6[250] = 14.449;
   mpChiSquare6[500] = 12.592;
   mpChiSquare6[750] = 11.466;
   mpChiSquare6[1000] = 10.645;
   mpChiSquare6[2000] = 8.558;
   mpChiSquare6[9000] = 2.204;
   mpChiSquare6[9500] = 1.635;
   mpChiSquare6[9750] = 1.237;
   mpChiSquare6[9990] = 0.872;
   mpChiSquare6[9995] = 0.676;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare7()
{
   mpChiSquare7[5] = 26.019;
   mpChiSquare7[10] = 24.322;
   mpChiSquare7[50] = 20.278;
   mpChiSquare7[100] = 18.475;
   mpChiSquare7[250] = 16.013;
   mpChiSquare7[500] = 14.067;
   mpChiSquare7[750] = 12.883;
   mpChiSquare7[1000] = 12.017;
   mpChiSquare7[2000] = 9.803;
   mpChiSquare7[9000] = 2.833;
   mpChiSquare7[9500] = 2.167;
   mpChiSquare7[9750] = 1.690;
   mpChiSquare7[9990] = 1.239;
   mpChiSquare7[9995] = 0.989;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare8()
{
   mpChiSquare8[5] = 27.869;
   mpChiSquare8[10] = 26.125;
   mpChiSquare8[50] = 21.955;
   mpChiSquare8[100] = 20.090;
   mpChiSquare8[250] = 17.535;
   mpChiSquare8[500] = 15.507;
   mpChiSquare8[750] = 14.270;
   mpChiSquare8[1000] = 13.362;
   mpChiSquare8[2000] = 11.030;
   mpChiSquare8[9000] = 3.490;
   mpChiSquare8[9500] = 2.733;
   mpChiSquare8[9750] = 2.180;
   mpChiSquare8[9990] = 1.646;
   mpChiSquare8[9995] = 1.344;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare9()
{
   mpChiSquare9[5] = 29.667;
   mpChiSquare9[10] = 27.878;
   mpChiSquare9[50] = 23.589;
   mpChiSquare9[100] = 21.666;
   mpChiSquare9[250] = 19.023;
   mpChiSquare9[500] = 16.919;
   mpChiSquare9[750] = 15.631;
   mpChiSquare9[1000] = 14.684;
   mpChiSquare9[2000] = 12.242;
   mpChiSquare9[9000] = 4.168;
   mpChiSquare9[9500] = 3.325;
   mpChiSquare9[9750] = 2.700;
   mpChiSquare9[9990] = 2.088;
   mpChiSquare9[9995] = 1.735;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare10()
{
   mpChiSquare10[5] = 31.421;
   mpChiSquare10[10] = 29.589;
   mpChiSquare10[50] = 25.188;
   mpChiSquare10[100] = 23.209;
   mpChiSquare10[250] = 20.483;
   mpChiSquare10[500] = 18.307;
   mpChiSquare10[750] = 16.971;
   mpChiSquare10[1000] = 15.987;
   mpChiSquare10[2000] = 13.442;
   mpChiSquare10[9000] = 4.865;
   mpChiSquare10[9500] = 3.940;
   mpChiSquare10[9750] = 3.247;
   mpChiSquare10[9990] = 2.558;
   mpChiSquare10[9995] = 2.156;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare11()
{
   mpChiSquare11[5] = 33.138;
   mpChiSquare11[10] = 31.265;
   mpChiSquare11[50] = 26.757;
   mpChiSquare11[100] = 24.725;
   mpChiSquare11[250] = 21.920;
   mpChiSquare11[500] = 19.675;
   mpChiSquare11[750] = 18.294;
   mpChiSquare11[1000] = 17.275;
   mpChiSquare11[2000] = 14.631;
   mpChiSquare11[9000] = 5.578;
   mpChiSquare11[9500] = 4.575;
   mpChiSquare11[9750] = 3.816;
   mpChiSquare11[9990] = 3.053;
   mpChiSquare11[9995] = 2.603;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare12()
{
   mpChiSquare12[5] = 34.822;
   mpChiSquare12[10] = 32.910;
   mpChiSquare12[50] = 28.300;
   mpChiSquare12[100] = 26.217;
   mpChiSquare12[250] = 23.337;
   mpChiSquare12[500] = 21.026;
   mpChiSquare12[750] = 19.602;
   mpChiSquare12[1000] = 18.549;
   mpChiSquare12[2000] = 15.812;
   mpChiSquare12[9000] = 6.304;
   mpChiSquare12[9500] = 5.226;
   mpChiSquare12[9750] = 4.404;
   mpChiSquare12[9990] = 3.571;
   mpChiSquare12[9995] = 3.074;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare13()
{
   mpChiSquare13[5] = 36.479;
   mpChiSquare13[10] = 34.529;
   mpChiSquare13[50] = 29.820;
   mpChiSquare13[100] = 27.688;
   mpChiSquare13[250] = 24.736;
   mpChiSquare13[500] = 22.362;
   mpChiSquare13[750] = 20.987;
   mpChiSquare13[1000] = 19.812;
   mpChiSquare13[2000] = 16.985;
   mpChiSquare13[9000] = 7.042;
   mpChiSquare13[9500] = 5.892;
   mpChiSquare13[9750] = 5.009;
   mpChiSquare13[9990] = 4.107;
   mpChiSquare13[9995] = 3.565;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare14()
{
   mpChiSquare14[5] = 38.111;
   mpChiSquare14[10] = 36.124;
   mpChiSquare14[50] = 31.319;
   mpChiSquare14[100] = 29.141;
   mpChiSquare14[250] = 26.119;
   mpChiSquare14[500] = 23.685;
   mpChiSquare14[750] = 22.180;
   mpChiSquare14[1000] = 21.064;
   mpChiSquare14[2000] = 18.151;
   mpChiSquare14[9000] = 7.790;
   mpChiSquare14[9500] = 6.571;
   mpChiSquare14[9750] = 5.629;
   mpChiSquare14[9990] = 4.660;
   mpChiSquare14[9995] = 4.075;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare15()
{
   mpChiSquare15[5] = 39.720;
   mpChiSquare15[10] = 37.698;
   mpChiSquare15[50] = 32.801;
   mpChiSquare15[100] = 30.578;
   mpChiSquare15[250] = 27.488;
   mpChiSquare15[500] = 24.996;
   mpChiSquare15[750] = 23.452;
   mpChiSquare15[1000] = 22.307;
   mpChiSquare15[2000] = 19.311;
   mpChiSquare15[9000] = 8.547;
   mpChiSquare15[9500] = 7.261;
   mpChiSquare15[9750] = 6.242;
   mpChiSquare15[9990] = 5.229;
   mpChiSquare15[9995] = 4.601;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare16()
{
   mpChiSquare16[5] = 41.309;
   mpChiSquare16[10] = 39.253;
   mpChiSquare16[50] = 34.267;
   mpChiSquare16[100] = 32.000;
   mpChiSquare16[250] = 28.845;
   mpChiSquare16[500] = 26.296;
   mpChiSquare16[750] = 24.716;
   mpChiSquare16[1000] = 23.542;
   mpChiSquare16[2000] = 20.465;
   mpChiSquare16[9000] = 9.312;
   mpChiSquare16[9500] = 7.962;
   mpChiSquare16[9750] = 6.908;
   mpChiSquare16[9990] = 5.812;
   mpChiSquare16[9995] = 5.142;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare17()
{
   mpChiSquare17[5] = 42.881;
   mpChiSquare17[10] = 40.791;
   mpChiSquare17[50] = 35.719;
   mpChiSquare17[100] = 33.409;
   mpChiSquare17[250] = 30.191;
   mpChiSquare17[500] = 27.587;
   mpChiSquare17[750] = 25.970;
   mpChiSquare17[1000] = 24.769;
   mpChiSquare17[2000] = 21.615;
   mpChiSquare17[9000] = 10.085;
   mpChiSquare17[9500] = 8.672;
   mpChiSquare17[9750] = 7.564;
   mpChiSquare17[9990] = 6.408;
   mpChiSquare17[9995] = 5.697;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare18()
{
   mpChiSquare18[5] = 44.435;
   mpChiSquare18[10] = 42.314;
   mpChiSquare18[50] = 37.157;
   mpChiSquare18[100] = 34.805;
   mpChiSquare18[250] = 31.526;
   mpChiSquare18[500] = 28.869;
   mpChiSquare18[750] = 27.218;
   mpChiSquare18[1000] = 25.989;
   mpChiSquare18[2000] = 22.760;
   mpChiSquare18[9000] = 10.865;
   mpChiSquare18[9500] = 9.390;
   mpChiSquare18[9750] = 8.231;
   mpChiSquare18[9990] = 7.015;
   mpChiSquare18[9995] = 6.265;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare19()
{
   mpChiSquare19[5] = 45.974;
   mpChiSquare19[10] = 43.821;
   mpChiSquare19[50] = 38.582;
   mpChiSquare19[100] = 36.191;
   mpChiSquare19[250] = 32.852;
   mpChiSquare19[500] = 30.144;
   mpChiSquare19[750] = 28.458;
   mpChiSquare19[1000] = 27.204;
   mpChiSquare19[2000] = 23.900;
   mpChiSquare19[9000] = 11.651;
   mpChiSquare19[9500] = 10.117;
   mpChiSquare19[9750] = 8.907;
   mpChiSquare19[9990] = 7.633;
   mpChiSquare19[9995] = 6.844;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare20()
{
   mpChiSquare20[5] = 47.501;
   mpChiSquare20[10] = 45.315;
   mpChiSquare20[50] = 39.997;
   mpChiSquare20[100] = 37.566;
   mpChiSquare20[250] = 34.170;
   mpChiSquare20[500] = 31.410;
   mpChiSquare20[750] = 29.692;
   mpChiSquare20[1000] = 28.412;
   mpChiSquare20[2000] = 25.038;
   mpChiSquare20[9000] = 12.443;
   mpChiSquare20[9500] = 10.851;
   mpChiSquare20[9750] = 9.591;
   mpChiSquare20[9990] = 8.260;
   mpChiSquare20[9995] = 7.434;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare21()
{
   mpChiSquare21[5] = 49.013;
   mpChiSquare21[10] = 46.798;
   mpChiSquare21[50] = 41.401;
   mpChiSquare21[100] = 38.932;
   mpChiSquare21[250] = 35.479;
   mpChiSquare21[500] = 32.671;
   mpChiSquare21[750] = 30.920;
   mpChiSquare21[1000] = 29.615;
   mpChiSquare21[2000] = 26.171;
   mpChiSquare21[9000] = 13.240;
   mpChiSquare21[9500] = 11.591;
   mpChiSquare21[9750] = 10.283;
   mpChiSquare21[9990] = 8.897;
   mpChiSquare21[9995] = 8.034;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare22()
{
   mpChiSquare22[5] = 50.512;
   mpChiSquare22[10] = 48.269;
   mpChiSquare22[50] = 42.796;
   mpChiSquare22[100] = 40.289;
   mpChiSquare22[250] = 36.781;
   mpChiSquare22[500] = 33.924;
   mpChiSquare22[750] = 32.142;
   mpChiSquare22[1000] = 30.813;
   mpChiSquare22[2000] = 27.301;
   mpChiSquare22[9000] = 14.041;
   mpChiSquare22[9500] = 12.338;
   mpChiSquare22[9750] = 10.982;
   mpChiSquare22[9990] = 9.542;
   mpChiSquare22[9995] = 8.643;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare23()
{
   mpChiSquare23[5] = 52.002;
   mpChiSquare23[10] = 49.729;
   mpChiSquare23[50] = 44.182;
   mpChiSquare23[100] = 41.639;
   mpChiSquare23[250] = 38.076;
   mpChiSquare23[500] = 35.172;
   mpChiSquare23[750] = 33.360;
   mpChiSquare23[1000] = 32.007;
   mpChiSquare23[2000] = 28.429;
   mpChiSquare23[9000] = 14.848;
   mpChiSquare23[9500] = 13.091;
   mpChiSquare23[9750] = 11.689;
   mpChiSquare23[9990] = 10.196;
   mpChiSquare23[9995] = 9.260;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare24()
{
   mpChiSquare24[5] = 53.480;
   mpChiSquare24[10] = 51.180;
   mpChiSquare24[50] = 45.559;
   mpChiSquare24[100] = 42.980;
   mpChiSquare24[250] = 39.364;
   mpChiSquare24[500] = 36.415;
   mpChiSquare24[750] = 34.572;
   mpChiSquare24[1000] = 33.196;
   mpChiSquare24[2000] = 29.553;
   mpChiSquare24[9000] = 15.659;
   mpChiSquare24[9500] = 13.848;
   mpChiSquare24[9750] = 12.401;
   mpChiSquare24[9990] = 10.856;
   mpChiSquare24[9995] = 9.886;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare25()
{
   mpChiSquare25[5] = 54.950;
   mpChiSquare25[10] = 52.620;
   mpChiSquare25[50] = 46.928;
   mpChiSquare25[100] = 44.314;
   mpChiSquare25[250] = 40.646;
   mpChiSquare25[500] = 37.653;
   mpChiSquare25[750] = 35.780;
   mpChiSquare25[1000] = 34.382;
   mpChiSquare25[2000] = 30.675;
   mpChiSquare25[9000] = 16.473;
   mpChiSquare25[9500] = 14.611;
   mpChiSquare25[9750] = 13.120;
   mpChiSquare25[9990] = 11.524;
   mpChiSquare25[9995] = 10.520;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare26()
{
   mpChiSquare26[5] = 56.409;
   mpChiSquare26[10] = 54.053;
   mpChiSquare26[50] = 48.290;
   mpChiSquare26[100] = 45.642;
   mpChiSquare26[250] = 41.923;
   mpChiSquare26[500] = 38.885;
   mpChiSquare26[750] = 36.984;
   mpChiSquare26[1000] = 35.563;
   mpChiSquare26[2000] = 31.795;
   mpChiSquare26[9000] = 17.292;
   mpChiSquare26[9500] = 15.379;
   mpChiSquare26[9750] = 13.844;
   mpChiSquare26[9990] = 12.198;
   mpChiSquare26[9995] = 11.160;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare27()
{
   mpChiSquare27[5] = 57.860;
   mpChiSquare27[10] = 55.477;
   mpChiSquare27[50] = 49.645;
   mpChiSquare27[100] = 46.963;
   mpChiSquare27[250] = 43.195;
   mpChiSquare27[500] = 40.113;
   mpChiSquare27[750] = 38.184;
   mpChiSquare27[1000] = 36.741;
   mpChiSquare27[2000] = 32.912;
   mpChiSquare27[9000] = 18.114;
   mpChiSquare27[9500] = 16.151;
   mpChiSquare27[9750] = 14.573;
   mpChiSquare27[9990] = 12.879;
   mpChiSquare27[9995] = 11.808;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare28()
{
   mpChiSquare28[5] = 59.302;
   mpChiSquare28[10] = 56.894;
   mpChiSquare28[50] = 50.994;
   mpChiSquare28[100] = 48.278;
   mpChiSquare28[250] = 44.461;
   mpChiSquare28[500] = 41.337;
   mpChiSquare28[750] = 39.380;
   mpChiSquare28[1000] = 37.916;
   mpChiSquare28[2000] = 34.027;
   mpChiSquare28[9000] = 18.939;
   mpChiSquare28[9500] = 16.928;
   mpChiSquare28[9750] = 15.308;
   mpChiSquare28[9990] = 13.565;
   mpChiSquare28[9995] = 12.461;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare29()
{
   mpChiSquare29[5] = 60.738;
   mpChiSquare29[10] = 58.302;
   mpChiSquare29[50] = 52.336;
   mpChiSquare29[100] = 49.588;
   mpChiSquare29[250] = 45.722;
   mpChiSquare29[500] = 42.557;
   mpChiSquare29[750] = 40.573;
   mpChiSquare29[1000] = 39.087;
   mpChiSquare29[2000] = 35.139;
   mpChiSquare29[9000] = 19.768;
   mpChiSquare29[9500] = 17.708;
   mpChiSquare29[9750] = 16.047;
   mpChiSquare29[9990] = 14.256;
   mpChiSquare29[9995] = 13.121;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare30()
{
   mpChiSquare30[5] = 62.164;
   mpChiSquare30[10] = 59.704;
   mpChiSquare30[50] = 53.672;
   mpChiSquare30[100] = 50.892;
   mpChiSquare30[250] = 46.979;
   mpChiSquare30[500] = 43.773;
   mpChiSquare30[750] = 41.762;
   mpChiSquare30[1000] = 40.256;
   mpChiSquare30[2000] = 36.250;
   mpChiSquare30[9000] = 20.599;
   mpChiSquare30[9500] = 18.493;
   mpChiSquare30[9750] = 16.791;
   mpChiSquare30[9990] = 14.953;
   mpChiSquare30[9995] = 13.187;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare40()
{
   mpChiSquare40[5] = 76.097;
   mpChiSquare40[10] = 73.403;
   mpChiSquare40[50] = 66.766;
   mpChiSquare40[100] = 63.691;
   mpChiSquare40[250] = 59.342;
   mpChiSquare40[500] = 55.759;
   mpChiSquare40[750] = 53.501;
   mpChiSquare40[1000] = 51.805;
   mpChiSquare40[2000] = 47.269;
   mpChiSquare40[9000] = 29.051;
   mpChiSquare40[9500] = 26.509;
   mpChiSquare40[9750] = 24.433;
   mpChiSquare40[9990] = 22.164;
   mpChiSquare40[9995] = 20.707;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare50()
{
   mpChiSquare50[5] = 89.564;
   mpChiSquare50[10] = 86.662;
   mpChiSquare50[50] = 79.490;
   mpChiSquare50[100] = 76.154;
   mpChiSquare50[250] = 71.420;
   mpChiSquare50[500] = 67.505;
   mpChiSquare50[750] = 65.030;
   mpChiSquare50[1000] = 63.167;
   mpChiSquare50[2000] = 58.164;
   mpChiSquare50[9000] = 37.689;
   mpChiSquare50[9500] = 34.764;
   mpChiSquare50[9750] = 32.357;
   mpChiSquare50[9990] = 29.707;
   mpChiSquare50[9995] = 27.991;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare60()
{
   mpChiSquare60[5] = 102.698;
   mpChiSquare60[10] = 99.609;
   mpChiSquare60[50] = 91.952;
   mpChiSquare60[100] = 88.380;
   mpChiSquare60[250] = 83.298;
   mpChiSquare60[500] = 79.082;
   mpChiSquare60[750] = 76.411;
   mpChiSquare60[1000] = 74.397;
   mpChiSquare60[2000] = 68.972;
   mpChiSquare60[9000] = 46.459;
   mpChiSquare60[9500] = 43.188;
   mpChiSquare60[9750] = 40.482;
   mpChiSquare60[9990] = 37.485;
   mpChiSquare60[9995] = 35.534;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare70()
{
   mpChiSquare70[5] = 115.582;
   mpChiSquare70[10] = 112.319;
   mpChiSquare70[50] = 104.215;
   mpChiSquare70[100] = 100.425;
   mpChiSquare70[250] = 95.023;
   mpChiSquare70[500] = 90.531;
   mpChiSquare70[750] = 87.680;
   mpChiSquare70[1000] = 85.527;
   mpChiSquare70[2000] = 79.715;
   mpChiSquare70[9000] = 55.329;
   mpChiSquare70[9500] = 51.739;
   mpChiSquare70[9750] = 48.758;
   mpChiSquare70[9990] = 45.442;
   mpChiSquare70[9995] = 43.275;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare80()
{
   mpChiSquare80[5] = 128.267;
   mpChiSquare80[10] = 124.842;
   mpChiSquare80[50] = 116.321;
   mpChiSquare80[100] = 112.329;
   mpChiSquare80[250] = 106.629;
   mpChiSquare80[500] = 101.880;
   mpChiSquare80[750] = 98.861;
   mpChiSquare80[1000] = 96.578;
   mpChiSquare80[2000] = 90.405;
   mpChiSquare80[9000] = 64.278;
   mpChiSquare80[9500] = 60.391;
   mpChiSquare80[9750] = 57.153;
   mpChiSquare80[9990] = 53.540;
   mpChiSquare80[9995] = 51.172;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare90()
{
   mpChiSquare90[5] = 140.789;
   mpChiSquare90[10] = 137.211;
   mpChiSquare90[50] = 128.300;
   mpChiSquare90[100] = 124.117;
   mpChiSquare90[250] = 118.136;
   mpChiSquare90[500] = 113.145;
   mpChiSquare90[750] = 109.969;
   mpChiSquare90[1000] = 107.565;
   mpChiSquare90[2000] = 101.054;
   mpChiSquare90[9000] = 73.291;
   mpChiSquare90[9500] = 69.126;
   mpChiSquare90[9750] = 65.647;
   mpChiSquare90[9990] = 61.754;
   mpChiSquare90[9995] = 59.196;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetChiSquare100()
{
   mpChiSquare100[5] = 153.174;
   mpChiSquare100[10] = 149.452;
   mpChiSquare100[50] = 140.170;
   mpChiSquare100[100] = 135.807;
   mpChiSquare100[250] = 129.561;
   mpChiSquare100[500] = 124.342;
   mpChiSquare100[750] = 121.017;
   mpChiSquare100[1000] = 118.498;
   mpChiSquare100[2000] = 111.667;
   mpChiSquare100[9000] = 82.358;
   mpChiSquare100[9500] = 77.929;
   mpChiSquare100[9750] = 74.222;
   mpChiSquare100[9990] = 70.065;
   mpChiSquare100[9995] = 67.328;
}

//------------------------------------------------------------------------------------------------
void Distributions::SetNormalDistribution()
{
   mpNormalDistribution[-399] = 0.00003;
   mpNormalDistribution[-398] = 0.00003;
   mpNormalDistribution[-397] = 0.00004;
   mpNormalDistribution[-396] = 0.00004;
   mpNormalDistribution[-395] = 0.00004;
   mpNormalDistribution[-394] = 0.00004;
   mpNormalDistribution[-393] = 0.00004;
   mpNormalDistribution[-392] = 0.00004;
   mpNormalDistribution[-391] = 0.00005;
   mpNormalDistribution[-390] = 0.00005;
   mpNormalDistribution[-389] = 0.00005;
   mpNormalDistribution[-388] = 0.00005;
   mpNormalDistribution[-387] = 0.00005;
   mpNormalDistribution[-386] = 0.00006;
   mpNormalDistribution[-385] = 0.00006;
   mpNormalDistribution[-384] = 0.00006;
   mpNormalDistribution[-383] = 0.00006;
   mpNormalDistribution[-382] = 0.00007;
   mpNormalDistribution[-381] = 0.00007;
   mpNormalDistribution[-380] = 0.00007;
   mpNormalDistribution[-379] = 0.00008;
   mpNormalDistribution[-378] = 0.00008;
   mpNormalDistribution[-377] = 0.00008;
   mpNormalDistribution[-376] = 0.00008;
   mpNormalDistribution[-375] = 0.00009;
   mpNormalDistribution[-374] = 0.00009;
   mpNormalDistribution[-373] = 0.00010;
   mpNormalDistribution[-372] = 0.00010;
   mpNormalDistribution[-371] = 0.00010;
   mpNormalDistribution[-370] = 0.00011;
   mpNormalDistribution[-369] = 0.00011;
   mpNormalDistribution[-368] = 0.00012;
   mpNormalDistribution[-367] = 0.00012;
   mpNormalDistribution[-366] = 0.00013;
   mpNormalDistribution[-365] = 0.00013;
   mpNormalDistribution[-364] = 0.00014;
   mpNormalDistribution[-363] = 0.00014;
   mpNormalDistribution[-362] = 0.00015;
   mpNormalDistribution[-361] = 0.00015;
   mpNormalDistribution[-360] = 0.00016;
   mpNormalDistribution[-359] = 0.00017;
   mpNormalDistribution[-358] = 0.00017;
   mpNormalDistribution[-357] = 0.00018;
   mpNormalDistribution[-356] = 0.00019;
   mpNormalDistribution[-355] = 0.00019;
   mpNormalDistribution[-354] = 0.00020;
   mpNormalDistribution[-353] = 0.00021;
   mpNormalDistribution[-352] = 0.00022;
   mpNormalDistribution[-351] = 0.00022;
   mpNormalDistribution[-350] = 0.00023;
   mpNormalDistribution[-349] = 0.00024;
   mpNormalDistribution[-348] = 0.00025;
   mpNormalDistribution[-347] = 0.00026;
   mpNormalDistribution[-346] = 0.00027;
   mpNormalDistribution[-345] = 0.00028;
   mpNormalDistribution[-344] = 0.00029;
   mpNormalDistribution[-343] = 0.00030;
   mpNormalDistribution[-342] = 0.00031;
   mpNormalDistribution[-341] = 0.00032;
   mpNormalDistribution[-340] = 0.00034;
   mpNormalDistribution[-339] = 0.00035;
   mpNormalDistribution[-338] = 0.00036;
   mpNormalDistribution[-337] = 0.00038;
   mpNormalDistribution[-336] = 0.00039;
   mpNormalDistribution[-335] = 0.00040;
   mpNormalDistribution[-334] = 0.00042;
   mpNormalDistribution[-333] = 0.00043;
   mpNormalDistribution[-332] = 0.00045;
   mpNormalDistribution[-331] = 0.00047;
   mpNormalDistribution[-330] = 0.00048;
   mpNormalDistribution[-329] = 0.00050;
   mpNormalDistribution[-328] = 0.00052;
   mpNormalDistribution[-327] = 0.00054;
   mpNormalDistribution[-326] = 0.00056;
   mpNormalDistribution[-325] = 0.00058;
   mpNormalDistribution[-324] = 0.00060;
   mpNormalDistribution[-323] = 0.00062;
   mpNormalDistribution[-322] = 0.00064;
   mpNormalDistribution[-321] = 0.00066;
   mpNormalDistribution[-320] = 0.00069;
   mpNormalDistribution[-319] = 0.00071;
   mpNormalDistribution[-318] = 0.00074;
   mpNormalDistribution[-317] = 0.00076;
   mpNormalDistribution[-316] = 0.00079;
   mpNormalDistribution[-315] = 0.00082;
   mpNormalDistribution[-314] = 0.00084;
   mpNormalDistribution[-313] = 0.00097;
   mpNormalDistribution[-312] = 0.00090;
   mpNormalDistribution[-311] = 0.00094;
   mpNormalDistribution[-310] = 0.00097;
   mpNormalDistribution[-309] = 0.00100;
   mpNormalDistribution[-308] = 0.00104;
   mpNormalDistribution[-307] = 0.00107;
   mpNormalDistribution[-306] = 0.00111;
   mpNormalDistribution[-305] = 0.00114;
   mpNormalDistribution[-304] = 0.00118;
   mpNormalDistribution[-303] = 0.00122;
   mpNormalDistribution[-302] = 0.00126;
   mpNormalDistribution[-301] = 0.00131;
   mpNormalDistribution[-300] = 0.00135;
   mpNormalDistribution[-299] = 0.00139;
   mpNormalDistribution[-298] = 0.00144;
   mpNormalDistribution[-297] = 0.00149;
   mpNormalDistribution[-296] = 0.00154;
   mpNormalDistribution[-295] = 0.00159;
   mpNormalDistribution[-294] = 0.00164;
   mpNormalDistribution[-293] = 0.00169;
   mpNormalDistribution[-292] = 0.00175;
   mpNormalDistribution[-291] = 0.00181;
   mpNormalDistribution[-290] = 0.00187;
   mpNormalDistribution[-289] = 0.00193;
   mpNormalDistribution[-288] = 0.00199;
   mpNormalDistribution[-287] = 0.00205;
   mpNormalDistribution[-286] = 0.00212;
   mpNormalDistribution[-285] = 0.00219;
   mpNormalDistribution[-284] = 0.00226;
   mpNormalDistribution[-283] = 0.00233;
   mpNormalDistribution[-282] = 0.00240;
   mpNormalDistribution[-281] = 0.00248;
   mpNormalDistribution[-280] = 0.00256;
   mpNormalDistribution[-279] = 0.00264;
   mpNormalDistribution[-278] = 0.00272;
   mpNormalDistribution[-277] = 0.00280;
   mpNormalDistribution[-276] = 0.00289;
   mpNormalDistribution[-275] = 0.00298;
   mpNormalDistribution[-274] = 0.00307;
   mpNormalDistribution[-273] = 0.00317;
   mpNormalDistribution[-272] = 0.00326;
   mpNormalDistribution[-271] = 0.00336;
   mpNormalDistribution[-270] = 0.00347;
   mpNormalDistribution[-269] = 0.00357;
   mpNormalDistribution[-268] = 0.00368;
   mpNormalDistribution[-267] = 0.00379;
   mpNormalDistribution[-266] = 0.00391;
   mpNormalDistribution[-265] = 0.00402;
   mpNormalDistribution[-264] = 0.00415;
   mpNormalDistribution[-263] = 0.00427;
   mpNormalDistribution[-262] = 0.00440;
   mpNormalDistribution[-261] = 0.00453;
   mpNormalDistribution[-260] = 0.00466;
   mpNormalDistribution[-259] = 0.00480;
   mpNormalDistribution[-258] = 0.00494;
   mpNormalDistribution[-257] = 0.00508;
   mpNormalDistribution[-256] = 0.00523;
   mpNormalDistribution[-255] = 0.00539;
   mpNormalDistribution[-254] = 0.00554;
   mpNormalDistribution[-253] = 0.00570;
   mpNormalDistribution[-252] = 0.00587;
   mpNormalDistribution[-251] = 0.00604;
   mpNormalDistribution[-250] = 0.00621;
   mpNormalDistribution[-249] = 0.00639;
   mpNormalDistribution[-248] = 0.00657;
   mpNormalDistribution[-247] = 0.00676;
   mpNormalDistribution[-246] = 0.00695;
   mpNormalDistribution[-245] = 0.00714;
   mpNormalDistribution[-244] = 0.00734;
   mpNormalDistribution[-243] = 0.00755;
   mpNormalDistribution[-242] = 0.00776;
   mpNormalDistribution[-241] = 0.00798;
   mpNormalDistribution[-240] = 0.00820;
   mpNormalDistribution[-239] = 0.00842;
   mpNormalDistribution[-238] = 0.00866;
   mpNormalDistribution[-237] = 0.00889;
   mpNormalDistribution[-236] = 0.00914;
   mpNormalDistribution[-235] = 0.00939;
   mpNormalDistribution[-234] = 0.00964;
   mpNormalDistribution[-233] = 0.00990;
   mpNormalDistribution[-232] = 0.01017;
   mpNormalDistribution[-231] = 0.01044;
   mpNormalDistribution[-230] = 0.01072;
   mpNormalDistribution[-229] = 0.01101;
   mpNormalDistribution[-228] = 0.01130;
   mpNormalDistribution[-227] = 0.01160;
   mpNormalDistribution[-226] = 0.01191;
   mpNormalDistribution[-225] = 0.01222;
   mpNormalDistribution[-224] = 0.01255;
   mpNormalDistribution[-223] = 0.01287;
   mpNormalDistribution[-222] = 0.01321;
   mpNormalDistribution[-221] = 0.01355;
   mpNormalDistribution[-220] = 0.01390;
   mpNormalDistribution[-219] = 0.01426;
   mpNormalDistribution[-218] = 0.01463;
   mpNormalDistribution[-217] = 0.01500;
   mpNormalDistribution[-216] = 0.01539;
   mpNormalDistribution[-215] = 0.01578;
   mpNormalDistribution[-214] = 0.01618;
   mpNormalDistribution[-213] = 0.01659;
   mpNormalDistribution[-212] = 0.01700;
   mpNormalDistribution[-211] = 0.01743;
   mpNormalDistribution[-210] = 0.01786;
   mpNormalDistribution[-209] = 0.01831;
   mpNormalDistribution[-208] = 0.01876;
   mpNormalDistribution[-207] = 0.01923;
   mpNormalDistribution[-206] = 0.01970;
   mpNormalDistribution[-205] = 0.02018;
   mpNormalDistribution[-204] = 0.02068;
   mpNormalDistribution[-203] = 0.02118;
   mpNormalDistribution[-202] = 0.02169;
   mpNormalDistribution[-201] = 0.02222;
   mpNormalDistribution[-200] = 0.02275;
   mpNormalDistribution[-199] = 0.02330;
   mpNormalDistribution[-198] = 0.02385;
   mpNormalDistribution[-197] = 0.02442;
   mpNormalDistribution[-196] = 0.02500;
   mpNormalDistribution[-195] = 0.02559;
   mpNormalDistribution[-194] = 0.02619;
   mpNormalDistribution[-193] = 0.02680;
   mpNormalDistribution[-192] = 0.02743;
   mpNormalDistribution[-191] = 0.02807;
   mpNormalDistribution[-190] = 0.02872;
   mpNormalDistribution[-189] = 0.02938;
   mpNormalDistribution[-188] = 0.03005;
   mpNormalDistribution[-187] = 0.03074;
   mpNormalDistribution[-186] = 0.03144;
   mpNormalDistribution[-185] = 0.03216;
   mpNormalDistribution[-184] = 0.03288;
   mpNormalDistribution[-183] = 0.03362;
   mpNormalDistribution[-182] = 0.03438;
   mpNormalDistribution[-181] = 0.03515;
   mpNormalDistribution[-180] = 0.03593;
   mpNormalDistribution[-179] = 0.03673;
   mpNormalDistribution[-178] = 0.03754;
   mpNormalDistribution[-177] = 0.03836;
   mpNormalDistribution[-176] = 0.03920;
   mpNormalDistribution[-175] = 0.04006;
   mpNormalDistribution[-174] = 0.04093;
   mpNormalDistribution[-173] = 0.04182;
   mpNormalDistribution[-172] = 0.04272;
   mpNormalDistribution[-171] = 0.04363;
   mpNormalDistribution[-170] = 0.04457;
   mpNormalDistribution[-169] = 0.04551;
   mpNormalDistribution[-168] = 0.04648;
   mpNormalDistribution[-167] = 0.04746;
   mpNormalDistribution[-166] = 0.04846;
   mpNormalDistribution[-165] = 0.04947;
   mpNormalDistribution[-164] = 0.05050;
   mpNormalDistribution[-163] = 0.05155;
   mpNormalDistribution[-162] = 0.05262;
   mpNormalDistribution[-161] = 0.05370;
   mpNormalDistribution[-160] = 0.05480;
   mpNormalDistribution[-159] = 0.05592;
   mpNormalDistribution[-158] = 0.05705;
   mpNormalDistribution[-157] = 0.05821;
   mpNormalDistribution[-156] = 0.05938;
   mpNormalDistribution[-155] = 0.06057;
   mpNormalDistribution[-154] = 0.06178;
   mpNormalDistribution[-153] = 0.06301;
   mpNormalDistribution[-152] = 0.06426;
   mpNormalDistribution[-151] = 0.06552;
   mpNormalDistribution[-150] = 0.06681;
   mpNormalDistribution[-149] = 0.06811;
   mpNormalDistribution[-148] = 0.06944;
   mpNormalDistribution[-147] = 0.07078;
   mpNormalDistribution[-146] = 0.07215;
   mpNormalDistribution[-145] = 0.07353;
   mpNormalDistribution[-144] = 0.07493;
   mpNormalDistribution[-143] = 0.07636;
   mpNormalDistribution[-142] = 0.07780;
   mpNormalDistribution[-141] = 0.07927;
   mpNormalDistribution[-140] = 0.08076;
   mpNormalDistribution[-139] = 0.08226;
   mpNormalDistribution[-138] = 0.08379;
   mpNormalDistribution[-137] = 0.08534;
   mpNormalDistribution[-136] = 0.08691;
   mpNormalDistribution[-135] = 0.08851;
   mpNormalDistribution[-134] = 0.09012;
   mpNormalDistribution[-133] = 0.09176;
   mpNormalDistribution[-132] = 0.09342;
   mpNormalDistribution[-131] = 0.09510;
   mpNormalDistribution[-130] = 0.09680;
   mpNormalDistribution[-129] = 0.09853;
   mpNormalDistribution[-128] = 0.10027;
   mpNormalDistribution[-127] = 0.10204;
   mpNormalDistribution[-126] = 0.10383;
   mpNormalDistribution[-125] = 0.10565;
   mpNormalDistribution[-124] = 0.10749;
   mpNormalDistribution[-123] = 0.10935;
   mpNormalDistribution[-122] = 0.11123;
   mpNormalDistribution[-121] = 0.11314;
   mpNormalDistribution[-120] = 0.11507;
   mpNormalDistribution[-119] = 0.11702;
   mpNormalDistribution[-118] = 0.11900;
   mpNormalDistribution[-117] = 0.12100;
   mpNormalDistribution[-116] = 0.12302;
   mpNormalDistribution[-115] = 0.12507;
   mpNormalDistribution[-114] = 0.12714;
   mpNormalDistribution[-113] = 0.12924;
   mpNormalDistribution[-112] = 0.13136;
   mpNormalDistribution[-111] = 0.13350;
   mpNormalDistribution[-110] = 0.13567;
   mpNormalDistribution[-109] = 0.13786;
   mpNormalDistribution[-108] = 0.14007;
   mpNormalDistribution[-107] = 0.14231;
   mpNormalDistribution[-106] = 0.14457;
   mpNormalDistribution[-105] = 0.14686;
   mpNormalDistribution[-104] = 0.14917;
   mpNormalDistribution[-103] = 0.15151;
   mpNormalDistribution[-102] = 0.15386;
   mpNormalDistribution[-101] = 0.15625;
   mpNormalDistribution[-100] = 0.15866;
   mpNormalDistribution[-99] = 0.16109;
   mpNormalDistribution[-98] = 0.16354;
   mpNormalDistribution[-97] = 0.16602;
   mpNormalDistribution[-96] = 0.16853;
   mpNormalDistribution[-95] = 0.17106;
   mpNormalDistribution[-94] = 0.17361;
   mpNormalDistribution[-93] = 0.17619;
   mpNormalDistribution[-92] = 0.17879;
   mpNormalDistribution[-91] = 0.18141;
   mpNormalDistribution[-90] = 0.18406;
   mpNormalDistribution[-89] = 0.18673;
   mpNormalDistribution[-88] = 0.18943;
   mpNormalDistribution[-87] = 0.19215;
   mpNormalDistribution[-86] = 0.19489;
   mpNormalDistribution[-85] = 0.19766;
   mpNormalDistribution[-84] = 0.20045;
   mpNormalDistribution[-83] = 0.20327;
   mpNormalDistribution[-82] = 0.20611;
   mpNormalDistribution[-81] = 0.20897;
   mpNormalDistribution[-80] = 0.21186;
   mpNormalDistribution[-79] = 0.21476;
   mpNormalDistribution[-78] = 0.21770;
   mpNormalDistribution[-77] = 0.22065;
   mpNormalDistribution[-76] = 0.22363;
   mpNormalDistribution[-75] = 0.22663;
   mpNormalDistribution[-74] = 0.22965;
   mpNormalDistribution[-73] = 0.23270;
   mpNormalDistribution[-72] = 0.23576;
   mpNormalDistribution[-71] = 0.23885;
   mpNormalDistribution[-70] = 0.24196;
   mpNormalDistribution[-69] = 0.24510;
   mpNormalDistribution[-68] = 0.24825;
   mpNormalDistribution[-67] = 0.25143;
   mpNormalDistribution[-66] = 0.25463;
   mpNormalDistribution[-65] = 0.25785;
   mpNormalDistribution[-64] = 0.26109;
   mpNormalDistribution[-63] = 0.26435;
   mpNormalDistribution[-62] = 0.26763;
   mpNormalDistribution[-61] = 0.27093;
   mpNormalDistribution[-60] = 0.27425;
   mpNormalDistribution[-59] = 0.27760;
   mpNormalDistribution[-58] = 0.28096;
   mpNormalDistribution[-57] = 0.28434;
   mpNormalDistribution[-56] = 0.28774;
   mpNormalDistribution[-55] = 0.29116;
   mpNormalDistribution[-54] = 0.29460;
   mpNormalDistribution[-53] = 0.29806;
   mpNormalDistribution[-52] = 0.30153;
   mpNormalDistribution[-51] = 0.30503;
   mpNormalDistribution[-50] = 0.30854;
   mpNormalDistribution[-49] = 0.31207;
   mpNormalDistribution[-48] = 0.31561;
   mpNormalDistribution[-47] = 0.31918;
   mpNormalDistribution[-46] = 0.32276;
   mpNormalDistribution[-45] = 0.32636;
   mpNormalDistribution[-44] = 0.32997;
   mpNormalDistribution[-43] = 0.33360;
   mpNormalDistribution[-42] = 0.33724;
   mpNormalDistribution[-41] = 0.34090;
   mpNormalDistribution[-40] = 0.34458;
   mpNormalDistribution[-39] = 0.34827;
   mpNormalDistribution[-38] = 0.35197;
   mpNormalDistribution[-37] = 0.35569;
   mpNormalDistribution[-36] = 0.35942;
   mpNormalDistribution[-35] = 0.36317;
   mpNormalDistribution[-34] = 0.36693;
   mpNormalDistribution[-33] = 0.37070;
   mpNormalDistribution[-32] = 0.37448;
   mpNormalDistribution[-31] = 0.37828;
   mpNormalDistribution[-30] = 0.38209;
   mpNormalDistribution[-29] = 0.38591;
   mpNormalDistribution[-28] = 0.38974;
   mpNormalDistribution[-27] = 0.39358;
   mpNormalDistribution[-26] = 0.39743;
   mpNormalDistribution[-25] = 0.40129;
   mpNormalDistribution[-24] = 0.40517;
   mpNormalDistribution[-23] = 0.40905;
   mpNormalDistribution[-22] = 0.41294;
   mpNormalDistribution[-21] = 0.41683;
   mpNormalDistribution[-20] = 0.42074;
   mpNormalDistribution[-19] = 0.42465;
   mpNormalDistribution[-18] = 0.42858;
   mpNormalDistribution[-17] = 0.43251;
   mpNormalDistribution[-16] = 0.43644;
   mpNormalDistribution[-15] = 0.44038;
   mpNormalDistribution[-14] = 0.44433;
   mpNormalDistribution[-13] = 0.44828;
   mpNormalDistribution[-12] = 0.45224;
   mpNormalDistribution[-11] = 0.45620;
   mpNormalDistribution[-10] = 0.46017;
   mpNormalDistribution[-9] = 0.46414;
   mpNormalDistribution[-8] = 0.46812;
   mpNormalDistribution[-7] = 0.47210;
   mpNormalDistribution[-6] = 0.47608;
   mpNormalDistribution[-5] = 0.48006;
   mpNormalDistribution[-4] = 0.48405;
   mpNormalDistribution[-3] = 0.48803;
   mpNormalDistribution[-2] = 0.49202;
   mpNormalDistribution[-1] = 0.49601;
   mpNormalDistribution[0] = 0.50000;
   mpNormalDistribution[1] = 0.50399;
   mpNormalDistribution[2] = 0.50798;
   mpNormalDistribution[3] = 0.51197;
   mpNormalDistribution[4] = 0.51595;
   mpNormalDistribution[5] = 0.51994;
   mpNormalDistribution[6] = 0.52392;
   mpNormalDistribution[7] = 0.52790;
   mpNormalDistribution[8] = 0.53188;
   mpNormalDistribution[9] = 0.53586;
   mpNormalDistribution[10] = 0.53983;
   mpNormalDistribution[11] = 0.54380;
   mpNormalDistribution[12] = 0.54776;
   mpNormalDistribution[13] = 0.55172;
   mpNormalDistribution[14] = 0.55567;
   mpNormalDistribution[15] = 0.55962;
   mpNormalDistribution[16] = 0.56356;
   mpNormalDistribution[17] = 0.56749;
   mpNormalDistribution[18] = 0.57142;
   mpNormalDistribution[19] = 0.57535;
   mpNormalDistribution[20] = 0.57926;
   mpNormalDistribution[21] = 0.58317;
   mpNormalDistribution[22] = 0.58706;
   mpNormalDistribution[23] = 0.59095;
   mpNormalDistribution[24] = 0.59483;
   mpNormalDistribution[25] = 0.59871;
   mpNormalDistribution[26] = 0.60257;
   mpNormalDistribution[27] = 0.60642;
   mpNormalDistribution[28] = 0.61026;
   mpNormalDistribution[29] = 0.61409;
   mpNormalDistribution[30] = 0.61791;
   mpNormalDistribution[31] = 0.62172;
   mpNormalDistribution[32] = 0.62552;
   mpNormalDistribution[33] = 0.62930;
   mpNormalDistribution[34] = 0.63307;
   mpNormalDistribution[35] = 0.63683;
   mpNormalDistribution[36] = 0.64058;
   mpNormalDistribution[37] = 0.64431;
   mpNormalDistribution[38] = 0.64803;
   mpNormalDistribution[39] = 0.65173;
   mpNormalDistribution[40] = 0.65542;
   mpNormalDistribution[41] = 0.65910;
   mpNormalDistribution[42] = 0.66276;
   mpNormalDistribution[43] = 0.66640;
   mpNormalDistribution[44] = 0.67003;
   mpNormalDistribution[45] = 0.67364;
   mpNormalDistribution[46] = 0.67724;
   mpNormalDistribution[47] = 0.68802;
   mpNormalDistribution[48] = 0.68439;
   mpNormalDistribution[49] = 0.68793;
   mpNormalDistribution[50] = 0.69146;
   mpNormalDistribution[51] = 0.69497;
   mpNormalDistribution[52] = 0.69847;
   mpNormalDistribution[53] = 0.70194;
   mpNormalDistribution[54] = 0.70540;
   mpNormalDistribution[55] = 0.70884;
   mpNormalDistribution[56] = 0.71226;
   mpNormalDistribution[57] = 0.71566;
   mpNormalDistribution[58] = 0.71904;
   mpNormalDistribution[59] = 0.72240;
   mpNormalDistribution[60] = 0.72575;
   mpNormalDistribution[61] = 0.72907;
   mpNormalDistribution[62] = 0.73237;
   mpNormalDistribution[63] = 0.73565;
   mpNormalDistribution[64] = 0.73891;
   mpNormalDistribution[65] = 0.74215;
   mpNormalDistribution[66] = 0.74537;
   mpNormalDistribution[67] = 0.74857;
   mpNormalDistribution[68] = 0.75175;
   mpNormalDistribution[69] = 0.75490;
   mpNormalDistribution[70] = 0.75804;
   mpNormalDistribution[71] = 0.76115;
   mpNormalDistribution[72] = 0.76424;
   mpNormalDistribution[73] = 0.76730;
   mpNormalDistribution[74] = 0.77035;
   mpNormalDistribution[75] = 0.77337;
   mpNormalDistribution[76] = 0.77637;
   mpNormalDistribution[77] = 0.77935;
   mpNormalDistribution[78] = 0.78230;
   mpNormalDistribution[79] = 0.78524;
   mpNormalDistribution[80] = 0.78814;
   mpNormalDistribution[81] = 0.79103;
   mpNormalDistribution[82] = 0.79389;
   mpNormalDistribution[83] = 0.79673;
   mpNormalDistribution[84] = 0.79955;
   mpNormalDistribution[85] = 0.80234;
   mpNormalDistribution[86] = 0.80511;
   mpNormalDistribution[87] = 0.80785;
   mpNormalDistribution[88] = 0.81057;
   mpNormalDistribution[89] = 0.81327;
   mpNormalDistribution[90] = 0.81594;
   mpNormalDistribution[91] = 0.81859;
   mpNormalDistribution[92] = 0.82121;
   mpNormalDistribution[93] = 0.82381;
   mpNormalDistribution[94] = 0.82639;
   mpNormalDistribution[95] = 0.82894;
   mpNormalDistribution[96] = 0.83147;
   mpNormalDistribution[97] = 0.83398;
   mpNormalDistribution[98] = 0.83646;
   mpNormalDistribution[99] = 0.83891;
   mpNormalDistribution[100] = 0.84134;
   mpNormalDistribution[101] = 0.84375;
   mpNormalDistribution[102] = 0.84614;
   mpNormalDistribution[103] = 0.84849;
   mpNormalDistribution[104] = 0.85083;
   mpNormalDistribution[105] = 0.85314;
   mpNormalDistribution[106] = 0.85543;
   mpNormalDistribution[107] = 0.85769;
   mpNormalDistribution[108] = 0.85993;
   mpNormalDistribution[109] = 0.86214;
   mpNormalDistribution[110] = 0.86433;
   mpNormalDistribution[111] = 0.86650;
   mpNormalDistribution[112] = 0.86864;
   mpNormalDistribution[113] = 0.87076;
   mpNormalDistribution[114] = 0.87286;
   mpNormalDistribution[115] = 0.87493;
   mpNormalDistribution[116] = 0.87698;
   mpNormalDistribution[117] = 0.87900;
   mpNormalDistribution[118] = 0.88100;
   mpNormalDistribution[119] = 0.88298;
   mpNormalDistribution[120] = 0.88493;
   mpNormalDistribution[121] = 0.88686;
   mpNormalDistribution[122] = 0.88877;
   mpNormalDistribution[123] = 0.89065;
   mpNormalDistribution[124] = 0.89251;
   mpNormalDistribution[125] = 0.89435;
   mpNormalDistribution[126] = 0.89617;
   mpNormalDistribution[127] = 0.89796;
   mpNormalDistribution[128] = 0.89973;
   mpNormalDistribution[129] = 0.90147;
   mpNormalDistribution[130] = 0.90320;
   mpNormalDistribution[131] = 0.90490;
   mpNormalDistribution[132] = 0.90658;
   mpNormalDistribution[133] = 0.90824;
   mpNormalDistribution[134] = 0.90988;
   mpNormalDistribution[135] = 0.91149;
   mpNormalDistribution[136] = 0.91309;
   mpNormalDistribution[137] = 0.91466;
   mpNormalDistribution[138] = 0.91621;
   mpNormalDistribution[139] = 0.91774;
   mpNormalDistribution[140] = 0.91924;
   mpNormalDistribution[141] = 0.92073;
   mpNormalDistribution[142] = 0.92220;
   mpNormalDistribution[143] = 0.92364;
   mpNormalDistribution[144] = 0.92507;
   mpNormalDistribution[145] = 0.92647;
   mpNormalDistribution[146] = 0.92785;
   mpNormalDistribution[147] = 0.92922;
   mpNormalDistribution[148] = 0.93056;
   mpNormalDistribution[149] = 0.93189;
   mpNormalDistribution[150] = 0.93319;
   mpNormalDistribution[151] = 0.93448;
   mpNormalDistribution[152] = 0.93574;
   mpNormalDistribution[153] = 0.93699;
   mpNormalDistribution[154] = 0.93822;
   mpNormalDistribution[155] = 0.93943;
   mpNormalDistribution[156] = 0.94062;
   mpNormalDistribution[157] = 0.94179;
   mpNormalDistribution[158] = 0.94295;
   mpNormalDistribution[159] = 0.94408;
   mpNormalDistribution[160] = 0.94520;
   mpNormalDistribution[161] = 0.94630;
   mpNormalDistribution[162] = 0.94738;
   mpNormalDistribution[163] = 0.94845;
   mpNormalDistribution[164] = 0.94950;
   mpNormalDistribution[165] = 0.95053;
   mpNormalDistribution[166] = 0.95154;
   mpNormalDistribution[167] = 0.95254;
   mpNormalDistribution[168] = 0.95352;
   mpNormalDistribution[169] = 0.95449;
   mpNormalDistribution[170] = 0.95543;
   mpNormalDistribution[171] = 0.95637;
   mpNormalDistribution[172] = 0.95728;
   mpNormalDistribution[173] = 0.95818;
   mpNormalDistribution[174] = 0.95907;
   mpNormalDistribution[175] = 0.95994;
   mpNormalDistribution[176] = 0.96080;
   mpNormalDistribution[177] = 0.96164;
   mpNormalDistribution[178] = 0.96246;
   mpNormalDistribution[179] = 0.96327;
   mpNormalDistribution[180] = 0.96407;
   mpNormalDistribution[181] = 0.96485;
   mpNormalDistribution[182] = 0.96562;
   mpNormalDistribution[183] = 0.96638;
   mpNormalDistribution[184] = 0.96712;
   mpNormalDistribution[185] = 0.96784;
   mpNormalDistribution[186] = 0.96856;
   mpNormalDistribution[187] = 0.96926;
   mpNormalDistribution[188] = 0.96995;
   mpNormalDistribution[189] = 0.97062;
   mpNormalDistribution[190] = 0.97128;
   mpNormalDistribution[191] = 0.97193;
   mpNormalDistribution[192] = 0.97257;
   mpNormalDistribution[193] = 0.97320;
   mpNormalDistribution[194] = 0.97381;
   mpNormalDistribution[195] = 0.97441;
   mpNormalDistribution[196] = 0.97500;
   mpNormalDistribution[197] = 0.97558;
   mpNormalDistribution[198] = 0.97615;
   mpNormalDistribution[199] = 0.97670;
   mpNormalDistribution[200] = 0.97725;
   mpNormalDistribution[201] = 0.97778;
   mpNormalDistribution[202] = 0.97831;
   mpNormalDistribution[203] = 0.97882;
   mpNormalDistribution[204] = 0.97932;
   mpNormalDistribution[205] = 0.97982;
   mpNormalDistribution[206] = 0.98030;
   mpNormalDistribution[207] = 0.98077;
   mpNormalDistribution[208] = 0.98124;
   mpNormalDistribution[209] = 0.98169;
   mpNormalDistribution[210] = 0.98214;
   mpNormalDistribution[211] = 0.98257;
   mpNormalDistribution[212] = 0.98300;
   mpNormalDistribution[213] = 0.98341;
   mpNormalDistribution[214] = 0.98382;
   mpNormalDistribution[215] = 0.98422;
   mpNormalDistribution[216] = 0.98461;
   mpNormalDistribution[217] = 0.98500;
   mpNormalDistribution[218] = 0.98537;
   mpNormalDistribution[219] = 0.98574;
   mpNormalDistribution[220] = 0.98610;
   mpNormalDistribution[221] = 0.98645;
   mpNormalDistribution[222] = 0.98679;
   mpNormalDistribution[223] = 0.98713;
   mpNormalDistribution[224] = 0.98745;
   mpNormalDistribution[225] = 0.98778;
   mpNormalDistribution[226] = 0.98809;
   mpNormalDistribution[227] = 0.98840;
   mpNormalDistribution[228] = 0.98870;
   mpNormalDistribution[229] = 0.98899;
   mpNormalDistribution[230] = 0.98928;
   mpNormalDistribution[231] = 0.98956;
   mpNormalDistribution[232] = 0.98983;
   mpNormalDistribution[233] = 0.99010;
   mpNormalDistribution[234] = 0.99036;
   mpNormalDistribution[235] = 0.99061;
   mpNormalDistribution[236] = 0.99086;
   mpNormalDistribution[237] = 0.99111;
   mpNormalDistribution[238] = 0.99134;
   mpNormalDistribution[239] = 0.99158;
   mpNormalDistribution[240] = 0.99180;
   mpNormalDistribution[241] = 0.99202;
   mpNormalDistribution[242] = 0.99224;
   mpNormalDistribution[243] = 0.99245;
   mpNormalDistribution[244] = 0.99266;
   mpNormalDistribution[245] = 0.99286;
   mpNormalDistribution[246] = 0.99305;
   mpNormalDistribution[247] = 0.99324;
   mpNormalDistribution[248] = 0.99343;
   mpNormalDistribution[249] = 0.99361;
   mpNormalDistribution[250] = 0.99379;
   mpNormalDistribution[251] = 0.99396;
   mpNormalDistribution[252] = 0.99413;
   mpNormalDistribution[253] = 0.99430;
   mpNormalDistribution[254] = 0.99446;
   mpNormalDistribution[255] = 0.99461;
   mpNormalDistribution[256] = 0.99477;
   mpNormalDistribution[257] = 0.99492;
   mpNormalDistribution[258] = 0.99506;
   mpNormalDistribution[259] = 0.99520;
   mpNormalDistribution[260] = 0.99534;
   mpNormalDistribution[261] = 0.99547;
   mpNormalDistribution[262] = 0.99560;
   mpNormalDistribution[263] = 0.99573;
   mpNormalDistribution[264] = 0.99585;
   mpNormalDistribution[265] = 0.99598;
   mpNormalDistribution[266] = 0.99609;
   mpNormalDistribution[267] = 0.99621;
   mpNormalDistribution[268] = 0.99632;
   mpNormalDistribution[269] = 0.99643;
   mpNormalDistribution[270] = 0.99653;
   mpNormalDistribution[271] = 0.99664;
   mpNormalDistribution[272] = 0.99674;
   mpNormalDistribution[273] = 0.99683;
   mpNormalDistribution[274] = 0.99693;
   mpNormalDistribution[275] = 0.99702;
   mpNormalDistribution[276] = 0.99711;
   mpNormalDistribution[277] = 0.99720;
   mpNormalDistribution[278] = 0.99728;
   mpNormalDistribution[279] = 0.99736;
   mpNormalDistribution[280] = 0.99744;
   mpNormalDistribution[281] = 0.99752;
   mpNormalDistribution[282] = 0.99760;
   mpNormalDistribution[283] = 0.99767;
   mpNormalDistribution[284] = 0.99774;
   mpNormalDistribution[285] = 0.99781;
   mpNormalDistribution[286] = 0.99788;
   mpNormalDistribution[287] = 0.99795;
   mpNormalDistribution[288] = 0.99801;
   mpNormalDistribution[289] = 0.99807;
   mpNormalDistribution[290] = 0.99813;
   mpNormalDistribution[291] = 0.99819;
   mpNormalDistribution[292] = 0.99825;
   mpNormalDistribution[293] = 0.99831;
   mpNormalDistribution[294] = 0.99836;
   mpNormalDistribution[295] = 0.99841;
   mpNormalDistribution[296] = 0.99846;
   mpNormalDistribution[297] = 0.99851;
   mpNormalDistribution[298] = 0.99856;
   mpNormalDistribution[299] = 0.99861;
   mpNormalDistribution[300] = 0.99865;
   mpNormalDistribution[301] = 0.99869;
   mpNormalDistribution[302] = 0.99874;
   mpNormalDistribution[303] = 0.99878;
   mpNormalDistribution[304] = 0.99882;
   mpNormalDistribution[305] = 0.99886;
   mpNormalDistribution[306] = 0.99889;
   mpNormalDistribution[307] = 0.99893;
   mpNormalDistribution[308] = 0.99896;
   mpNormalDistribution[309] = 0.99900;
   mpNormalDistribution[310] = 0.99903;
   mpNormalDistribution[311] = 0.99906;
   mpNormalDistribution[312] = 0.99910;
   mpNormalDistribution[313] = 0.99913;
   mpNormalDistribution[314] = 0.99916;
   mpNormalDistribution[315] = 0.99918;
   mpNormalDistribution[316] = 0.99921;
   mpNormalDistribution[317] = 0.99924;
   mpNormalDistribution[318] = 0.99926;
   mpNormalDistribution[319] = 0.99929;
   mpNormalDistribution[320] = 0.99931;
   mpNormalDistribution[321] = 0.99934;
   mpNormalDistribution[322] = 0.99936;
   mpNormalDistribution[323] = 0.99938;
   mpNormalDistribution[324] = 0.99940;
   mpNormalDistribution[325] = 0.99942;
   mpNormalDistribution[326] = 0.99944;
   mpNormalDistribution[327] = 0.99946;
   mpNormalDistribution[328] = 0.99948;
   mpNormalDistribution[329] = 0.99950;
   mpNormalDistribution[330] = 0.99952;
   mpNormalDistribution[331] = 0.99953;
   mpNormalDistribution[332] = 0.99955;
   mpNormalDistribution[333] = 0.99957;
   mpNormalDistribution[334] = 0.99958;
   mpNormalDistribution[335] = 0.99960;
   mpNormalDistribution[336] = 0.99961;
   mpNormalDistribution[337] = 0.99962;
   mpNormalDistribution[338] = 0.99964;
   mpNormalDistribution[339] = 0.99965;
   mpNormalDistribution[340] = 0.99966;
   mpNormalDistribution[341] = 0.99968;
   mpNormalDistribution[342] = 0.99969;
   mpNormalDistribution[343] = 0.99970;
   mpNormalDistribution[344] = 0.99971;
   mpNormalDistribution[345] = 0.99972;
   mpNormalDistribution[346] = 0.99973;
   mpNormalDistribution[347] = 0.99974;
   mpNormalDistribution[348] = 0.99975;
   mpNormalDistribution[349] = 0.99976;
   mpNormalDistribution[350] = 0.99977;
   mpNormalDistribution[351] = 0.99978;
   mpNormalDistribution[352] = 0.99978;
   mpNormalDistribution[353] = 0.99979;
   mpNormalDistribution[354] = 0.99980;
   mpNormalDistribution[355] = 0.99981;
   mpNormalDistribution[356] = 0.99981;
   mpNormalDistribution[357] = 0.99982;
   mpNormalDistribution[358] = 0.99983;
   mpNormalDistribution[359] = 0.99983;
   mpNormalDistribution[360] = 0.99984;
   mpNormalDistribution[361] = 0.99985;
   mpNormalDistribution[362] = 0.99985;
   mpNormalDistribution[363] = 0.99986;
   mpNormalDistribution[364] = 0.99986;
   mpNormalDistribution[365] = 0.99987;
   mpNormalDistribution[366] = 0.99987;
   mpNormalDistribution[367] = 0.99988;
   mpNormalDistribution[368] = 0.99988;
   mpNormalDistribution[369] = 0.99989;
   mpNormalDistribution[370] = 0.99989;
   mpNormalDistribution[371] = 0.99990;
   mpNormalDistribution[372] = 0.99990;
   mpNormalDistribution[373] = 0.99990;
   mpNormalDistribution[374] = 0.99991;
   mpNormalDistribution[375] = 0.99991;
   mpNormalDistribution[376] = 0.99992;
   mpNormalDistribution[377] = 0.99992;
   mpNormalDistribution[378] = 0.99992;
   mpNormalDistribution[379] = 0.99992;
   mpNormalDistribution[380] = 0.99993;
   mpNormalDistribution[381] = 0.99993;
   mpNormalDistribution[382] = 0.99993;
   mpNormalDistribution[383] = 0.99994;
   mpNormalDistribution[384] = 0.99994;
   mpNormalDistribution[385] = 0.99994;
   mpNormalDistribution[386] = 0.99994;
   mpNormalDistribution[387] = 0.99995;
   mpNormalDistribution[388] = 0.99995;
   mpNormalDistribution[389] = 0.99995;
   mpNormalDistribution[390] = 0.99995;
   mpNormalDistribution[391] = 0.99995;
   mpNormalDistribution[392] = 0.99996;
   mpNormalDistribution[393] = 0.99996;
   mpNormalDistribution[394] = 0.99996;
   mpNormalDistribution[395] = 0.99996;
   mpNormalDistribution[396] = 0.99996;
   mpNormalDistribution[397] = 0.99996;
   mpNormalDistribution[398] = 0.99997;
   mpNormalDistribution[399] = 0.99997;
}