// Estimation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

// Truth position values
DOUBLE dTheTruthLatitude = 51.0799537305556;
DOUBLE dTheTruthLongitude = -114.133848205556;
DOUBLE dTheTruthHeight = 1118.792;
// Default point of expansion
DOUBLE dTheLatitude = 51.1225;
DOUBLE dTheLongitude = -114.0133;
DOUBLE dTheHeight = 1084;
DOUBLE dTheTShift = 0;
DOUBLE dTheVelEast = 0;
DOUBLE dTheVelNorth = 0;
DOUBLE dTheVelVert = 0;
DOUBLE dTheTDrift = 0;
DOUBLE dTheHeightConstrained = 1118.792;
DOUBLE dTheInitPosStd = 1000;
DOUBLE dTheInitVelStd = 100;
DOUBLE dTheInitRadius = 1;
// Default position error standard deviation
DOUBLE dTheMeasStd = 10;
DOUBLE dTheDopplerMeasStd = 1;
DOUBLE dTheHeightObsStd = 0.1;
// String with full path to binary data file
string stTheFullFilePath("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO620\\Labs\\Lab1\\Data.bin");
string sTheTruthReference("C:\\Users\\Juanchelsea\\Desktop\\MASTERS\\ENGO620\\Labs\\Lab1\\Data.txt");
// Least squares algorithm
BOOLEANO bTheECEFLS = FALSE;
BOOLEANO bTheGeoidLS = FALSE;
BOOLEANO bTheGeoidConstrainedLS = FALSE;
BOOLEANO bTheVRW = FALSE;
BOOLEANO bThePRW = FALSE;
BOOLEANO bTheUseDoppler = FALSE;
BOOLEANO bTheBlunderEnable = FALSE;
BOOLEANO bTheKalmanPolar = FALSE;
string stTheObsCovKey("CONSTANT");
// Statistical Tests
DOUBLE dTheConfidenceLevel = 0.05;
DOUBLE dTheTestPower = 0.1;
Distributions clTheProbabilityDistributions;
// Power spectral noise densities
DOUBLE dThePRWQ = 1;
DOUBLE dTheVRWQ = 1;
DOUBLE dTheh0 = 1e-21;
DOUBLE dThehm2 = 1e-20;
DOUBLE dTheClkBiasQ = dTheh0 / 2 * pow(SPEED_OF_LIGHT, 2);
DOUBLE dTheClkDriftQ = dThehm2 * pow(SPEED_OF_LIGHT, 2) * 2 * COORDTRANS_PI;
// Control switch options
const CHAR* acOptions[ESTIMATION_NUM_OPTIONS] = {"THELAT", "THELNG", "THEHGT", "THETSH", "FILE", "ECEF", "GEOID", "GEOCNST", "CONSTANT", "ELEVATION", "ALPHA", "BETA", "HCNST", "TRUTH",
                                                 "VEAST", "VNORTH", "VVERT", "TDRIFT", "ELEVVEL", "PRWQ", "VRWQ", "TSHQ", "TDRIFTQ", "PINITP", "VINITP", "VRW", "USEDOPPLER", "PRW", 
                                                 "BLUNDER", "POLAR", "RADIUS"};
vector<string> vOptions(acOptions, acOptions + ESTIMATION_NUM_OPTIONS);
map<string, void*> mpTheSwitchOptions;
// Error codes
map<INT, string> mpTheErrorCodes;

INT main(INT iargc_, CHAR** ppcargv_)
{
   // Truth trajectory string pointer
   string* psTruth = NULL;

   // Initialize control switch options
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGELATITUDE_INDEX)] = &dTheLatitude;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGELONGITUDE_INDEX)] = &dTheLongitude;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGEHEIGHT_INDEX)] = &dTheHeight;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGETIMESHIFT_INDEX)] = &dTheTShift;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGEFILE_INDEX)] = &stTheFullFilePath;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_ECEFENABLE_INDEX)] = &bTheECEFLS;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_GEOIDENABLE_INDEX)] = &bTheGeoidLS;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_GEOIDCONSTRAINEDENABLE_INDEX)] = &dTheHeightObsStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_COVCONST_INDEX)] = &dTheMeasStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_COVELEV_INDEX)] = &dTheMeasStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CONFIDENCELEVEL_INDEX)] = &dTheConfidenceLevel;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_TESTPOWER_INDEX)] = &dTheTestPower;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_HEIGHTCNST_INDEX)] = &dTheHeightConstrained;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_TRUTHFILE_INDEX)] = &sTheTruthReference;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGEVEAST_INDEX)] = &dTheVelEast;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGEVNORTH_INDEX)] = &dTheVelNorth;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGEVVERT_INDEX)] = &dTheVelVert;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_CHANGECLOCKDRIFT_INDEX)] = &dTheTDrift;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_COVVELEV_INDEX)] = &dTheDopplerMeasStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_SPECTRALPRWQ_INDEX)] = &dThePRWQ;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_SPECTRALVRWQ_INDEX)] = &dTheVRWQ;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_SPECTRALTSHQ_INDEX)] = &dTheh0;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_SPECTRALTDRIFTQ_INDEX)] = &dThehm2;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITPOSSTD_INDEX)] = &dTheInitPosStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_INITVELSTD_INDEX)] = &dTheInitVelStd;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_VRWENABLE_INDEX)] = &bTheVRW;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_DOPPLERENABLE_INDEX)] = &bTheUseDoppler;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_PRWENABLE_INDEX)] = &bThePRW;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_BLUNDERENABLE_INDEX)] = &bTheBlunderEnable;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POLARENABLE_INDEX)] = &bTheKalmanPolar;
   mpTheSwitchOptions[vOptions.at(ESTIMATION_POLARRADIUS_INDEX)] = &dTheInitRadius;

   // Error map table
   mpTheErrorCodes[1]  = "Invalid latitude vector dimensions (CalcN)";
   mpTheErrorCodes[2]  = "Invalid latitude vector dimensions (CalcM)";
   mpTheErrorCodes[3]  = "Invalid Geoid and ECEF vector dimensions (GeoidToECEF)";
   mpTheErrorCodes[4]  = "Invalid Geoid and ECEF vector dimensions (ECEFToGeoid)";
   mpTheErrorCodes[5]  = "Invalid dimensions of point of expansion (CalculatePseudorangeResiduals)";
   mpTheErrorCodes[6]  = "Invalid dimensions of point of expansion (ComputeDesignECEFMatrixH)";
   mpTheErrorCodes[7]  = "Invalid dimensions of point of expansion (ComputeDesignLLFMatrixH)";
   mpTheErrorCodes[8]  = "Invalid dimensions of point of expansion (ComputePseudorangeGeoidT)";
   mpTheErrorCodes[9]  = "Invalid dimensions of point of expansion (ComputeGeoidToLlfS)";
   mpTheErrorCodes[10] = "Invalid observation covariance matrix key type (ComputeECEFPosition)";
   mpTheErrorCodes[11] = "Invalid observation covariance matrix key type (ComputeGeoidPosition)";
   mpTheErrorCodes[12] = "Invalid dimensions of point of expansion (SetPointOfExpansionGeoid)";
   mpTheErrorCodes[13] = "Invalid dimensions of point of expansion (SetPointOfExpansionEcef)";
   mpTheErrorCodes[14] = "Invalid option selected from LAT, LNG, HGT, TSH, FILE, ECEF, GEOID, GEOCNST, CONSTANT, ELEVATION";
   mpTheErrorCodes[15] = "Invalid switch delimiter from / or -";
   mpTheErrorCodes[16] = "Invalid ECEF vector dimensions (FillECEFCoordinates)";
   mpTheErrorCodes[17] = "Invalid vector dimensions (KGeoidToLLFS)";
   mpTheErrorCodes[18] = "Invalid vector dimensions (KPseudorangeGeoidT)";
   mpTheErrorCodes[19] = "Invalid vector dimensions (KPseudorangeMatrixH)";
   mpTheErrorCodes[20] = "Invalid vector dimensions (KPseudorateMatrixH)";
   mpTheErrorCodes[21] = "Invalid vector dimensions (KPseudorateGeoidTS)";
   mpTheErrorCodes[22] = "Invalid vector dimensions (KPseudorangeResiduals)";
   mpTheErrorCodes[23] = "Invalid vector dimensions (KPseudorateResiduals)";
   mpTheErrorCodes[24] = "Invalid vector dimensions (KPolarJacobian)";
   mpTheErrorCodes[25] = "Invalid vector dimensions (KPolarTransformH)";

   try
   {
      // Parse arguments
      for (INT i = 1; i < iargc_; i++)
      {
         // Transform one string to uppercase and keep original string
         string stOption(ppcargv_[i]);
         string stOptionOrig(ppcargv_[i]);
         transform(stOption.begin(), stOption.end(), stOption.begin(), toupper);
         //-----------------------------------------------------------------------------------
         // Control switches are:
         // /-lat<Latitude in Degrees>             : Change latitude default value
         // /-lng<Longitude in Degrees>            : Change longitude default value
         // /-hgt<Height in Degrees>               : Change height default value
         // /-tsh<Time shift in meters>            : Change time shift default value
         // /-file<Full path to binary file>       : Change file path default value
         // /-ECEF                                 : Enable ECEF Least squares computations
         // /-GEOID                                : Enable GEOID Least squares computations
         // /-GEOCNST<Height Std>                  : Enable GEOID Least squares, height constrained
         // /-CONSTANT<Standard Deviation>         : Constant diagonal pseudorange Observation Covariance
         // /-ELEVATION<Standard Deviation>        : Elevation dependant pseudorange Observation Covariance
         // /-ALPHA<Significance Level>            : Change test confidence level
         // /-BETA<Power of the Test>              : Change test power or reach
         // /-HCNST<Height constrained value>      : Change the default Height constrain
         // /-TRUTH<Full path to truth file>       : Change the default path to truth trajectory file
         // /-VEAST<Velocity east m/s>             : Change of east velocity default value
         // /-VNORTH<Velocity north m/s>           : Change of north velocity default value
         // /-VVERT<Velocity vertical m/s>         : Change of vertical velocity default value
         // /-TDRIFT<Clock drift m/s>              : Change of clock drift default value
         // /-ELEVVEL<Standard deviation>          : Elevation dependant doppler Observation Covariance
         // /-PRWQ<Noise density m/s-sqrt(Hz)>     : Position random walk spectral density
         // /-VRWQ<Noise density m/s^2-sqrt(Hz>    : Velocity random walk spectral density
         // /-TSHQ<h0 oscillator parameter>        : h0 oscillator parameter
         // /-TDRIFTQ<h-2 oscillator parameter>    : h-2 oscillator parameter
         // /-PINITP<Initial position Std Dev>     : Initial position standard deviation
         // /-VINITP<Initial velocity Std Dev>     : Initial velocity standard deviation
         // /-VRW                                  : Enable Kalman velocity random walk
         // /-USEDOPPLER                           : Enalbe using doppler measurements in Kalman filter
         // /-PRW                                  : Enable Kalman position random walk
         // /-BLUNDER                              : Enable Blunder detection algorithm
         // /-POLAR                                : Enable Kalman filter with polar coordinates
         //-----------------------------------------------------------------------------------
         switch (stOption.at(0))
         {
            // Allowed switch delimiters
         case '/':
         case '-':
            for (vector<string>::iterator itOptions = vOptions.begin(); itOptions != vOptions.end(); itOptions++)
            {
               // Option found
               if (stOption.find(*itOptions) != string::npos)
               {
                  // Full file path assignment or Truth file path assignment
                  if (*itOptions == vOptions.at(ESTIMATION_CHANGEFILE_INDEX) || *itOptions == vOptions.at(ESTIMATION_TRUTHFILE_INDEX))
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     *static_cast<string*>(mpTheSwitchOptions.at(*itOptions)) = stChangeValue;
                     // Truth trajectory file pointer
                     if (*itOptions == vOptions.at(ESTIMATION_TRUTHFILE_INDEX))
                        psTruth = static_cast<string*>(mpTheSwitchOptions.at(*itOptions));
                  }
                  // ECEF or Geoid least squares
                  else if ( (*itOptions == vOptions.at(ESTIMATION_ECEFENABLE_INDEX)) || 
                            (*itOptions == vOptions.at(ESTIMATION_GEOIDENABLE_INDEX)) ||
                            (*itOptions == vOptions.at(ESTIMATION_VRWENABLE_INDEX)) ||
                            (*itOptions == vOptions.at(ESTIMATION_DOPPLERENABLE_INDEX)) ||
                            (*itOptions == vOptions.at(ESTIMATION_PRWENABLE_INDEX)) ||
                            (*itOptions == vOptions.at(ESTIMATION_BLUNDERENABLE_INDEX)) ||
                            (*itOptions == vOptions.at(ESTIMATION_POLARENABLE_INDEX)) )
                  {
                     *((BOOLEANO*)mpTheSwitchOptions.at(*itOptions)) = TRUE;
                  }
                  // Latitude, Longitude, Height, Time shift or Measurement standard deviation assignment
                  else
                  {
                     string stChangeValue = stOptionOrig.substr(stOption.find(*itOptions) + itOptions->length());
                     // Observation covariance matrix method selection
                     if (*itOptions == vOptions.at(ESTIMATION_COVCONST_INDEX) || *itOptions == vOptions.at(ESTIMATION_COVELEV_INDEX))
                     {
                        stTheObsCovKey = *itOptions;
                     }
                     // Implement Height constrained estimation
                     else if (*itOptions == vOptions.at(ESTIMATION_GEOIDCONSTRAINEDENABLE_INDEX) || *itOptions == vOptions.at(ESTIMATION_HEIGHTCNST_INDEX))
                     {
                        bTheGeoidConstrainedLS = TRUE;
                     }
                     *((DOUBLE*)mpTheSwitchOptions.at(*itOptions)) = stod(stChangeValue);
                  }
                  break;
               }
               // Option not found
               if (itOptions == vOptions.end())
               {
                  throw 14;
               }
            }
            break;
            // Invalid switch delimiter
         default:
            throw 15;
         }
      }

      if (bTheECEFLS || bTheGeoidLS || bTheGeoidConstrainedLS)
      {
         // Compute least squares solution
         LeastSquaresSolution(bTheECEFLS, bTheGeoidLS, bTheGeoidConstrainedLS,
                              &dTheLatitude, &dTheLongitude, &dTheHeight, &dTheTShift,
                              &dTheTruthLatitude, &dTheTruthLongitude, &dTheTruthHeight,
                              psTruth, &clTheProbabilityDistributions, &stTheObsCovKey,
                              &stTheFullFilePath, &dTheMeasStd, &dTheConfidenceLevel, &dTheTestPower,
                              &dTheHeightConstrained, &dTheHeightObsStd, bTheBlunderEnable);
      }

      if (bTheVRW || bThePRW)
      {
         // Compute Kalman filter solution
         KalmanSolution(bTheVRW, bThePRW, bTheUseDoppler, &dTheLatitude, &dTheLongitude, &dTheHeight,
                        &dTheVelEast, &dTheVelNorth, &dTheVelVert, &dTheTShift, &dTheTDrift, &dTheTruthLatitude,
                        &dTheTruthLongitude, &dTheTruthHeight, psTruth, &stTheFullFilePath, &dTheMeasStd, &dTheDopplerMeasStd, 
                        &dThePRWQ, &dTheVRWQ, &dTheClkBiasQ, &dTheClkDriftQ, &dTheInitPosStd, &dTheInitVelStd);
      }

      if (bTheKalmanPolar)
      {
         // Compute Kalman filter solution using polar coordinates
         KPolarSolution(&dTheLatitude, &dTheLongitude, &dTheHeight, &dTheVelEast, &dTheVelNorth, &dTheVelVert, 
                        &dTheTShift, &dTheTDrift, &dTheTruthLatitude, &dTheTruthLongitude, &dTheTruthHeight, 
                        psTruth, &stTheFullFilePath, &dTheMeasStd, &dTheDopplerMeasStd, &dThePRWQ, &dTheVRWQ, 
                        &dTheClkBiasQ, &dTheClkDriftQ, &dTheInitPosStd, &dTheInitVelStd, bTheUseDoppler, &dTheInitRadius);
      }
   }
   // Catch errors discovered in the program
   catch (INT iError)
   {
      cout << mpTheErrorCodes.at(iError) << endl;
      return EXIT_FAILURE;
   }
   // Catch unknown errors discovered in sub-functions and dependencies
   catch (...)
   {
      cout << "Unknown Error" << endl;
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}

