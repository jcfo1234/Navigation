//----------------------------------------------------------------------
// LeastSquares.cpp
// Source file of LeastSquares.h
// Source code for least squares computations
//----------------------------------------------------------------------
#include "..\inc\LeastSquares.h"

//----------------------------------------------------------------------------
void LeastSquaresSolution(BOOLEANO bTheECEFLS_, BOOLEANO bTheGeoidLS_, BOOLEANO bTheGeoidConstrainedLS_,
                          DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_, DOUBLE* pdTheHeight_, DOUBLE* pdTheTShift_,
                          DOUBLE* pdTheTruthLatitude_, DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_,
                          string* psTheTruthFullFilePath_, Distributions* pclTheProbabilityDistributions_, string* psTheObsCovKey_,
                          string* psTheFullFilePath_, DOUBLE* pdTheMeasStd_, DOUBLE* pdTheConfidenceLevel_, DOUBLE* pdTheTestPower_,
                          DOUBLE* pdTheHeightConstrained_, DOUBLE* pdTheHeightObsStd_, BOOLEANO bTheBlunderEnable_)
{
   // Truth Trajectory file
   ifstream clTruthFile;
   // Container with satellite data
   vector<EpochInfo> vEpochData;

   PopulateSVEpochData(*psTheFullFilePath_, &vEpochData);
   SendSVEpochDatatoText(*psTheFullFilePath_, *psTheFullFilePath_);

   if (psTheTruthFullFilePath_ != NULL)
      clTruthFile.open(*psTheTruthFullFilePath_);

   // Do ECEF Least Squares Analysis
   if (bTheECEFLS_)
   {
      cout << "\n\nECEF Least squares processing..." << endl;
      cout << "\n\nGPS Time, Latitude, Longitude, Height" << endl;
      // Point of expansion position
      CMatrix clECEFPos("ECEF Point of Expansion", 4, 1);
      CMatrix clErrorStateCovariance("P", 4, 4);
      // Objects and containers used for statistic computations
      CMatrix clECEFCoordinates("ECEF Coordinates", 3, 1);
      CMatrix clGeoidCoordinates("Geoid Coordinates", 3, 1);
      CMatrix clRle("ECEF to LLF", 4, 4);
      TrajectoryInfo stTrajInfo;
      vector<ErrorInfo> vPosErrorTimeSeries;
      vector<ReliablityTest> vReliabilityTimeSeries;
      // Initialize point of expansion
      SetPointOfExpansionEcef(*pdTheLatitude_, *pdTheLongitude_, *pdTheHeight_, *pdTheTShift_, &clECEFPos);
      // Compute solution, epoch by epoch
      for (vector<EpochInfo>::iterator itEpoch = vEpochData.begin(); itEpoch != vEpochData.end(); itEpoch++)
      {
         clErrorStateCovariance = ComputeECEFPosition(&(*itEpoch), pclTheProbabilityDistributions_, &clECEFPos, *pdTheMeasStd_, *psTheObsCovKey_, *pdTheConfidenceLevel_, *pdTheTestPower_, &vReliabilityTimeSeries, bTheBlunderEnable_);
         FillCoordinates(&clECEFPos, &clECEFCoordinates);
         ECEFToGeoid(&clECEFCoordinates, &clGeoidCoordinates);
         // Truth trajectory file provided
         if (psTheTruthFullFilePath_ != NULL)
         {
            AlignEpochs(clTruthFile, itEpoch, &stTrajInfo, &vEpochData);
         }
         // Truth trajectory file not provided
         else
         {
            FillTrajectory(&stTrajInfo, *pdTheTruthLatitude_, *pdTheTruthLongitude_, *pdTheTruthHeight_, itEpoch->dGPSTime);
         }
         clRle = ECEFToLLF(&clGeoidCoordinates);
         FillPositionErrorVector(&vPosErrorTimeSeries, &clGeoidCoordinates, &stTrajInfo, &clErrorStateCovariance, &clRle, itEpoch->cNumSV, *pdTheMeasStd_);
         cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidCoordinates.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidCoordinates.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(9) << clGeoidCoordinates.GetComponent(2, 0) << endl;
      }
      // Print coordinates
      PrintPlots(*psTheFullFilePath_, &vPosErrorTimeSeries, &vReliabilityTimeSeries, "ECEF");
   }

   // Do LLF Least Squares Analysis
   if (bTheGeoidLS_)
   {
      cout << "\n\nCurvilinear Least squares processing..." << endl;
      cout << "\n\nGPS Time, Latitude, Longitude, Height" << endl;
      // Truth Data variables
      string sFileRead;
      // Initial position
      CMatrix clGeoidPos("Geoid Initial Point of Expansion", 4, 1);
      CMatrix clErrorStateCovariance("P", 4, 4);
      // Objects and containers used for statistic computations
      CMatrix clIdentity("Identity", 4, 4);
      TrajectoryInfo stTrajInfo;
      vector<ErrorInfo> vPosErrorTimeSeries;
      vector<ReliablityTest> vReliabilityTimeSeries;
      // Initialize point of expansion
      clIdentity.SetIdentity();
      SetPointOfExpansionGeoid(*pdTheLatitude_, *pdTheLongitude_, *pdTheHeight_, *pdTheTShift_, &clGeoidPos);
      // Compute solution, epoch by epoch
      for (vector<EpochInfo>::iterator itEpoch = vEpochData.begin(); itEpoch != vEpochData.end(); itEpoch++)
      {
         clErrorStateCovariance = ComputeGeoidPosition(&(*itEpoch), pclTheProbabilityDistributions_, &clGeoidPos, *pdTheMeasStd_, *psTheObsCovKey_, *pdTheConfidenceLevel_, *pdTheTestPower_, &vReliabilityTimeSeries, bTheBlunderEnable_);
         // Truth trajectory file provided
         if (psTheTruthFullFilePath_ != NULL)
         {
            AlignEpochs(clTruthFile, itEpoch, &stTrajInfo, &vEpochData);
         }
         else
         {
            FillTrajectory(&stTrajInfo, *pdTheTruthLatitude_, *pdTheTruthLongitude_, *pdTheTruthHeight_, itEpoch->dGPSTime);
         }
         FillPositionErrorVector(&vPosErrorTimeSeries, &clGeoidPos, &stTrajInfo, &clErrorStateCovariance, &clIdentity, itEpoch->cNumSV, *pdTheMeasStd_);
         cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidPos.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidPos.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(9) << clGeoidPos.GetComponent(2, 0) << endl;
      }
      // Print coordinates
      PrintPlots(*psTheFullFilePath_, &vPosErrorTimeSeries, &vReliabilityTimeSeries, "CURVILINEAR");
   }

   // Do LLF Least Squares Analysis with Height constrained
   if (bTheGeoidConstrainedLS_)
   {
      cout << "\n\nCurvilinear-Height Constrained Least squares processing..." << endl;
      cout << "\n\nGPS Time, Latitude, Longitude, Height" << endl;
      // Initial position
      CMatrix clGeoidPos("Geoid Initial Point of Expansion", 4, 1);
      CMatrix clErrorStateCovariance("P", 4, 4);
      // Objects and containers used for statistic computations
      CMatrix clIdentity("Identity", 4, 4);
      TrajectoryInfo stTrajInfo;
      vector<ErrorInfo> vPosErrorTimeSeries;
      vector<ReliablityTest> vReliabilityTimeSeries;
      // Initialize point of expansion
      clIdentity.SetIdentity();
      SetPointOfExpansionGeoid(*pdTheLatitude_, *pdTheLongitude_, *pdTheHeight_, *pdTheTShift_, &clGeoidPos);
      // Compute solution, epoch by epoch
      for (vector<EpochInfo>::iterator itEpoch = vEpochData.begin(); itEpoch != vEpochData.end(); itEpoch++)
      {
         clErrorStateCovariance = ComputeGeoidPositionConstrained(&(*itEpoch), pclTheProbabilityDistributions_, &clGeoidPos, *pdTheMeasStd_, *psTheObsCovKey_, *pdTheHeightConstrained_, *pdTheConfidenceLevel_, *pdTheTestPower_, *pdTheHeightObsStd_, &vReliabilityTimeSeries, bTheBlunderEnable_);
         // Truth trajectory file provided
         if (psTheTruthFullFilePath_ != NULL)
         {
            AlignEpochs(clTruthFile, itEpoch, &stTrajInfo, &vEpochData);
         }
         else
         {
            FillTrajectory(&stTrajInfo, *pdTheTruthLatitude_, *pdTheTruthLongitude_, *pdTheTruthHeight_, itEpoch->dGPSTime);
         }
         FillPositionErrorVector(&vPosErrorTimeSeries, &clGeoidPos, &stTrajInfo, &clErrorStateCovariance, &clIdentity, itEpoch->cNumSV, *pdTheMeasStd_);
         cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidPos.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidPos.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(9) << clGeoidPos.GetComponent(2, 0) << endl;
      }
      // Print coordinates
      PrintPlots(*psTheFullFilePath_, &vPosErrorTimeSeries, &vReliabilityTimeSeries, "CURVILINEARHGTCONSTRAINT");
   }

   // Close Truth trajectory file
   if (psTheTruthFullFilePath_ != NULL)
      clTruthFile.close();
}
