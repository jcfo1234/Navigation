//----------------------------------------------------------------------
// Kalman.cpp
// Source file for Kalman.h
// Computes position based on Kalman Filter algorithm
//----------------------------------------------------------------------
#include "..\inc\Kalman.h"

//----------------------------------------------------------------------------
void KalmanSolution(BOOLEANO bTheVRW_, BOOLEANO bThePRW_, BOOLEANO bUsePseudorate_,
                    DOUBLE* pdTheInitLatitude, DOUBLE* pdTheInitLongitude, DOUBLE* pdTheInitHeight,
                    DOUBLE* pdTheInitVelEast, DOUBLE* pdTheInitVelNorth, DOUBLE* pdTheInitVelUp,
                    DOUBLE* pdTheInitClkBias, DOUBLE* pdTheInitClkDrift, DOUBLE* pdTheTruthLatitude_,
                    DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_, string* psTheTruthFilePath_,
                    string* psTheFullFilePath_, DOUBLE* pdThePSRStd_, DOUBLE* pdThePSRRateStd_,
                    DOUBLE* pdThePRWSpectralQ_, DOUBLE* pdTheVRWSpectralQ_, DOUBLE* pdTheClkBiasQ_,
                    DOUBLE* pdTheClkDriftQ_, DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_)
{
   // Truth Trajectory file
   ifstream clTruthFile;
   // Container with satellite data
   vector<EpochInfo> vEpochData;

   PopulateSVEpochData(*psTheFullFilePath_, &vEpochData);
   SendSVEpochDatatoText(*psTheFullFilePath_, *psTheFullFilePath_);

   if (psTheTruthFilePath_ != NULL)
      clTruthFile.open(*psTheTruthFilePath_);

   // Position Random Walk model
   if (bThePRW_)
   {
      cout << "\n\nPosition Random Walk Kalman processing..." << endl;
      cout << "\n\nGPS Time, Latitude, Longitude, Height" << endl;
      // State space matrix
      CMatrix clStateSpaceF(PRWStateSpaceMatrixF());
      // Noise shaping matrix
      CMatrix clNoiseShapeG(PRWNoiseShapeMatrixG());
      // Spectral density matrix
      CMatrix clSpectralDensityQ(KSpectralDensityQ(pdThePRWSpectralQ_, pdTheClkBiasQ_, pdTheClkDriftQ_));
      // System state vector
      CMatrix clGeoidSystemState(PRWInitializeState(pdTheInitLatitude, pdTheInitLongitude, pdTheInitHeight, pdTheInitClkBias, pdTheInitClkDrift));
      // System state covariance matrix
      CMatrix clSystemCovarianceP(PRWInitializeStateCovarianceP(pdTheInitPosStd_));
      // Trajectory epoch
      TrajectoryInfo stTrajInfo;
      TrajectoryInfo stKalmanTraj;
      // Save trajectories in containers
      vector<TrajectoryInfo> vTruthTrajectory;
      vector<TrajectoryInfo> vKalmanTrajectory;
      vector<PRWErrorInfo> vPRWError;
      // GPS Time from epoch before
      DOUBLE dPrevGPSTime = 0;
      // Perform Kalman filter operations
      for (vector<EpochInfo>::iterator itEpoch = vEpochData.begin(); itEpoch != vEpochData.end(); itEpoch++)
      {
         // Truth trajectory file provided
         if (psTheTruthFilePath_ != NULL)
         {
            AlignEpochs(clTruthFile, itEpoch, &stTrajInfo, &vEpochData);
         }
         else
         {
            FillTrajectory(&stTrajInfo, *pdTheTruthLatitude_, *pdTheTruthLongitude_, *pdTheTruthHeight_, itEpoch->dGPSTime);
         }
         // Kalman prediction
         if (itEpoch != vEpochData.begin())
         {
            clStateSpaceF = PRWStateSpaceMatrixF();
            clNoiseShapeG = PRWNoiseShapeMatrixG();
            KPredictionLoop(&clStateSpaceF, &clNoiseShapeG, &clSpectralDensityQ, &clGeoidSystemState, &clSystemCovarianceP, itEpoch->dGPSTime - dPrevGPSTime);
         }
         // Kalman update
         PRWUpdateLoop(&clGeoidSystemState, &clSystemCovarianceP, &(*itEpoch), *pdThePSRStd_);
         // Save previous GPS Time for step size integration
         dPrevGPSTime = itEpoch->dGPSTime;
         // Save estimated trajectory into a container
         FillTrajectory(&stKalmanTraj, RadToDeg(clGeoidSystemState.GetComponent(0, 0)), RadToDeg(clGeoidSystemState.GetComponent(1, 0)), clGeoidSystemState.GetComponent(2, 0), itEpoch->dGPSTime);
         // Fill Kalman state error structure
         PRWFillPositionError(&stTrajInfo, &stKalmanTraj, &clSystemCovarianceP, &vPRWError, itEpoch->cNumSV);
         vTruthTrajectory.push_back(stTrajInfo);
         vKalmanTrajectory.push_back(stKalmanTraj);
         // Output data on console window
         cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidSystemState.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidSystemState.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(9) << clGeoidSystemState.GetComponent(2, 0) << endl;
      }
      KalmanPosErrorsPlots(&vTruthTrajectory, &vKalmanTrajectory, *psTheFullFilePath_, "PRW");
      KalmanCrossCorrPlots(&vPRWError, *psTheFullFilePath_, "PRW");
      KalmanTrajectory(&vTruthTrajectory, *psTheFullFilePath_, "PRW", "TruthTrajectory");
      KalmanTrajectory(&vKalmanTrajectory, *psTheFullFilePath_, "PRW", "KalmanTrajectory");
   }

   // Velocity random walk model
   if (bTheVRW_)
   {
      cout << "\n\nVelocity Random Walk Kalman processing..." << endl;
      cout << "\n\nGPS Time, Latitude, Longitude, Height, VelocityEast, VelocityNorth, VelocityUp" << endl;
      // State space matrix
      CMatrix clStateSpaceF(VRWStateSpaceMatrixF());
      // Noise shaping matrix
      CMatrix clNoiseShapeG(VRWNoiseShapeMatrixG());
      // Spectral density matrix
      CMatrix clSpectralDensityQ(KSpectralDensityQ(pdTheVRWSpectralQ_, pdTheClkBiasQ_, pdTheClkDriftQ_));
      // System state vector
      CMatrix clGeoidSystemState(VRWInitializeState(pdTheInitLatitude, pdTheInitLongitude, pdTheInitHeight, pdTheInitVelEast, pdTheInitVelNorth, pdTheInitVelUp, pdTheInitClkBias, pdTheInitClkDrift));
      // System state covariance matrix
      CMatrix clSystemCovarianceP(VRWInitializeStateCovarianceP(pdTheInitPosStd_, pdTheInitVelStd_));
      // Trajectory epoch
      TrajectoryInfo stTrajInfo;
      TrajectoryInfo stKalmanTraj;
      // Save trajectories in containers
      vector<TrajectoryInfo> vTruthTrajectory;
      vector<TrajectoryInfo> vKalmanTrajectory;
      vector<VRWErrorInfo> vVRWError;
      // GPS Time from epoch before
      DOUBLE dPrevGPSTime = 0;
      // Perform Kalman filter operations
      for (vector<EpochInfo>::iterator itEpoch = vEpochData.begin(); itEpoch != vEpochData.end(); itEpoch++)
      {
         // Truth trajectory file provided
         if (psTheTruthFilePath_ != NULL)
         {
            AlignEpochs(clTruthFile, itEpoch, &stTrajInfo, &vEpochData);
         }
         else
         {
            FillTrajectory(&stTrajInfo, *pdTheTruthLatitude_, *pdTheTruthLongitude_, *pdTheTruthHeight_, itEpoch->dGPSTime);
         }

         // Kalman prediction
         if (itEpoch != vEpochData.begin())
         {
            clStateSpaceF = VRWStateSpaceMatrixF();
            clNoiseShapeG = VRWNoiseShapeMatrixG();
            KPredictionLoop(&clStateSpaceF, &clNoiseShapeG, &clSpectralDensityQ, &clGeoidSystemState, &clSystemCovarianceP, itEpoch->dGPSTime - dPrevGPSTime);
         }
         // Kalman update
         VRWUpdateLoop(&clGeoidSystemState, &clSystemCovarianceP, &(*itEpoch), *pdThePSRStd_, *pdThePSRRateStd_, bUsePseudorate_);
         // Save previous GPS Time for step size integration
         dPrevGPSTime = itEpoch->dGPSTime;
         // Save estimated trajectory into a container
         FillTrajectory(&stKalmanTraj, RadToDeg(clGeoidSystemState.GetComponent(0, 0)), RadToDeg(clGeoidSystemState.GetComponent(1, 0)), clGeoidSystemState.GetComponent(2, 0), itEpoch->dGPSTime,
                        clGeoidSystemState.GetComponent(4, 0), clGeoidSystemState.GetComponent(3, 0), clGeoidSystemState.GetComponent(5, 0));
         // Kalman error state structure
         VRWFillPositionError(&stTrajInfo, &stKalmanTraj, &clSystemCovarianceP, &vVRWError, itEpoch->cNumSV);
         vTruthTrajectory.push_back(stTrajInfo);
         vKalmanTrajectory.push_back(stKalmanTraj);
         // Output data on console window
         cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidSystemState.GetComponent(0, 0)) << ",";
         cout << fixed << setprecision(9) << RadToDeg(clGeoidSystemState.GetComponent(1, 0)) << ",";
         cout << fixed << setprecision(3) << clGeoidSystemState.GetComponent(2, 0) << ","; 
         cout << fixed << setprecision(3) << clGeoidSystemState.GetComponent(3, 0) << ",";
         cout << fixed << setprecision(3) << clGeoidSystemState.GetComponent(4, 0) << ",";
         cout << fixed << setprecision(3) << clGeoidSystemState.GetComponent(5, 0) << endl;
      }
      KalmanPosErrorsPlots(&vTruthTrajectory, &vKalmanTrajectory, *psTheFullFilePath_, "VRW");
      KalmanVelErrorsPlots(&vTruthTrajectory, &vKalmanTrajectory, *psTheFullFilePath_, "VRW");
      KalmanCrossCorrPlots(&vVRWError, *psTheFullFilePath_, "VRW");
      KalmanTrajectory(&vTruthTrajectory, *psTheFullFilePath_, "VRW", "TruthTrajectory");
      KalmanTrajectory(&vKalmanTrajectory, *psTheFullFilePath_, "VRW", "KalmanTrajectory");
   }
   // Close the file
   if (psTheTruthFilePath_ != NULL)
      clTruthFile.close();
}