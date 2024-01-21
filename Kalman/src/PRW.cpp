//----------------------------------------------------------------------
// PRW.cpp
// Source file of PRW.h
// Source code for Position Random Walk system model
//----------------------------------------------------------------------
#include "..\inc\PRW.h"
#include <iomanip>

//-------------------------------------------------------------------
CMatrix PRWStateSpaceMatrixF()
{
   CMatrix clStateF("State Matrix", PRW_STATE_DIMENSION, PRW_STATE_DIMENSION);
   // Clock drift
   clStateF.SetComponent(3, 4, 1);

   return clStateF;
}

//-------------------------------------------------------------------
CMatrix PRWNoiseShapeMatrixG()
{
   CMatrix clNoiseG("Noise Shaping Matrix", PRW_STATE_DIMENSION, PRW_STATE_DIMENSION);
   clNoiseG.SetIdentity();

   return clNoiseG;
}

//-------------------------------------------------------------------
CMatrix PRWInitializeState(DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_,
                           DOUBLE* pdTheHeight_, DOUBLE* pdTheClkBias_,
                           DOUBLE* pdTheClkDrift_)
{
   CMatrix clSystemState("System State", PRW_STATE_DIMENSION, 1);

   clSystemState.SetComponent(0, 0, DegToRad(*pdTheLatitude_));
   clSystemState.SetComponent(1, 0, DegToRad(*pdTheLongitude_));
   clSystemState.SetComponent(2, 0, *pdTheHeight_);
   clSystemState.SetComponent(3, 0, *pdTheClkBias_);
   clSystemState.SetComponent(4, 0, *pdTheClkDrift_);

   return clSystemState;
}

//-------------------------------------------------------------------
CMatrix PRWInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_)
{
   CMatrix clStateCovP("State covariance", PRW_STATE_DIMENSION, PRW_STATE_DIMENSION);

   for (INT i = 0; i < clStateCovP.GetNumRows(); i++)
   {
      clStateCovP.SetComponent(i, i, pow(*pdTheInitPosStd_, 2));
   }

   return clStateCovP;
}

//-------------------------------------------------------------------
void PRWUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_, 
                   EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_)
{
   // Curvilinear coordinates and Clock bias term
   CMatrix clGeoidVector("Curvilinear + Clock Bias", 4, 1);
   clGeoidVector.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   clGeoidVector.SetComponent(1, 0, pclSysState_->GetComponent(1, 0));
   clGeoidVector.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   clGeoidVector.SetComponent(3, 0, pclSysState_->GetComponent(3, 0));
   // Curvilinear coordinates only
   CMatrix clGeoidCoordinates("Curvilinear coordinates", 3, 1);
   clGeoidCoordinates.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   clGeoidCoordinates.SetComponent(1, 0, pclSysState_->GetComponent(1, 0));
   clGeoidCoordinates.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   // ECEF coordinates and Clock bias term
   CMatrix clECEFVector("ECEF + Clock Bias", 4, 1);
   // Find current state ECEF coordinates
   GeoidToECEF(&clGeoidVector, &clECEFVector);
   // Find pseudo-range residuals with current epoch coordinates
   CMatrix clDeltaPseudorange(KPseudorangeResiduals(pstEpochInfo_, &clECEFVector));
   // Find pseudo-range geometry matrix
   CMatrix clPseudorangeH(KPseudorangeMatrixH(pstEpochInfo_, &clGeoidCoordinates));
   // Find system geometry matrix at current epoch
   CMatrix clGeometryH(KGeometryMatrixH(&clPseudorangeH));
   // Find observation covariance matrix
   CMatrix clObsCovR(KObsCovSubMatrixR(pstEpochInfo_, dStdPSR_));
   // Find Kalman Gain
   CMatrix clTemp(clGeometryH * (*pclSysCovP_) * clGeometryH.Transpose() + clObsCovR);
   CMatrix clKalmanGain(*pclSysCovP_ * clGeometryH.Transpose() * clTemp.NumericInverse2(1e-6));
   // Find correction vector
   CMatrix clDeltaState(clKalmanGain * clDeltaPseudorange);
   // Update system state
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   CMatrix clN(CalcN(&clLatSingleton));
   CMatrix clM(CalcM(&clLatSingleton));
   // Correction values of latitude and longitude
   DOUBLE dLatCorrection = clDeltaState.GetComponent(1, 0) / (clM.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0));
   DOUBLE dLngCorrection = clDeltaState.GetComponent(0, 0) / ((clN.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0)) * cos(pclSysState_->GetComponent(0, 0)));
   // Apply correction values to latitude, longitude, height, clock bias, clock drift states
   pclSysState_->SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatCorrection);
   pclSysState_->SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngCorrection);
   pclSysState_->SetComponent(2, 0, pclSysState_->GetComponent(2, 0) + clDeltaState.GetComponent(2, 0));
   pclSysState_->SetComponent(3, 0, pclSysState_->GetComponent(3, 0) + clDeltaState.GetComponent(3, 0));
   pclSysState_->SetComponent(4, 0, pclSysState_->GetComponent(4, 0) + clDeltaState.GetComponent(4, 0));
   // Update state covariance matrix
   CMatrix clIdentity("Identity", pclSysState_->GetNumRows(), pclSysState_->GetNumRows());
   clIdentity.SetIdentity();
   CMatrix clTemp1(clIdentity - clKalmanGain * clGeometryH);
   *pclSysCovP_ = (clTemp1 * (*pclSysCovP_));
}

//-------------------------------------------------------------------
void PRWFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_,
                          CMatrix* pclStateCovP_, vector<PRWErrorInfo>* pvPRWKalmanError_, CHAR cNumSV_)
{
   PRWErrorInfo stPosError;

   stPosError.dGPSTime = pstKalmanTraj_->dGPSTime;
   stPosError.cNumSV = cNumSV_;
   // Compute curvilinear errors
   DOUBLE dLatError = pstTruthTraj_->dLatitude - pstKalmanTraj_->dLatitude;
   DOUBLE dLngError = pstTruthTraj_->dLongitude - pstKalmanTraj_->dLongitude;
   // Compute meridial radii
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pstTruthTraj_->dLatitude);
   CMatrix N(CalcN(&clLatSingleton));
   CMatrix M(CalcM(&clLatSingleton));
   // Compute East, North, Up errors
   stPosError.dEastError = dLngError * (N.GetComponent(0, 0) + pstTruthTraj_->dHeight) * cos(clLatSingleton.GetComponent(0, 0));
   stPosError.dNorthError = dLatError * (M.GetComponent(0, 0) + pstTruthTraj_->dHeight);
   stPosError.dHeightError = pstTruthTraj_->dHeight - pstKalmanTraj_->dHeight;
   // Set East, North and Up estimated standard deviation
   stPosError.dEastStd = sqrt(pclStateCovP_->GetComponent(0, 0));
   stPosError.dNorthStd = sqrt(pclStateCovP_->GetComponent(1, 1));
   stPosError.dHeightStd = sqrt(pclStateCovP_->GetComponent(2, 2));
   // Compute correlation coefficients
   stPosError.dENCorr = pclStateCovP_->GetComponent(0, 1) / (stPosError.dEastStd * stPosError.dNorthStd);
   stPosError.dEUCorr = pclStateCovP_->GetComponent(0, 2) / (stPosError.dEastStd * stPosError.dHeightStd);
   stPosError.dEClkCorr = pclStateCovP_->GetComponent(0, 3) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dEClkDriftCorr = pclStateCovP_->GetComponent(0, 4) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dNUCorr = pclStateCovP_->GetComponent(1, 2) / (stPosError.dNorthStd * stPosError.dHeightStd);
   stPosError.dNClkCorr = pclStateCovP_->GetComponent(1, 3) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dNClkDriftCorr = pclStateCovP_->GetComponent(1, 4) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dUClkCorr = pclStateCovP_->GetComponent(2, 3) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dUClkDriftCorr = pclStateCovP_->GetComponent(2, 4) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dClkClkDriftCorr = pclStateCovP_->GetComponent(3, 4) / (sqrt(pclStateCovP_->GetComponent(3, 3)) * sqrt(pclStateCovP_->GetComponent(4, 4)));

   pvPRWKalmanError_->push_back(stPosError);
}