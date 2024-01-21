//----------------------------------------------------------------------
// VRW.cpp
// Source file of VRW.h
// Source code for Position Random Walk system model
//----------------------------------------------------------------------
#include "..\inc\VRW.h"

//-------------------------------------------------------------------
CMatrix VRWStateSpaceMatrixF()
{
   CMatrix clStateF("State Matrix", VRW_STATE_DIMENSION, VRW_STATE_DIMENSION);
   // Velocity East Error
   clStateF.SetComponent(0, 3, 1);
   // Velocity North Error
   clStateF.SetComponent(1, 4, 1);
   // Velocity Height Error
   clStateF.SetComponent(2, 5, 1);
   // Clock drift
   clStateF.SetComponent(6, 7, 1);

   return clStateF;
}

//-------------------------------------------------------------------
CMatrix VRWNoiseShapeMatrixG()
{
   CMatrix clNoiseG("Noise Shaping Matrix", VRW_STATE_DIMENSION, VRW_NUM_NOISE_SOURCES);
   // Velocity error east noise
   clNoiseG.SetComponent(3, 0, 1);
   // Velocity error north noise
   clNoiseG.SetComponent(4, 1, 1);
   // Velocity error up noise
   clNoiseG.SetComponent(5, 2, 1);
   // Clock bias noise
   clNoiseG.SetComponent(6, 3, 1);
   // Clock drift noise
   clNoiseG.SetComponent(7, 4, 1);

   return clNoiseG;
}

//-------------------------------------------------------------------
CMatrix VRWInitializeState(DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_,
                           DOUBLE* pdTheHeight_, DOUBLE* pdTheVelEast_,
                           DOUBLE* pdTheVelNorth_, DOUBLE* pdTheVelUp_,
                           DOUBLE* pdTheClkBias_, DOUBLE* pdTheClkDrift_)
{
   CMatrix clSystemState("System State", VRW_STATE_DIMENSION, 1);

   clSystemState.SetComponent(0, 0, DegToRad(*pdTheLatitude_));
   clSystemState.SetComponent(1, 0, DegToRad(*pdTheLongitude_));
   clSystemState.SetComponent(2, 0, *pdTheHeight_);
   clSystemState.SetComponent(3, 0, *pdTheVelEast_);
   clSystemState.SetComponent(4, 0, *pdTheVelNorth_);
   clSystemState.SetComponent(5, 0, *pdTheVelUp_);
   clSystemState.SetComponent(6, 0, *pdTheClkBias_);
   clSystemState.SetComponent(7, 0, *pdTheClkDrift_);

   return clSystemState;
}

//-------------------------------------------------------------------
CMatrix VRWInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_)
{
   CMatrix clStateCovP("State covariance", VRW_STATE_DIMENSION, VRW_STATE_DIMENSION);

   for (INT i = 0; i < clStateCovP.GetNumRows(); i++)
   {
      // Velocity states covariance
      (i > 2 && i < 6) ? clStateCovP.SetComponent(i, i, pow(*pdTheInitVelStd_, 2)) : clStateCovP.SetComponent(i, i, pow(*pdTheInitPosStd_, 2));
   }

   return clStateCovP;
}

//-------------------------------------------------------------------
void VRWUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_,
                   EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_, DOUBLE dStdPSRRate_, BOOLEANO bUsePseudorate_)
{
   CMatrix* pclDeltaPseudorate;
   CMatrix* pclPseudorateH;
   CMatrix* pclObsCovVelR;
   // Curvilinear coordinates and Clock bias term
   CMatrix clGeoidVector("Curvilinear + Clock Bias", 4, 1);
   clGeoidVector.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   clGeoidVector.SetComponent(1, 0, pclSysState_->GetComponent(1, 0));
   clGeoidVector.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   clGeoidVector.SetComponent(3, 0, pclSysState_->GetComponent(6, 0));
   // Curvilinear coordinates only
   CMatrix clGeoidCoordinates("Curvilinear coordinates", 3, 1);
   clGeoidCoordinates.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   clGeoidCoordinates.SetComponent(1, 0, pclSysState_->GetComponent(1, 0));
   clGeoidCoordinates.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   // Local Level Frame velocity vector
   CMatrix clLLFVelocity("LLF Velocity", 4, 1);
   clLLFVelocity.SetComponent(0, 0, pclSysState_->GetComponent(3, 0));
   clLLFVelocity.SetComponent(1, 0, pclSysState_->GetComponent(4, 0));
   clLLFVelocity.SetComponent(2, 0, pclSysState_->GetComponent(5, 0));
   clLLFVelocity.SetComponent(3, 0, pclSysState_->GetComponent(6, 0));
   // LLF to ECEF rotation matrix
   CMatrix clRel(ECEFToLLF(&clGeoidVector).Transpose());
   // ECEF coordinates and Clock bias term
   CMatrix clECEFVector("ECEF + Clock Bias", 4, 1);
   // Find current state ECEF coordinates
   GeoidToECEF(&clGeoidVector, &clECEFVector);
   // ECEF Velocity vector
   CMatrix clECEFVelVector(clRel * clLLFVelocity);
   CMatrix clSysStateECEF(clECEFVector.AugmentDim(clECEFVector, 3, 0));
   // Re-organize vector (X, Y, Z, Vx, Vy, Vz, clock-drift)
   clSysStateECEF.SetComponent(3, 0, clECEFVelVector.GetComponent(0, 0));
   clSysStateECEF.SetComponent(4, 0, clECEFVelVector.GetComponent(1, 0));
   clSysStateECEF.SetComponent(5, 0, clECEFVelVector.GetComponent(2, 0));
   clSysStateECEF.SetComponent(6, 0, pclSysState_->GetComponent(7, 0));
   // Find pseudo-range residuals with current epoch coordinates
   CMatrix clDeltaPseudorange(KPseudorangeResiduals(pstEpochInfo_, &clECEFVector));
   // Find pseudorange rate residuals with current epoch coordinates
   CMatrix clDeltaPseudorate(KPseudorateResiduals(pstEpochInfo_, &clSysStateECEF));
   // Find the complete measurement residuals
   (bUsePseudorate_) ? pclDeltaPseudorate = &clDeltaPseudorate : pclDeltaPseudorate = NULL;
   CMatrix clResiduals(KResiudals(&clDeltaPseudorange, pclDeltaPseudorate));
   // Find pseudo-range geometry matrix
   CMatrix clPseudorangeH(KPseudorangeMatrixH(pstEpochInfo_, &clGeoidCoordinates));
   // Find pseudorange rate geometry matrix
   CMatrix clPseudorateH(KPseudorateMatrixH(pstEpochInfo_, &clGeoidCoordinates));
   // Find system geometry matrix at current epoch
   (bUsePseudorate_) ? pclPseudorateH = &clPseudorateH : pclPseudorateH = NULL;
   CMatrix clGeometryH(KGeometryMatrixH(&clPseudorangeH, pclPseudorateH, 3));
   // Find pseudorange observation covariance matrix
   CMatrix clObsCovR1(KObsCovSubMatrixR(pstEpochInfo_, dStdPSR_));
   // Find pseudorange rate observation covaraince matrix
   CMatrix clObsCovR2(KObsCovSubMatrixR(pstEpochInfo_, dStdPSRRate_));
   // Find the observation covariance matrix
   (bUsePseudorate_) ? pclObsCovVelR = &clObsCovR2 : pclObsCovVelR = NULL;
   CMatrix clObsCovR(KObsCovMatrixR(&clObsCovR1, pclObsCovVelR));
   // Find Kalman Gain
   CMatrix clTemp(clGeometryH * (*pclSysCovP_) * clGeometryH.Transpose() + clObsCovR);
   CMatrix clKalmanGain(*pclSysCovP_ * clGeometryH.Transpose() * clTemp.NumericInverse2(1e-6));
   // Find correction vector
   CMatrix clDeltaState(clKalmanGain * clResiduals);
   // Update system state
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   CMatrix clN(CalcN(&clLatSingleton));
   CMatrix clM(CalcM(&clLatSingleton));
   // Correction values of latitude and longitude
   DOUBLE dLatCorrection = clDeltaState.GetComponent(1, 0) / (clM.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0));
   DOUBLE dLngCorrection = clDeltaState.GetComponent(0, 0) / ((clN.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0)) * cos(pclSysState_->GetComponent(0, 0)));
   // Apply correction values to latitude, longitude, height, velocity east, velocity north, velocity up, clock bias, clock drift states
   pclSysState_->SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatCorrection);
   pclSysState_->SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngCorrection);
   pclSysState_->SetComponent(2, 0, pclSysState_->GetComponent(2, 0) + clDeltaState.GetComponent(2, 0));
   pclSysState_->SetComponent(3, 0, pclSysState_->GetComponent(3, 0) + clDeltaState.GetComponent(3, 0));
   pclSysState_->SetComponent(4, 0, pclSysState_->GetComponent(4, 0) + clDeltaState.GetComponent(4, 0));
   pclSysState_->SetComponent(5, 0, pclSysState_->GetComponent(5, 0) + clDeltaState.GetComponent(5, 0));
   pclSysState_->SetComponent(6, 0, pclSysState_->GetComponent(6, 0) + clDeltaState.GetComponent(6, 0));
   pclSysState_->SetComponent(7, 0, pclSysState_->GetComponent(7, 0) + clDeltaState.GetComponent(7, 0));
   // Update state covariance matrix
   CMatrix clIdentity("Identity", pclSysState_->GetNumRows(), pclSysState_->GetNumRows());
   clIdentity.SetIdentity();
   CMatrix clTemp1(clIdentity - clKalmanGain * clGeometryH);
   *pclSysCovP_ = (clTemp1 * (*pclSysCovP_));
}

//-------------------------------------------------------------------
void VRWFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_,
                          CMatrix* pclStateCovP_, vector<VRWErrorInfo>* pvVRWKalmanError_, CHAR cNumSV_)
{
   VRWErrorInfo stPosError;

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
   stPosError.dVelEastError = pstTruthTraj_->dEastVel - pstKalmanTraj_->dEastVel;
   stPosError.dVelNorthError = pstTruthTraj_->dNorthVel - pstKalmanTraj_->dNorthVel;
   stPosError.dVelUpError = pstTruthTraj_->dHeightVel - pstKalmanTraj_->dHeightVel;
   // Set East, North and Up estimated standard deviation
   stPosError.dEastStd = sqrt(pclStateCovP_->GetComponent(0, 0));
   stPosError.dNorthStd = sqrt(pclStateCovP_->GetComponent(1, 1));
   stPosError.dHeightStd = sqrt(pclStateCovP_->GetComponent(2, 2));
   stPosError.dVelEastStd = sqrt(pclStateCovP_->GetComponent(3, 3));
   stPosError.dVelNorthStd = sqrt(pclStateCovP_->GetComponent(4, 4));
   stPosError.dVelUpStd = sqrt(pclStateCovP_->GetComponent(5, 5));
   // Compute correlation coefficients
   // East error cross-correlation
   stPosError.dENCorr = pclStateCovP_->GetComponent(0, 1) / (stPosError.dEastStd * stPosError.dNorthStd);
   stPosError.dEUCorr = pclStateCovP_->GetComponent(0, 2) / (stPosError.dEastStd * stPosError.dHeightStd);
   stPosError.dEVECorr = pclStateCovP_->GetComponent(0, 3) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dEVNCorr = pclStateCovP_->GetComponent(0, 4) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dEVUCorr = pclStateCovP_->GetComponent(0, 5) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(5, 5)));
   stPosError.dEClkCorr = pclStateCovP_->GetComponent(0, 6) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dEClkDriftCorr = pclStateCovP_->GetComponent(0, 7) / (stPosError.dEastStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // North error cross-correlation
   stPosError.dNUCorr = pclStateCovP_->GetComponent(1, 2) / (stPosError.dNorthStd * stPosError.dHeightStd);
   stPosError.dNVECorr = pclStateCovP_->GetComponent(1, 3) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dNVNCorr = pclStateCovP_->GetComponent(1, 4) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dNVUCorr = pclStateCovP_->GetComponent(1, 5) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(5, 5)));
   stPosError.dNClkCorr = pclStateCovP_->GetComponent(1, 6) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dNClkDriftCorr = pclStateCovP_->GetComponent(1, 7) / (stPosError.dNorthStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // Vertical error cross-correlation
   stPosError.dUVECorr = pclStateCovP_->GetComponent(2, 3) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(3, 3)));
   stPosError.dUVNCorr = pclStateCovP_->GetComponent(2, 4) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dUVUCorr = pclStateCovP_->GetComponent(2, 5) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(5, 5)));
   stPosError.dUClkCorr = pclStateCovP_->GetComponent(2, 6) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dUClkDriftCorr = pclStateCovP_->GetComponent(2, 7) / (stPosError.dHeightStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // East velocity error cross-correlation
   stPosError.dVEVNCorr = pclStateCovP_->GetComponent(3, 4) / (stPosError.dVelEastStd * sqrt(pclStateCovP_->GetComponent(4, 4)));
   stPosError.dVEVUCorr = pclStateCovP_->GetComponent(3, 5) / (stPosError.dVelEastStd * sqrt(pclStateCovP_->GetComponent(5, 5)));
   stPosError.dVEClkCorr = pclStateCovP_->GetComponent(3, 6) / (stPosError.dVelEastStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dVEClkDriftCorr = pclStateCovP_->GetComponent(3, 7) / (stPosError.dVelEastStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // North velocity error cross-correlation
   stPosError.dVNVUCorr = pclStateCovP_->GetComponent(4, 5) / (stPosError.dVelNorthStd * sqrt(pclStateCovP_->GetComponent(5, 5)));
   stPosError.dVNClkCorr = pclStateCovP_->GetComponent(4, 6) / (stPosError.dVelNorthStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dVNClkDriftCorr = pclStateCovP_->GetComponent(4, 7) / (stPosError.dVelNorthStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // Up velocity error cross-correlation
   stPosError.dVUClkCorr = pclStateCovP_->GetComponent(5, 6) / (stPosError.dVelUpStd * sqrt(pclStateCovP_->GetComponent(6, 6)));
   stPosError.dVUClkDriftCorr = pclStateCovP_->GetComponent(5, 7) / (stPosError.dVelUpStd * sqrt(pclStateCovP_->GetComponent(7, 7)));
   // Clock bias error cross-correlation
   stPosError.dClkClkDriftCorr = pclStateCovP_->GetComponent(6, 7) / (sqrt(pclStateCovP_->GetComponent(6, 6)) * sqrt(pclStateCovP_->GetComponent(7, 7)));

   pvVRWKalmanError_->push_back(stPosError);
}