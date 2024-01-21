//----------------------------------------------------------------------
// VRW.cpp
// Source file of VRW.h
// Source code for Position Random Walk system model
//----------------------------------------------------------------------
#include "..\inc\KPolar.h"

//-------------------------------------------------------------------
CMatrix KPolarStateSpaceMatrixF()
{
   CMatrix clStateF("State Matrix", KPOLAR_STATE_DIMENSION, KPOLAR_STATE_DIMENSION);

   // Error polar angle time derivative
   clStateF.SetComponent(4, 5, 1);
   // Error clock bias time derivative
   clStateF.SetComponent(6, 7, 1);

   return clStateF;
}

//-------------------------------------------------------------------
CMatrix KPolarNoiseShapeMatrixG()
{
   CMatrix clNoiseG("Noise Shaping Matrix", KPOLAR_STATE_DIMENSION, KPOLAR_NUM_NOISE_SOURCES);

   // Radial polar error noise
   clNoiseG.SetComponent(3, 0, 1);
   // Polar angle noise
   clNoiseG.SetComponent(4, 1, 1);
   // Angular velocity polar error noise
   clNoiseG.SetComponent(5, 2, 1);
   // Clock bias noise
   clNoiseG.SetComponent(6, 3, 1);
   // Clock drift noise
   clNoiseG.SetComponent(7, 4, 1);

   return clNoiseG;
}

//-------------------------------------------------------------------
CMatrix KPolarJacobian(CMatrix* pclPolar_)
{
   // Verify vector has the right coordinates
   if (pclPolar_->GetNumRows() != 3 || pclPolar_->GetNumCols() != 1)
      throw 24;

   DOUBLE dRadius = pclPolar_->GetComponent(0, 0);
   DOUBLE ddRadius;
   (abs(dRadius) > 0) ? ddRadius = dRadius / abs(dRadius) : ddRadius = 0;

   CMatrix clJacobian("Jacobian", 3, 3);
   // J(0, 0)
   clJacobian.SetComponent(0, 0, -ddRadius * pclPolar_->GetComponent(2, 0) * sin(pclPolar_->GetComponent(1, 0)));
   // J(0, 1)
   clJacobian.SetComponent(0, 1, -abs(dRadius) * pclPolar_->GetComponent(2, 0) * cos(pclPolar_->GetComponent(1, 0)));
   // J(0, 2)
   clJacobian.SetComponent(0, 2, -abs(dRadius) * sin(pclPolar_->GetComponent(1, 0)));
   // J(1, 0)
   clJacobian.SetComponent(1, 0, ddRadius * pclPolar_->GetComponent(2, 0) * cos(pclPolar_->GetComponent(1, 0)));
   // J(1, 1)
   clJacobian.SetComponent(1, 1, -abs(dRadius) * pclPolar_->GetComponent(2, 0) * sin(pclPolar_->GetComponent(1, 0)));
   // J(1, 2)
   clJacobian.SetComponent(1, 2, abs(dRadius) * cos(pclPolar_->GetComponent(1, 0)));
   // J(2, 0)
   clJacobian.SetComponent(2, 0, 0);
   // J(2, 1)
   clJacobian.SetComponent(2, 1, 0);
   // J(2, 2)
   clJacobian.SetComponent(2, 2, 0);

   return clJacobian;
}

//-------------------------------------------------------------------
CMatrix KPolarTransformH(CMatrix* pclPolar_)
{
   // Verify vector has the right coordinates
   if (pclPolar_->GetNumRows() != 3 || pclPolar_->GetNumCols() != 1)
      throw 25;

   CMatrix clTransformH("Transform Measurement Matrix", KPOLAR_STATE_DIMENSION, KPOLAR_STATE_DIMENSION);

   DOUBLE dRadius = pclPolar_->GetComponent(0, 0);
   DOUBLE ddRadius;
   (abs(dRadius) > 0) ? ddRadius = dRadius / abs(dRadius) : ddRadius = 0;

   // Row 1 settings
   clTransformH.SetComponent(0, 0, 1);
   clTransformH.SetComponent(0, 3, ddRadius * cos(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(0, 4, -abs(dRadius) * sin(pclPolar_->GetComponent(1, 0)));
   // Row 2 settings
   clTransformH.SetComponent(1, 1, 1);
   clTransformH.SetComponent(1, 3, ddRadius * sin(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(1, 4, abs(dRadius) * cos(pclPolar_->GetComponent(1, 0)));
   // Row 3 settings
   clTransformH.SetComponent(2, 2, 1);
   // Row 4 settings
   clTransformH.SetComponent(3, 3, -ddRadius * pclPolar_->GetComponent(2, 0) * sin(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(3, 4, -abs(dRadius) * pclPolar_->GetComponent(2, 0) * cos(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(3, 5, -abs(dRadius) * sin(pclPolar_->GetComponent(1, 0)));
   // Row 5 settings
   clTransformH.SetComponent(4, 3, ddRadius * pclPolar_->GetComponent(2, 0) * cos(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(4, 4, -abs(dRadius) * pclPolar_->GetComponent(2, 0) * sin(pclPolar_->GetComponent(1, 0)));
   clTransformH.SetComponent(4, 5, abs(dRadius) * cos(pclPolar_->GetComponent(1, 0)));
   // Row 7 settings
   clTransformH.SetComponent(6, 6, 1);
   // Row 8 settings
   clTransformH.SetComponent(7, 7, 1);

   return clTransformH;
}

//-------------------------------------------------------------------
CMatrix KPolarInitializeState(DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_,
                              DOUBLE* pdTheHeight_, DOUBLE* pdTheVelEast_,
                              DOUBLE* pdTheVelNorth_, DOUBLE* pdTheVelUp_,
                              DOUBLE* pdTheClkBias_, DOUBLE* pdTheClkDrift_,
                              DOUBLE* pdTheRadius_)
{
   CMatrix clSystemState("System State", KPOLAR_STATE_DIMENSION, 1);

   // Assign any reasonable arbitrary value to the radius [m]
   DOUBLE dRadius = *pdTheRadius_;
   // Find out initial polar angle based on the east and north velocities
   DOUBLE dTheta = atan2(*pdTheVelEast_ * (-1), *pdTheVelNorth_);
   // Compute angular velocity based on radius, polar angle and velocity
   DOUBLE dOmega;
   (fmod(dTheta, COORDTRANS_PI) < COORDTRANS_PI / 4) ? (dOmega = *pdTheVelNorth_ / (dRadius * cos(dTheta))) : (dOmega = *pdTheVelEast_ / (dRadius * sin(dTheta)));

   // Center of circle latitude
   clSystemState.SetComponent(0, 0, DegToRad(*pdTheLatitude_));
   // Center of circle longitude
   clSystemState.SetComponent(1, 0, DegToRad(*pdTheLongitude_));
   // Center of circle height
   clSystemState.SetComponent(2, 0, *pdTheHeight_);
   // Radius of circle
   clSystemState.SetComponent(3, 0, dRadius);
   // Polar angle of trajectory
   clSystemState.SetComponent(4, 0, dTheta);
   // Polar angular velocity of trajectory
   clSystemState.SetComponent(5, 0, dOmega);
   // Clock bias
   clSystemState.SetComponent(6, 0, *pdTheClkBias_);
   // Clock drift
   clSystemState.SetComponent(7, 0, *pdTheClkDrift_);

   return clSystemState;
}

//-------------------------------------------------------------------
CMatrix KPolarInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_, CMatrix* pclPolar_)
{
   CMatrix clStateCovP("State covariance", KPOLAR_STATE_DIMENSION, KPOLAR_STATE_DIMENSION);
   CMatrix clJacob(KPolarJacobian(pclPolar_));
   CMatrix clJacobTranspose(clJacob.Transpose());
   CMatrix clPolarVariance("Polar Variance", 3, 3);
   CMatrix clVelVariance("Diagonal", 3, 3);
   clVelVariance.SetIdentity();

   clVelVariance = clVelVariance * pow(*pdTheInitVelStd_, 2);
   clPolarVariance = clJacob.NumericPseudoInverse() * clVelVariance * clJacobTranspose.NumericPseudoInverse();

   // Set covariance for the position and clock states
   for (INT i = 0; i < clStateCovP.GetNumRows(); i++)
   {
      // Polar states covariance
      (i > 2 && i < 6) ? clStateCovP.SetComponent(i, i, pow(*pdTheInitVelStd_, 2)) : clStateCovP.SetComponent(i, i, pow(*pdTheInitPosStd_, 2));
   }

   // Set covariance for the relative polar coordinate states
   for (INT i = 0; i < clPolarVariance.GetNumRows(); i++)
   {
      for (INT j = 0; j < clPolarVariance.GetNumCols(); j++)
      {
         clStateCovP.SetComponent(i + 3, j + 3, clPolarVariance.GetComponent(i, j));
      }
   }

   return clStateCovP;
}

//-------------------------------------------------------------------
CMatrix KPolarSpectralDensityQ(DOUBLE* pdRWSpectralQ_, DOUBLE* pdClkBiasQ_, DOUBLE* pdClkDriftQ_, CMatrix* pclPolar_)
{
   CMatrix clQ("Spectral Densities", 5, 5);
   CMatrix clJacob(KPolarJacobian(pclPolar_));
   CMatrix clJacobTranspose(clJacob.Transpose());
   CMatrix clPolarQ("Polar Noise Density", 3, 3);
   CMatrix clVelQ("Diagonal", 3, 3);
   clVelQ.SetIdentity();

   clVelQ = clVelQ * pow(*pdRWSpectralQ_, 2);
   clPolarQ = clJacob.NumericPseudoInverse(1e-6) * clVelQ * clJacobTranspose.NumericPseudoInverse(1e-6);

   // Set spectral noise submatrix corresponding to radius, polar angle, angular velocity noise
   for (INT i = 0; i < clPolarQ.GetNumRows(); i++)
   {
      for (INT j = 0; j < clPolarQ.GetNumCols(); j++)
      {
         clQ.SetComponent(i, j, clPolarQ.GetComponent(i, j));
      }
   }
   clQ.SetComponent(3, 3, *pdClkBiasQ_);
   clQ.SetComponent(4, 4, *pdClkDriftQ_);

   return clQ;
}

//-------------------------------------------------------------------
void KPolarPredictionLoop(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
                          CMatrix* pclSysState_, CMatrix* pclSysCovP_, DOUBLE dStepSize_)
{
   // Local level frame state
   CMatrix clLLFState(*pclSysState_);
   // Set East, North and Up coordinates to origin
   clLLFState.SetComponent(0, 0, 0);
   clLLFState.SetComponent(1, 0, 0);
   clLLFState.SetComponent(2, 0, 0);

   // Transition state matrix
   CMatrix clTransitionPhi("State Transition Matrix", pclSysState_->GetNumRows(), pclSysState_->GetNumRows());
   // Process noise matrix
   CMatrix clProcessNoiseQ(KProcessNoiseQ(pclF_, pclG_, pclQ_, dStepSize_));

   // Find state transition matrix
   clTransitionPhi = (*pclF_ * dStepSize_).MatrixExponential2();

   // Next a-priori state prediction
   clLLFState = clTransitionPhi * clLLFState;
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   CMatrix clN(CalcN(&clLatSingleton));
   CMatrix clM(CalcM(&clLatSingleton));
   // Update values of latitude and longitude
   DOUBLE dLatCorrection = clLLFState.GetComponent(1, 0) / (clM.GetComponent(0, 0) + clLLFState.GetComponent(2, 0));
   DOUBLE dLngCorrection = clLLFState.GetComponent(0, 0) / ((clN.GetComponent(0, 0) + clLLFState.GetComponent(2, 0)) * cos(pclSysState_->GetComponent(0, 0)));

   // Local level frame updates, add corrections since origin corresponds to the curvilinear coordinates (latitude, longitude, height)
   pclSysState_->SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatCorrection);
   pclSysState_->SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngCorrection);
   pclSysState_->SetComponent(2, 0, pclSysState_->GetComponent(2, 0) + clLLFState.GetComponent(2, 0));
   // Apply updates
   for (INT i = 3; i < pclSysState_->GetNumRows(); i++)
   {
      pclSysState_->SetComponent(i, 0, clLLFState.GetComponent(i, 0));
   }
   pclSysState_->SetComponent(4, 0, fmod(pclSysState_->GetComponent(4, 0), 2 * COORDTRANS_PI));

   // Update state covariance matrix
   *pclSysCovP_ = clTransitionPhi * (*pclSysCovP_) * clTransitionPhi.Transpose() + clProcessNoiseQ;
}

//-------------------------------------------------------------------
void KPolarUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_,
                      EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_,
                      DOUBLE dStdPSRRate_, BOOLEANO bUsePseudorate_)
{
   CMatrix* pclDeltaPseudorate;
   CMatrix* pclPseudorateH;
   CMatrix* pclObsCovVelR;
   // Extract polar coordinates in a 3x1 vector
   CMatrix clPolar("LLF Polar Coordinates", 3, 1);
   clPolar.SetComponent(0, 0, pclSysState_->GetComponent(3, 0));
   clPolar.SetComponent(1, 0, pclSysState_->GetComponent(4, 0));
   clPolar.SetComponent(2, 0, pclSysState_->GetComponent(5, 0));
   // Curvilinear coordinates and Clock bias term
   CMatrix clGeoidVector("Curvilinear + Clock Bias", 4, 1);
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   CMatrix clN(CalcN(&clLatSingleton));
   CMatrix clM(CalcM(&clLatSingleton));
   DOUBLE dLatShift = abs(clPolar.GetComponent(0, 0)) * sin(clPolar.GetComponent(1, 0)) / (clM.GetComponent(0, 0) + pclSysState_->GetComponent(2, 0));
   DOUBLE dLngShift = abs(clPolar.GetComponent(0, 0)) * cos(clPolar.GetComponent(1, 0)) / ((clN.GetComponent(0, 0) + pclSysState_->GetComponent(2, 0)) * cos(pclSysState_->GetComponent(0, 0)));
   clGeoidVector.SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatShift);
   clGeoidVector.SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngShift);
   clGeoidVector.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   clGeoidVector.SetComponent(3, 0, pclSysState_->GetComponent(6, 0));
   // Curvilinear coordinates only
   CMatrix clGeoidCoordinates("Curvilinear coordinates", 3, 1);
   clGeoidCoordinates.SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatShift);
   clGeoidCoordinates.SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngShift);
   clGeoidCoordinates.SetComponent(2, 0, pclSysState_->GetComponent(2, 0));
   // Local Level Frame velocity vector
   CMatrix clLLFVelocity("LLF Velocity", 4, 1);
   DOUBLE dVelEast = -abs(clPolar.GetComponent(0, 0)) * clPolar.GetComponent(2, 0) * sin(clPolar.GetComponent(1, 0));
   DOUBLE dVelNorth = abs(clPolar.GetComponent(0, 0)) * clPolar.GetComponent(2, 0) * cos(clPolar.GetComponent(1, 0));
   clLLFVelocity.SetComponent(0, 0, dVelEast);
   clLLFVelocity.SetComponent(1, 0, dVelNorth);
   clLLFVelocity.SetComponent(2, 0, 0);
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
   clGeometryH = clGeometryH * KPolarTransformH(&clPolar);
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
   clLatSingleton.SetComponent(0, 0, pclSysState_->GetComponent(0, 0));
   clN = CalcN(&clLatSingleton);
   clM = CalcM(&clLatSingleton);
   // Correction values of latitude and longitude
   DOUBLE dLatCorrection = clDeltaState.GetComponent(1, 0) / (clM.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0));
   DOUBLE dLngCorrection = clDeltaState.GetComponent(0, 0) / ((clN.GetComponent(0, 0) + clDeltaState.GetComponent(2, 0)) * cos(pclSysState_->GetComponent(0, 0)));
   // Apply correction values to latitude, longitude, height, velocity east, velocity north, velocity up, clock bias, clock drift states
   pclSysState_->SetComponent(0, 0, pclSysState_->GetComponent(0, 0) + dLatCorrection);
   pclSysState_->SetComponent(1, 0, pclSysState_->GetComponent(1, 0) + dLngCorrection);
   pclSysState_->SetComponent(2, 0, pclSysState_->GetComponent(2, 0) + clDeltaState.GetComponent(2, 0));
   pclSysState_->SetComponent(3, 0, pclSysState_->GetComponent(3, 0) + clDeltaState.GetComponent(3, 0));
   pclSysState_->SetComponent(4, 0, fmod(pclSysState_->GetComponent(4, 0) + clDeltaState.GetComponent(4, 0), 2 * COORDTRANS_PI));
   pclSysState_->SetComponent(5, 0, pclSysState_->GetComponent(5, 0) + clDeltaState.GetComponent(5, 0));
   pclSysState_->SetComponent(6, 0, pclSysState_->GetComponent(6, 0) + clDeltaState.GetComponent(6, 0));
   pclSysState_->SetComponent(7, 0, pclSysState_->GetComponent(7, 0) + clDeltaState.GetComponent(7, 0));
   // Update state covariance matrix
   CMatrix clIdentity("Identity", pclSysState_->GetNumRows(), pclSysState_->GetNumRows());
   clIdentity.SetIdentity();
   CMatrix clTemp1(clIdentity - clKalmanGain * clGeometryH);
   *pclSysCovP_ = (clTemp1 * (*pclSysCovP_));
}

//----------------------------------------------------------------------------
void KPolarFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_,
                             CMatrix* pclStateCovP_, vector<VRWErrorInfo>* pvVRWKalmanError_,
                             CHAR cNumSV_, CMatrix* pclPolar_)
{
   VRWErrorInfo stPosError;

   // Transform state covariance matrix to LLF covariance matrix in trajectory points
   CMatrix clJacobi(KPolarTransformH(pclPolar_));
   CMatrix clStateCovP(clJacobi * (*pclStateCovP_) * clJacobi.Transpose());

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
   stPosError.dEastStd = sqrt(clStateCovP.GetComponent(0, 0));
   stPosError.dNorthStd = sqrt(clStateCovP.GetComponent(1, 1));
   stPosError.dHeightStd = sqrt(clStateCovP.GetComponent(2, 2));
   stPosError.dVelEastStd = sqrt(clStateCovP.GetComponent(3, 3));
   stPosError.dVelNorthStd = sqrt(clStateCovP.GetComponent(4, 4));
   stPosError.dVelUpStd = sqrt(clStateCovP.GetComponent(5, 5));
   // Compute correlation coefficients
   // East error cross-correlation
   stPosError.dENCorr = clStateCovP.GetComponent(0, 1) / (stPosError.dEastStd * stPosError.dNorthStd);
   stPosError.dEUCorr = clStateCovP.GetComponent(0, 2) / (stPosError.dEastStd * stPosError.dHeightStd);
   stPosError.dEVECorr = clStateCovP.GetComponent(0, 3) / (stPosError.dEastStd * sqrt(clStateCovP.GetComponent(3, 3)));
   stPosError.dEVNCorr = clStateCovP.GetComponent(0, 4) / (stPosError.dEastStd * sqrt(clStateCovP.GetComponent(4, 4)));
   stPosError.dEVUCorr = clStateCovP.GetComponent(0, 5) / (stPosError.dEastStd * sqrt(clStateCovP.GetComponent(5, 5)));
   stPosError.dEClkCorr = clStateCovP.GetComponent(0, 6) / (stPosError.dEastStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dEClkDriftCorr = clStateCovP.GetComponent(0, 7) / (stPosError.dEastStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // North error cross-correlation
   stPosError.dNUCorr = clStateCovP.GetComponent(1, 2) / (stPosError.dNorthStd * stPosError.dHeightStd);
   stPosError.dNVECorr = clStateCovP.GetComponent(1, 3) / (stPosError.dNorthStd * sqrt(clStateCovP.GetComponent(3, 3)));
   stPosError.dNVNCorr = clStateCovP.GetComponent(1, 4) / (stPosError.dNorthStd * sqrt(clStateCovP.GetComponent(4, 4)));
   stPosError.dNVUCorr = clStateCovP.GetComponent(1, 5) / (stPosError.dNorthStd * sqrt(clStateCovP.GetComponent(5, 5)));
   stPosError.dNClkCorr = clStateCovP.GetComponent(1, 6) / (stPosError.dNorthStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dNClkDriftCorr = clStateCovP.GetComponent(1, 7) / (stPosError.dNorthStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // Vertical error cross-correlation
   stPosError.dUVECorr = clStateCovP.GetComponent(2, 3) / (stPosError.dHeightStd * sqrt(clStateCovP.GetComponent(3, 3)));
   stPosError.dUVNCorr = clStateCovP.GetComponent(2, 4) / (stPosError.dHeightStd * sqrt(clStateCovP.GetComponent(4, 4)));
   stPosError.dUVUCorr = clStateCovP.GetComponent(2, 5) / (stPosError.dHeightStd * sqrt(clStateCovP.GetComponent(5, 5)));
   stPosError.dUClkCorr = clStateCovP.GetComponent(2, 6) / (stPosError.dHeightStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dUClkDriftCorr = clStateCovP.GetComponent(2, 7) / (stPosError.dHeightStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // East velocity error cross-correlation
   stPosError.dVEVNCorr = clStateCovP.GetComponent(3, 4) / (stPosError.dVelEastStd * sqrt(clStateCovP.GetComponent(4, 4)));
   stPosError.dVEVUCorr = clStateCovP.GetComponent(3, 5) / (stPosError.dVelEastStd * sqrt(clStateCovP.GetComponent(5, 5)));
   stPosError.dVEClkCorr = clStateCovP.GetComponent(3, 6) / (stPosError.dVelEastStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dVEClkDriftCorr = clStateCovP.GetComponent(3, 7) / (stPosError.dVelEastStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // North velocity error cross-correlation
   stPosError.dVNVUCorr = clStateCovP.GetComponent(4, 5) / (stPosError.dVelNorthStd * sqrt(clStateCovP.GetComponent(5, 5)));
   stPosError.dVNClkCorr = clStateCovP.GetComponent(4, 6) / (stPosError.dVelNorthStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dVNClkDriftCorr = clStateCovP.GetComponent(4, 7) / (stPosError.dVelNorthStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // Up velocity error cross-correlation
   stPosError.dVUClkCorr = clStateCovP.GetComponent(5, 6) / (stPosError.dVelUpStd * sqrt(clStateCovP.GetComponent(6, 6)));
   stPosError.dVUClkDriftCorr = clStateCovP.GetComponent(5, 7) / (stPosError.dVelUpStd * sqrt(clStateCovP.GetComponent(7, 7)));
   // Clock bias error cross-correlation
   stPosError.dClkClkDriftCorr = clStateCovP.GetComponent(6, 7) / (sqrt(clStateCovP.GetComponent(6, 6)) * sqrt(clStateCovP.GetComponent(7, 7)));

   pvVRWKalmanError_->push_back(stPosError);
}

//----------------------------------------------------------------------------
void KPolarSolution(DOUBLE* pdTheInitLatitude, DOUBLE* pdTheInitLongitude, DOUBLE* pdTheInitHeight,
                    DOUBLE* pdTheInitVelEast, DOUBLE* pdTheInitVelNorth, DOUBLE* pdTheInitVelUp,
                    DOUBLE* pdTheInitClkBias, DOUBLE* pdTheInitClkDrift, DOUBLE* pdTheTruthLatitude_,
                    DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_, string* psTheTruthFilePath_,
                    string* psTheFullFilePath_, DOUBLE* pdThePSRStd_, DOUBLE* pdThePSRRateStd_,
                    DOUBLE* pdThePRWSpectralQ_, DOUBLE* pdTheVRWSpectralQ_, DOUBLE* pdTheClkBiasQ_,
                    DOUBLE* pdTheClkDriftQ_, DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_,
                    BOOLEANO bTheUseDoppler, DOUBLE* pdTheInitRadius_)
{
   // Truth Trajectory file
   ifstream clTruthFile;
   // Container with satellite data
   vector<EpochInfo> vEpochData;

   PopulateSVEpochData(*psTheFullFilePath_, &vEpochData);
   SendSVEpochDatatoText(*psTheFullFilePath_, *psTheFullFilePath_);

   if (psTheTruthFilePath_ != NULL)
      clTruthFile.open(*psTheTruthFilePath_);

   cout << "\n\nPolar coordinates Kalman processing..." << endl;
   cout << "\n\nGPS Time, Latitude, Longitude, Height, VelocityEast, VelocityNorth, VelocityUp" << endl;
   // State space matrix
   CMatrix clStateSpaceF(KPolarStateSpaceMatrixF());
   // Noise shaping matrix
   CMatrix clNoiseShapeG(KPolarNoiseShapeMatrixG());
   // System state vector
   CMatrix clSystemState(KPolarInitializeState(pdTheInitLatitude, pdTheInitLongitude, pdTheInitHeight, pdTheInitVelEast, pdTheInitVelNorth, pdTheInitVelUp, pdTheInitClkBias, pdTheInitClkDrift, pdTheInitRadius_));
   CMatrix clPolarStates("Polar coordinate states", 3, 1);
   clPolarStates.SetComponent(0, 0, clSystemState.GetComponent(3, 0));
   clPolarStates.SetComponent(1, 0, clSystemState.GetComponent(4, 0));
   clPolarStates.SetComponent(2, 0, clSystemState.GetComponent(5, 0));
   // Spectral density matrix
   CMatrix clSpectralDensityQ(KPolarSpectralDensityQ(pdTheVRWSpectralQ_, pdTheClkBiasQ_, pdTheClkDriftQ_, &clPolarStates));
   // System state covariance matrix
   CMatrix clSystemCovarianceP(KPolarInitializeStateCovarianceP(pdTheInitPosStd_, pdTheInitVelStd_, &clPolarStates));
   // Geoid trajectory state
   CMatrix clGeoidState("Trajectory point", 6, 1);
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, clSystemState.GetComponent(0, 0));
   CMatrix clN(CalcN(&clLatSingleton));
   CMatrix clM(CalcM(&clLatSingleton));
   DOUBLE dNorthFactor = 0;
   DOUBLE dEastFactor = 0;
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
         // Get last epoch polar coordinates for computing spectral density matrix
         clPolarStates.SetComponent(0, 0, clSystemState.GetComponent(3, 0));
         clPolarStates.SetComponent(1, 0, clSystemState.GetComponent(4, 0));
         clPolarStates.SetComponent(2, 0, clSystemState.GetComponent(5, 0));
         // State space matrix
         clStateSpaceF = KPolarStateSpaceMatrixF();
         // Shape noise matrix
         clNoiseShapeG = KPolarNoiseShapeMatrixG();
         // Spectral density matrix
         clSpectralDensityQ = KPolarSpectralDensityQ(pdTheVRWSpectralQ_, pdTheClkBiasQ_, pdTheClkDriftQ_, &clPolarStates);
         // Compute next epoch a-priori estimates
         KPolarPredictionLoop(&clStateSpaceF, &clNoiseShapeG, &clSpectralDensityQ, &clSystemState, &clSystemCovarianceP, itEpoch->dGPSTime - dPrevGPSTime);

         //---------------------------------------------------------------------------------
         cout << "\n\nA-priori state" << endl;
         for (INT i = 0; i < clSystemState.GetNumRows(); i++)
         {
            cout << clSystemState.GetComponent(i, 0) << ",";
         }
         cout << "\n" << endl;
         //---------------------------------------------------------------------------------
      }
      // Kalman update
      KPolarUpdateLoop(&clSystemState, &clSystemCovarianceP, &(*itEpoch), *pdThePSRStd_, *pdThePSRRateStd_, bTheUseDoppler);

      //---------------------------------------------------------------------------------
      cout << "\n\nA-posteriori state" << endl;
      for (INT i = 0; i < clSystemState.GetNumRows(); i++)
      {
         cout << clSystemState.GetComponent(i, 0) << ",";
      }
      cout << "\n" << endl;
      //---------------------------------------------------------------------------------

      // Save previous GPS Time for step size integration
      dPrevGPSTime = itEpoch->dGPSTime;
      // Get last epoch polar coordinates for computing spectral density matrix
      clPolarStates.SetComponent(0, 0, clSystemState.GetComponent(3, 0));
      clPolarStates.SetComponent(1, 0, clSystemState.GetComponent(4, 0));
      clPolarStates.SetComponent(2, 0, clSystemState.GetComponent(5, 0));
      // Set Trajectory point estimates based on system state
      clLatSingleton.SetComponent(0, 0, clSystemState.GetComponent(0, 0));
      clN = CalcN(&clLatSingleton);
      clM = CalcM(&clLatSingleton);
      dNorthFactor = clM.GetComponent(0, 0) + clSystemState.GetComponent(2, 0);
      dEastFactor = (clN.GetComponent(0, 0) + clSystemState.GetComponent(2, 0)) * cos(clSystemState.GetComponent(0, 0));
      clGeoidState.SetComponent(0, 0, clSystemState.GetComponent(0, 0) + abs(clPolarStates.GetComponent(0, 0)) * sin(clPolarStates.GetComponent(1, 0)) / dNorthFactor);
      clGeoidState.SetComponent(1, 0, clSystemState.GetComponent(1, 0) + abs(clPolarStates.GetComponent(0, 0)) * cos(clPolarStates.GetComponent(1, 0)) / dEastFactor);
      clGeoidState.SetComponent(2, 0, clSystemState.GetComponent(2, 0));
      clGeoidState.SetComponent(3, 0, -abs(clPolarStates.GetComponent(0, 0)) * clPolarStates.GetComponent(2, 0) * sin(clPolarStates.GetComponent(1, 0)));
      clGeoidState.SetComponent(4, 0, abs(clPolarStates.GetComponent(0, 0)) * clPolarStates.GetComponent(2, 0) * cos(clPolarStates.GetComponent(1, 0)));
      clGeoidState.SetComponent(5, 0, 0);
      // Save estimated trajectory into a container
      FillTrajectory(&stKalmanTraj, RadToDeg(clGeoidState.GetComponent(0, 0)), RadToDeg(clGeoidState.GetComponent(1, 0)), clGeoidState.GetComponent(2, 0), 
                     itEpoch->dGPSTime, clGeoidState.GetComponent(4, 0), clGeoidState.GetComponent(3, 0), clGeoidState.GetComponent(5, 0));
      // Kalman error state structure
      KPolarFillPositionError(&stTrajInfo, &stKalmanTraj, &clSystemCovarianceP, &vVRWError, itEpoch->cNumSV, &clPolarStates);
      vTruthTrajectory.push_back(stTrajInfo);
      vKalmanTrajectory.push_back(stKalmanTraj);
      // Output data on console window
      cout << fixed << setprecision(3) << itEpoch->dGPSTime << ",";
      cout << fixed << setprecision(9) << RadToDeg(clGeoidState.GetComponent(0, 0)) << ",";
      cout << fixed << setprecision(9) << RadToDeg(clGeoidState.GetComponent(1, 0)) << ",";
      cout << fixed << setprecision(3) << clGeoidState.GetComponent(2, 0) << ",";
      cout << fixed << setprecision(3) << clGeoidState.GetComponent(3, 0) << ",";
      cout << fixed << setprecision(3) << clGeoidState.GetComponent(4, 0) << ",";
      cout << fixed << setprecision(3) << clGeoidState.GetComponent(5, 0) << endl;
   }
   KalmanPosErrorsPlots(&vTruthTrajectory, &vKalmanTrajectory, *psTheFullFilePath_, "POLAR");
   KalmanVelErrorsPlots(&vTruthTrajectory, &vKalmanTrajectory, *psTheFullFilePath_, "POLAR");
   KalmanCrossCorrPlots(&vVRWError, *psTheFullFilePath_, "POLAR");
   KalmanTrajectory(&vTruthTrajectory, *psTheFullFilePath_, "POLAR", "TruthTrajectory");
   KalmanTrajectory(&vKalmanTrajectory, *psTheFullFilePath_, "POLAR", "KalmanTrajectory");
   // Close the file
   if (psTheTruthFilePath_ != NULL)
      clTruthFile.close();
}