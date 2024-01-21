//----------------------------------------------------------------------
// Kmeasurements.cpp
// Source file of Kmeasurements.h
// Source code for Kalman measurements
//----------------------------------------------------------------------
#include "..\inc\Kmeasurements.h"

//-------------------------------------------------------------------
CMatrix KGeometryMatrixH(CMatrix* pclPseudorangeH_, CMatrix* pclPseudorateH_, INT iNumExtraDim_)
{
   // Output geometry matrix
   CMatrix clGeometryH(pclPseudorangeH_->AugmentDim(*pclPseudorangeH_, 0, 2 + iNumExtraDim_));

   // Augment matrix with velocity states
   if (pclPseudorateH_ != NULL)
   {
      clGeometryH.SetWarningEnable(FALSE);
      clGeometryH = pclPseudorangeH_->AugmentDim(*pclPseudorangeH_, pclPseudorateH_->GetNumRows(), 5);
   }

   // Fill in velocity, clock bias and clock drift components
   for (INT i = 0; i < clGeometryH.GetNumRows(); i++)
   {
      if (pclPseudorateH_ != NULL)
      {
         (i < pclPseudorangeH_->GetNumRows()) ? clGeometryH.SetComponent(i, 6, 1) : clGeometryH.SetComponent(i, 6, 0);
         // Fill the geometry matrix with pseudo-range rates measurements
         if (i >= pclPseudorangeH_->GetNumRows())
         {
            clGeometryH.SetComponent(i, 0, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 0));
            clGeometryH.SetComponent(i, 1, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 1));
            clGeometryH.SetComponent(i, 2, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 2));
            clGeometryH.SetComponent(i, 3, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 3));
            clGeometryH.SetComponent(i, 4, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 4));
            clGeometryH.SetComponent(i, 5, pclPseudorateH_->GetComponent(i - pclPseudorangeH_->GetNumRows(), 5));
            clGeometryH.SetComponent(i, 7, 1);
         }
      }
      else
      {
         clGeometryH.SetComponent(i, 3 + iNumExtraDim_, 1);
      }
   }

   clGeometryH.SetWarningEnable(TRUE);
   return clGeometryH;
}

//-------------------------------------------------------------------
CMatrix KPseudorangeMatrixH(EpochInfo* pstEpochData_, CMatrix* pclSysStateGeoid_)
{
   // Verify the point of expansion is a 3x1 vector
   // <latitude,longitude,height>
   if (pclSysStateGeoid_->GetNumCols() != 1 && pclSysStateGeoid_->GetNumRows() != 3)
      throw 19;

   // Hi = CTS
   CMatrix clAzimuthElevationC("Azimuth Elevation C", 1, 3);
   CMatrix clPseudorangeGeoidT("T Matrix", 3, 3);
   CMatrix clGeoidToLlfS(KGeoidToLLFS(pclSysStateGeoid_));
   CMatrix clDesignMatrixH("Design Matrix", (INT)pstEpochData_->cNumSV, 3);
   CMatrix clRowH("Row of Design Matrix", 1, 3);

   // Build the design matrix H
   for (INT i = 0; i < clDesignMatrixH.GetNumRows(); i++)
   {
      clAzimuthElevationC = KAzimuthElevationC(&pstEpochData_->vstSVInfo.at(i));
      clPseudorangeGeoidT = KPseudorangeGeoidT(&pstEpochData_->vstSVInfo.at(i), pclSysStateGeoid_);
      clRowH = clAzimuthElevationC * clPseudorangeGeoidT * clGeoidToLlfS;
      // Set components of design matrix H
      clDesignMatrixH.SetComponent(i, 0, clRowH.GetComponent(0, 0));
      clDesignMatrixH.SetComponent(i, 1, clRowH.GetComponent(0, 1));
      clDesignMatrixH.SetComponent(i, 2, clRowH.GetComponent(0, 2));
   }

   return clDesignMatrixH;
}

//-------------------------------------------------------------------
CMatrix KPseudorateMatrixH(EpochInfo* pstEpochData_, CMatrix* pclSysStateGeoid_)
{
   // Verify the point of expansion is a 3x1 vector
   // <latitude,longitude,height>
   if (pclSysStateGeoid_->GetNumCols() != 1 && pclSysStateGeoid_->GetNumRows() != 3)
      throw 20;

   // Hi = CTS
   CMatrix clAzimuthElevationC("Azimuth Elevation C", 1, 3);
   CMatrix clPseudorateGeoidT("T Matrix", 3, 6);
   CMatrix clDesignMatrixH("Design Matrix", (INT)pstEpochData_->cNumSV, 6);
   CMatrix clRowH("Row of Design Matrix", 1, 6);

   // Build the design matrix H
   for (INT i = 0; i < clDesignMatrixH.GetNumRows(); i++)
   {
      clAzimuthElevationC = KAzimuthElevationC(&pstEpochData_->vstSVInfo.at(i));
      clPseudorateGeoidT = KPseudorateGeoidTS(&pstEpochData_->vstSVInfo.at(i), pclSysStateGeoid_);
      clRowH = clAzimuthElevationC * clPseudorateGeoidT;
      // Set components of design matrix H
      clDesignMatrixH.SetComponent(i, 0, clRowH.GetComponent(0, 0));
      clDesignMatrixH.SetComponent(i, 1, clRowH.GetComponent(0, 1));
      clDesignMatrixH.SetComponent(i, 2, clRowH.GetComponent(0, 2));
      clDesignMatrixH.SetComponent(i, 3, clRowH.GetComponent(0, 3));
      clDesignMatrixH.SetComponent(i, 4, clRowH.GetComponent(0, 4));
      clDesignMatrixH.SetComponent(i, 5, clRowH.GetComponent(0, 5));
   }

   return clDesignMatrixH;
}

//-------------------------------------------------------------------
CMatrix KAzimuthElevationC(SVInfo* pstSVData_)
{
   CMatrix C("Azimuth Elevation C", 1, 3);
   C.SetComponent(0, 0, cos(DegToRad(pstSVData_->dElevation)) * sin(DegToRad(pstSVData_->dAzimuth)));
   C.SetComponent(0, 1, cos(DegToRad(pstSVData_->dElevation)) * cos(DegToRad(pstSVData_->dAzimuth)));
   C.SetComponent(0, 2, sin(DegToRad(pstSVData_->dElevation)));

   return C;
}

//-------------------------------------------------------------------
CMatrix KPseudorangeGeoidT(SVInfo* pstSVData_, CMatrix* pclSysStateGeoid_)
{
   // Verify the point of expansion is a 3x1 vector
   // <latitude,longitude,height>
   if (pclSysStateGeoid_->GetNumCols() != 1 && pclSysStateGeoid_->GetNumRows() != 3)
      throw 18;

   // Extract latitude from point of expansion to compute N and M
   CMatrix dLatSingleton("Latitude", 1, 1);
   dLatSingleton.SetComponent(0, 0, pclSysStateGeoid_->GetComponent(0, 0));
   CMatrix NSingleton(CalcN(&dLatSingleton));
   // Derivative of vertical meridian with respect to latitude
   DOUBLE dNdLatitude = CalcdNdLatitude(dLatSingleton.GetComponent(0, 0));
   // Eccentricity
   DOUBLE dEcc = CalcEccentricity();
   // Pseudo-range Geoid Matrix
   CMatrix T("T Matrix", 3, 3);
   // Set components
   // T(0,0)
   T.SetComponent(0, 0, 0);
   // T(0,1)
   T.SetComponent(0, 1, -pstSVData_->dXPosition * cos(pclSysStateGeoid_->GetComponent(1, 0)) - pstSVData_->dYPosition * sin(pclSysStateGeoid_->GetComponent(1, 0)));
   // T(0,2)
   T.SetComponent(0, 2, 0);
   // T(1,0)
   T.SetComponent(1, 0, -pstSVData_->dXPosition * cos(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dYPosition * cos(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dZPosition * sin(pclSysStateGeoid_->GetComponent(0, 0)) +
                         pow(dEcc, 2) * dNdLatitude * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0)) +
                         pow(dEcc, 2) * NSingleton.GetComponent(0, 0) * (pow(cos(pclSysStateGeoid_->GetComponent(0, 0)), 2) - pow(sin(pclSysStateGeoid_->GetComponent(0, 0)), 2)));
   // T(1,1)
   T.SetComponent(1, 1, pstSVData_->dXPosition * sin(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) -
                        pstSVData_->dYPosition * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)));
   // T(1,2)
   T.SetComponent(1, 2, 0);
   // T(2,0)
   T.SetComponent(2, 0, -pstSVData_->dXPosition * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dYPosition * sin(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) +
                         pstSVData_->dZPosition * cos(pclSysStateGeoid_->GetComponent(0, 0)) -
                         dNdLatitude * (1 - pow(dEcc, 2) * pow(sin(pclSysStateGeoid_->GetComponent(0, 0)), 2)) +
                         2 * pow(dEcc, 2) * NSingleton.GetComponent(0, 0) * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0)));
   // T(2,1)
   T.SetComponent(2, 1, -pstSVData_->dXPosition * cos(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) +
                         pstSVData_->dYPosition * cos(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)));
   // T(2,2)
   T.SetComponent(2, 2, -1);

   return T;
}

//-------------------------------------------------------------------
CMatrix KPseudorateGeoidTS(SVInfo* pstSVData_, CMatrix* pclSysStateGeoid_)
{
   // Verify the point of expansion is a 3x1 vector
   // <latitude,longitude,height>
   if (pclSysStateGeoid_->GetNumCols() != 1 && pclSysStateGeoid_->GetNumRows() != 3)
      throw 21;
   
   // Earth prime vertical and prime meridial radii
   CMatrix clN("Radius N", 1, 1);
   CMatrix clLatSingleton("Latitude", 1, 1);
   CMatrix clM("Radius M", 1, 1);
   CMatrix clTS("Geoid Velocity Matrix", 3, 6);
   
   // Compute prime vertical and prime meridial radii
   clLatSingleton.SetComponent(0, 0, pclSysStateGeoid_->GetComponent(0, 0));
   clN = CalcN(&clLatSingleton);
   clM = CalcM(&clLatSingleton);

   //# Set the velocity matrix
   // T(0, 0)
   clTS.SetComponent(0, 0, (-pstSVData_->dXVelocity * cos(pclSysStateGeoid_->GetComponent(1, 0)) - pstSVData_->dYVelocity * sin(pclSysStateGeoid_->GetComponent(1, 0))) /
                           ((clN.GetComponent(0, 0) + pclSysStateGeoid_->GetComponent(2, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0))));
   // T(0, 1)
   clTS.SetComponent(0, 1, 0);
   // T(0, 2)
   clTS.SetComponent(0, 2, 0);
   // T(0, 3)
   clTS.SetComponent(0, 3, -1);
   // T(0, 4)
   clTS.SetComponent(0, 4, 0);
   // T(0, 5)
   clTS.SetComponent(0, 5, 0);
   // T(1, 0)
   clTS.SetComponent(1, 0, (pstSVData_->dXVelocity * sin(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) -
                            pstSVData_->dYVelocity * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0))) /
                           ((clN.GetComponent(0, 0) + pclSysStateGeoid_->GetComponent(2, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0))));
   // T(1, 1)
   clTS.SetComponent(1, 1, (-pstSVData_->dXVelocity * cos(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)) -
                             pstSVData_->dYVelocity * cos(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) -
                             pstSVData_->dZVelocity * sin(pclSysStateGeoid_->GetComponent(0, 0))) / (clM.GetComponent(0, 0) + pclSysStateGeoid_->GetComponent(2, 0)));
   // T(1, 2)
   clTS.SetComponent(1, 2, 0);
   // T(1, 3)
   clTS.SetComponent(1, 3, 0);
   // T(1, 4)
   clTS.SetComponent(1, 4, -1);
   // T(1, 5)
   clTS.SetComponent(1, 5, 0);
   // T(2, 0)
   clTS.SetComponent(2, 0, (-pstSVData_->dXVelocity * cos(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) +
                             pstSVData_->dYVelocity * cos(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0))) /
                           ((clN.GetComponent(0, 0) + pclSysStateGeoid_->GetComponent(2, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0))));
   // T(2, 1)
   clTS.SetComponent(2, 1, (-pstSVData_->dXVelocity * sin(pclSysStateGeoid_->GetComponent(0, 0)) * cos(pclSysStateGeoid_->GetComponent(1, 0)) -
                             pstSVData_->dYVelocity * sin(pclSysStateGeoid_->GetComponent(0, 0)) * sin(pclSysStateGeoid_->GetComponent(1, 0)) +
                             pstSVData_->dZVelocity * cos(pclSysStateGeoid_->GetComponent(0, 0))) / (clM.GetComponent(0, 0) + pclSysStateGeoid_->GetComponent(2, 0)));
   // T(2, 2)
   clTS.SetComponent(2, 2, 0);
   // T(2, 3)
   clTS.SetComponent(2, 3, 0);
   // T(2, 4)
   clTS.SetComponent(2, 4, 0);
   // T(2, 5)
   clTS.SetComponent(2, 5, -1);

   return clTS;
}

//------------------------------------------------------------------------------------------------
CMatrix KGeoidToLLFS(CMatrix* pclSysStateGeoid_)
{
   // Verify the point of expansion is a 4x1 vector
   // <latitude,longitude,height,c*deltat>
   if (pclSysStateGeoid_->GetNumCols() != 1 && pclSysStateGeoid_->GetNumRows() != 3)
      throw 17;

   // Extract latitude from point of expansion to compute N and M
   CMatrix dLatSingleton("Latitude", 1, 1);
   dLatSingleton.SetComponent(0, 0, pclSysStateGeoid_->GetComponent(0, 0));
   DOUBLE dN = CalcN(&dLatSingleton).GetComponent(0, 0);
   DOUBLE dM = CalcM(&dLatSingleton).GetComponent(0, 0);
   CMatrix S("S Matrix", 3, 3);
   // Set components of matrix
   S.SetComponent(0, 0, 0);
   S.SetComponent(0, 1, 1 / (dM + pclSysStateGeoid_->GetComponent(2, 0)));
   S.SetComponent(0, 2, 0);
   S.SetComponent(1, 0, 1 / ((dN + pclSysStateGeoid_->GetComponent(2, 0)) * cos(pclSysStateGeoid_->GetComponent(0, 0))));
   S.SetComponent(1, 1, 0);
   S.SetComponent(2, 0, 0);
   S.SetComponent(2, 1, 0);
   S.SetComponent(2, 2, 1);

   return S;
}

//------------------------------------------------------------------------------------------------
CMatrix KResiudals(CMatrix* pclPseudoRangeResidual_, CMatrix* pclPseudoRateResidual_)
{
   CMatrix clResidual(*pclPseudoRangeResidual_);

   // Pseudo-range rate residuals present
   if (pclPseudoRateResidual_ != NULL)
   {
      clResidual.SetWarningEnable(FALSE);
      clResidual = pclPseudoRangeResidual_->AugmentDim(*pclPseudoRangeResidual_, pclPseudoRateResidual_->GetNumRows(), 0);
      // Populate residual vector with pseudorange rate residuals
      for (INT i = 0; i < pclPseudoRateResidual_->GetNumRows(); i++)
      {
         clResidual.SetComponent(pclPseudoRangeResidual_->GetNumRows() + i, 0, pclPseudoRateResidual_->GetComponent(i, 0));
      }
   }

   clResidual.SetWarningEnable(TRUE);
   return clResidual;
}

//------------------------------------------------------------------------------------------------
CMatrix KPseudorangeResiduals(EpochInfo* pstEpochData_, CMatrix* pclSysStateEcef_)
{
   // Verify the point of expansion is a 4x1 vector
   // <X,Y,Z,c*deltat>
   if (pclSysStateEcef_->GetNumCols() != 1 && pclSysStateEcef_->GetNumRows() != 4)
      throw 22;
   // Pseudo-range residual vector
   CMatrix PseudoRes("Pseudo-range residual vector", (INT)pstEpochData_->cNumSV, 1);
   DOUBLE dRange = 0;

   // Compute pseudo-range residual vector
   for (INT i = 0; i < PseudoRes.GetNumRows(); i++)
   {
      // Compute pseudo-range of point of expansion to satellite "i"
      dRange = sqrt(pow(pclSysStateEcef_->GetComponent(0, 0) - pstEpochData_->vstSVInfo.at(i).dXPosition, 2) +
                    pow(pclSysStateEcef_->GetComponent(1, 0) - pstEpochData_->vstSVInfo.at(i).dYPosition, 2) +
                    pow(pclSysStateEcef_->GetComponent(2, 0) - pstEpochData_->vstSVInfo.at(i).dZPosition, 2)) +
               pclSysStateEcef_->GetComponent(3, 0);
      // Set pseudo-range residual component
      PseudoRes.SetComponent(i, 0, pstEpochData_->vstSVInfo.at(i).dPSR - dRange);
   }

   return PseudoRes;
}

//------------------------------------------------------------------------------------------------
CMatrix KPseudorateResiduals(EpochInfo* pstEpochData_, CMatrix* pclSysStateEcef_)
{
   // Verify the point of expansion is a 7x1 vector
   // <X,Y,Z,VX,VY,VZ,cdeltat^dot>
   if (pclSysStateEcef_->GetNumCols() != 1 && pclSysStateEcef_->GetNumRows() != 7)
      throw 23;
   // Pseudo-range residual vector
   CMatrix clRateRes("Pseudorange rate residual vector", (INT)pstEpochData_->cNumSV, 1);
   DOUBLE dRange = 0;
   DOUBLE dVelDotPos = 0;
   DOUBLE dPseudoRate = 0;

   // Compute pseudo-range residual vector
   for (INT i = 0; i < clRateRes.GetNumRows(); i++)
   {
      // Compute pseudo-range of point of expansion to satellite "i"
      dRange = sqrt(pow(pclSysStateEcef_->GetComponent(0, 0) - pstEpochData_->vstSVInfo.at(i).dXPosition, 2) +
                    pow(pclSysStateEcef_->GetComponent(1, 0) - pstEpochData_->vstSVInfo.at(i).dYPosition, 2) +
                    pow(pclSysStateEcef_->GetComponent(2, 0) - pstEpochData_->vstSVInfo.at(i).dZPosition, 2));
      dVelDotPos = (pstEpochData_->vstSVInfo.at(i).dXVelocity - pclSysStateEcef_->GetComponent(3, 0)) * (pstEpochData_->vstSVInfo.at(i).dXPosition - pclSysStateEcef_->GetComponent(0, 0)) +
                   (pstEpochData_->vstSVInfo.at(i).dYVelocity - pclSysStateEcef_->GetComponent(4, 0)) * (pstEpochData_->vstSVInfo.at(i).dYPosition - pclSysStateEcef_->GetComponent(1, 0)) +
                   (pstEpochData_->vstSVInfo.at(i).dZVelocity - pclSysStateEcef_->GetComponent(5, 0)) * (pstEpochData_->vstSVInfo.at(i).dZPosition - pclSysStateEcef_->GetComponent(2, 0));
      dPseudoRate = dVelDotPos / dRange + pclSysStateEcef_->GetComponent(6, 0);
      // Set pseudo-range residual component
      clRateRes.SetComponent(i, 0, pstEpochData_->vstSVInfo.at(i).dPSRRate - dPseudoRate);
   }

   return clRateRes;
}

//------------------------------------------------------------------------------------------------
CMatrix KObsCovMatrixR(CMatrix* pclObsCovPseudorangeR_, CMatrix* pclObsCovPseudorateR_)
{
   CMatrix KObsCovMatrixR(*pclObsCovPseudorangeR_);

   if (pclObsCovPseudorateR_ != NULL)
   {
      KObsCovMatrixR.SetWarningEnable(FALSE);
      KObsCovMatrixR = pclObsCovPseudorangeR_->AugmentDim(*pclObsCovPseudorangeR_, pclObsCovPseudorateR_->GetNumRows(), pclObsCovPseudorateR_->GetNumCols());
      // Copy the pseudorange rate observation covariance matrix into 
      // the augmented observation covariance matrix
      for (INT i = 0; i < pclObsCovPseudorateR_->GetNumRows(); i++)
      {
         KObsCovMatrixR.SetComponent(pclObsCovPseudorangeR_->GetNumRows() + i, pclObsCovPseudorangeR_->GetNumRows() + i, pclObsCovPseudorateR_->GetComponent(i, i));
      }
   }

   KObsCovMatrixR.SetWarningEnable(TRUE);
   return KObsCovMatrixR;
}

//------------------------------------------------------------------------------------------------
CMatrix KObsCovSubMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_)
{
   CMatrix ObsCovarianceR("Observation Covariance Matrix", (INT)pstEpochData_->cNumSV, (INT)pstEpochData_->cNumSV);

   // Guess of observation covariance matrix
   for (INT i = 0; i < ObsCovarianceR.GetNumCols(); i++)
   {
      // Consider only satellites with elevation greater than 5 degrees
      if (pstEpochData_->vstSVInfo.at(i).dElevation > 5)
         ObsCovarianceR.SetComponent(i, i, pow(dStd_, 2) / pow(sin(DegToRad(pstEpochData_->vstSVInfo.at(i).dElevation)), 2));
      // Truncate elevation to a minimum of 5 degrees
      else
         ObsCovarianceR.SetComponent(i, i, pow(dStd_, 2) / pow(sin(DegToRad(5)), 2));
   }

   return ObsCovarianceR;
}

//-------------------------------------------------------------------
void KPredictionLoop(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
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

   // Update state covariance matrix
   *pclSysCovP_ = clTransitionPhi * (*pclSysCovP_) * clTransitionPhi.Transpose() + clProcessNoiseQ;
}

//-------------------------------------------------------------------
CMatrix KProcessNoiseQ(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
                       DOUBLE dStepSize_, ULONG ulNumSubIntervals_)
{
   DOUBLE dDeltaT = dStepSize_ / (DOUBLE)ulNumSubIntervals_;
   CMatrix clCurrScaledF("F*m*DeltaT", pclF_->GetNumRows(), pclF_->GetNumRows());
   CMatrix clPrevScaledF("F*(m-1)*DeltaT", pclF_->GetNumRows(), pclF_->GetNumRows());
   CMatrix clTopSubinterval("Top", pclF_->GetNumRows(), pclF_->GetNumRows());
   CMatrix clBottomSubinterval("Bottom", pclF_->GetNumRows(), pclF_->GetNumRows());
   CMatrix clQ("Process Noise", pclF_->GetNumRows(), pclF_->GetNumCols());

   // Perform the Riemann summation
   for (ULONG m = 1; m <= ulNumSubIntervals_; m++)
   {
      clCurrScaledF = *pclF_ * ((DOUBLE)m * dDeltaT);
      clPrevScaledF = *pclF_ * ((DOUBLE)(m - 1) * dDeltaT);
      clTopSubinterval = clCurrScaledF.MatrixExponential2() * (*pclG_) * (*pclQ_) * pclG_->Transpose() * clCurrScaledF.MatrixExponential2().Transpose() * 0.5;
      clBottomSubinterval = clPrevScaledF.MatrixExponential2() * (*pclG_) * (*pclQ_) * pclG_->Transpose() * clPrevScaledF.MatrixExponential2().Transpose() * 0.5;
      clQ = clQ + clTopSubinterval * dDeltaT + clBottomSubinterval * dDeltaT;
   }

   return clQ;
}

//-------------------------------------------------------------------
CMatrix KSpectralDensityQ(DOUBLE* pdRWSpectralQ_, DOUBLE* pdClkBiasQ_, DOUBLE* pdClkDriftQ_)
{
   CMatrix clQ("Spectral Densities", 5, 5);
   clQ.SetComponent(0, 0, pow(*pdRWSpectralQ_, 2));
   clQ.SetComponent(1, 1, pow(*pdRWSpectralQ_, 2));
   clQ.SetComponent(2, 2, pow(*pdRWSpectralQ_, 2));
   clQ.SetComponent(3, 3, *pdClkBiasQ_);
   clQ.SetComponent(4, 4, *pdClkDriftQ_);

   return clQ;
}