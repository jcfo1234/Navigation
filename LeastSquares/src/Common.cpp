//----------------------------------------------------------------------
// Common.cpp
// Source file of Common.h
// Source code common to ECEF and curvilinear least squares computations
//----------------------------------------------------------------------
#include "..\inc\Common.h"

//------------------------------------------------------------------------------------------------
CMatrix CalculatePseudorangeResiduals(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionEcef_)
{
   // Verify the point of expansion is a 4x1 vector
   // <X,Y,Z,c*deltat>
   if (pclPointExpansionEcef_->GetNumCols() != 1 && pclPointExpansionEcef_->GetNumRows() != 4)
      throw 5;
   // Pseudo-range residual vector
   CMatrix PseudoRes("Pseudo-range residual vector", (INT)pstEpochData_->cNumSV, 1);
   DOUBLE dRange = 0;

   // Compute pseudo-range residual vector
   for (INT i = 0; i < PseudoRes.GetNumRows(); i++)
   {
      // Compute pseudo-range of point of expansion to satellite "i"
      dRange = sqrt(pow(pclPointExpansionEcef_->GetComponent(0, 0) - pstEpochData_->vstSVInfo.at(i).dXPosition, 2) +
                    pow(pclPointExpansionEcef_->GetComponent(1, 0) - pstEpochData_->vstSVInfo.at(i).dYPosition, 2) +
                    pow(pclPointExpansionEcef_->GetComponent(2, 0) - pstEpochData_->vstSVInfo.at(i).dZPosition, 2)) +
               pclPointExpansionEcef_->GetComponent(3, 0);
      // Set pseudo-range residual component
      PseudoRes.SetComponent(i, 0, pstEpochData_->vstSVInfo.at(i).dPSR - dRange);
   }

   return PseudoRes;
}

//------------------------------------------------------------------------------------------------
CMatrix ElevObsCovariancMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_)
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

//------------------------------------------------------------------------------------------------
CMatrix ConstObsCovariancMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_)
{
   CMatrix ObsCovarianceR("Observation Covariance Matrix", (INT) pstEpochData_->cNumSV, (INT) pstEpochData_->cNumSV);

   // Guess of observation covariance matrix
   for (INT i = 0; i < ObsCovarianceR.GetNumCols(); i++)
   {
      ObsCovarianceR.SetComponent(i, i, pow(dStd_, 2));
   }

   return ObsCovarianceR;
}

//------------------------------------------------------------------------------------------------
CMatrix ErrorStateCovarianceP(CMatrix* pclObsCovarianceMatrixR_, CMatrix* pclDesignMatrixH_)
{
   CMatrix clDesignMatrixTransposeH(pclDesignMatrixH_->Transpose());
   CMatrix clObsCovMatrixInverseRINV(pclObsCovarianceMatrixR_->NumericInverse());
   CMatrix clTemp(clDesignMatrixTransposeH * clObsCovMatrixInverseRINV * (*pclDesignMatrixH_));
   return clTemp.NumericInverse2();
}

//------------------------------------------------------------------------------------------------
CMatrix ErrorCorrectionVector(CMatrix* pclObsCovarianceMatrixR_, CMatrix* pclDesignMatrixH_, CMatrix* pclMeasurementResidual_)
{
   // Compute the error correction vector
   CMatrix clDesignMatrixTransposeH(pclDesignMatrixH_->Transpose());
   CMatrix clObsCovMatrixInverseRINV(pclObsCovarianceMatrixR_->NumericInverse());
   CMatrix clTemp(clDesignMatrixTransposeH * clObsCovMatrixInverseRINV * (*pclDesignMatrixH_));
   CMatrix cldeltax(clTemp.NumericInverse2() * clDesignMatrixTransposeH * clObsCovMatrixInverseRINV * (*pclMeasurementResidual_));

   return cldeltax;
}

//------------------------------------------------------------------------------------------------
void SetPointOfExpansionGeoid(DOUBLE dLat_, DOUBLE dLng_, DOUBLE dHgt_, DOUBLE dTshift_, CMatrix* pclPointExpansionGeoid_)
{
   // Verify size of vector
   if (pclPointExpansionGeoid_->GetNumCols() != 1 && pclPointExpansionGeoid_->GetNumRows() != 4)
      throw 12;

   // Initialize point of expansion
   pclPointExpansionGeoid_->SetComponent(0, 0, DegToRad(dLat_));
   pclPointExpansionGeoid_->SetComponent(1, 0, DegToRad(dLng_));
   pclPointExpansionGeoid_->SetComponent(2, 0, dHgt_);
   pclPointExpansionGeoid_->SetComponent(3, 0, dTshift_);
}

//------------------------------------------------------------------------------------------------
void SetPointOfExpansionEcef(DOUBLE dLat_, DOUBLE dLng_, DOUBLE dHgt_, DOUBLE dTshift_, CMatrix* pclPointExpansionEcef_)
{
   // Verify size of vector
   if (pclPointExpansionEcef_->GetNumCols() != 1 && pclPointExpansionEcef_->GetNumRows() != 4)
      throw 13;

   // Initialize point of expansion
   CMatrix clPointExpansionGeoid("Geoid", 4, 1);
   clPointExpansionGeoid.SetComponent(0, 0, DegToRad(dLat_));
   clPointExpansionGeoid.SetComponent(1, 0, DegToRad(dLng_));
   clPointExpansionGeoid.SetComponent(2, 0, dHgt_);
   clPointExpansionGeoid.SetComponent(3, 0, dTshift_);

   GeoidToECEF(&clPointExpansionGeoid, pclPointExpansionEcef_);
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeDOP(CMatrix* pclTransformationMatrix_, CMatrix* pclErrorCovMatrix_)
{
   return ( (*pclTransformationMatrix_) * (*pclErrorCovMatrix_) * pclTransformationMatrix_->Transpose() );
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeResidualCovarianceMatrix(CMatrix* pclDesignMatrixH_, CMatrix* pclObsCovarianceMatrixR_)
{
   // (H_Transpose * R_Inverse * H)
   CMatrix clTemp(pclDesignMatrixH_->Transpose() * pclObsCovarianceMatrixR_->NumericInverse() * (*pclDesignMatrixH_));
   // R - H * (H_Transpose * R_Inverse * H)_Inverse * H_Transpose
   CMatrix clCr((*pclObsCovarianceMatrixR_) - (*pclDesignMatrixH_) * clTemp.NumericInverse2() * pclDesignMatrixH_->Transpose());
   return clCr;
}

//------------------------------------------------------------------------------------------------
BOOLEANO BlunderDetection(EpochInfo* pstEpochData_, Distributions* pclProbDist_, CMatrix* pclPseudoRes_, 
                          CMatrix* pclResCovMatrixCr_, DOUBLE dAlpha_, INT iNumBlunders_)
{
   DOUBLE dMaxBlunder = 0;
   DOUBLE dTestValue = pclProbDist_->GetNormalAbscise(1.0 - dAlpha_ / 2.0) / (DOUBLE) iNumBlunders_;
   INT iBlunderIndex = 0;
   BOOLEANO bBlunderPass = TRUE;

   // Test each residual
   for (INT i = 0; i < pclPseudoRes_->GetNumRows(); i++)
   {
      // Record the maximum blunder
      if ( abs(pclPseudoRes_->GetComponent(i, 0) / sqrt(pclResCovMatrixCr_->GetComponent(i, i))) > dMaxBlunder )
      {
         dMaxBlunder = pclPseudoRes_->GetComponent(i, 0) / sqrt(pclResCovMatrixCr_->GetComponent(i, i));
         iBlunderIndex = i;
      }
   }

   // Remove item with highest blunder if blunder test fails
   if (dMaxBlunder > dTestValue)
   {
      // If blunder is in pseudo-range measurements
      if (iBlunderIndex < pstEpochData_->cNumSV)
      {
         pstEpochData_->vstSVInfo.erase(pstEpochData_->vstSVInfo.begin() + iBlunderIndex);
         pstEpochData_->cNumSV = pstEpochData_->cNumSV - 1;
      }
      bBlunderPass = FALSE;
   }

   return bBlunderPass;
}

//------------------------------------------------------------------------------------------------
void FillCoordinates(CMatrix* pclCoord4D_, CMatrix* pclCoord3D_)
{
	// Check for vector dimensions
	if ((pclCoord4D_->GetNumCols() != 1 && pclCoord4D_->GetNumRows() != 4) ||
		(pclCoord3D_->GetNumCols() != 1 && pclCoord3D_->GetNumRows() != 3))
	{
		throw 16;
	}
	// Set ECEF coordinates 3D vector
	for (INT i = 0; i < pclCoord3D_->GetNumRows(); i++)
		pclCoord3D_->SetComponent(i, 0, pclCoord4D_->GetComponent(i, 0));
}

//------------------------------------------------------------------------------------------------
void FillPositionErrorVector(vector<ErrorInfo>* pvPosError_, CMatrix* pclGeoidCoordinates_, 
                             TrajectoryInfo* pstTrajEpoch_, CMatrix* pclErrCovP_, 
                             CMatrix* pclTransformation_, CHAR cNumSV_, DOUBLE dMeasStd_)
{
	ErrorInfo stPosError;

	stPosError.dGPSTime = pstTrajEpoch_->dGPSTime;
   stPosError.cNumSV = cNumSV_;
	// Compute curvilinear errors
	DOUBLE dLatError = pstTrajEpoch_->dLatitude - pclGeoidCoordinates_->GetComponent(0, 0);
	DOUBLE dLngError = pstTrajEpoch_->dLongitude - pclGeoidCoordinates_->GetComponent(1, 0);
	// Compute meridial radii
	CMatrix clLatSingleton("Latitude", 1, 1);
	clLatSingleton.SetComponent(0, 0, pstTrajEpoch_->dLatitude);
	CMatrix N(CalcN(&clLatSingleton));
	CMatrix M(CalcM(&clLatSingleton));
	// Compute East, North, Up errors
	stPosError.dEastError = dLngError * (N.GetComponent(0, 0) + pstTrajEpoch_->dHeight) * cos(clLatSingleton.GetComponent(0, 0));
	stPosError.dNorthError = dLatError * (M.GetComponent(0, 0) + pstTrajEpoch_->dHeight);
	stPosError.dHeightError = pstTrajEpoch_->dHeight - pclGeoidCoordinates_->GetComponent(2, 0);
	// Rotation to Local level frame
	CMatrix clDOP(ComputeDOP(pclTransformation_, pclErrCovP_));
	// Set East, North and Up estimated standard deviation
	stPosError.dEastStd = sqrt(clDOP.GetComponent(0, 0));
	stPosError.dNorthStd = sqrt(clDOP.GetComponent(1, 1));
	stPosError.dHeightStd = sqrt(clDOP.GetComponent(2, 2));
	// Compute Dilution of precision
	stPosError.dHDop = sqrt(clDOP.GetComponent(0, 0) + clDOP.GetComponent(1, 1)) / dMeasStd_;
	stPosError.dVDop = stPosError.dHeightStd / dMeasStd_;
	stPosError.dPDop = sqrt(clDOP.GetComponent(0, 0) + clDOP.GetComponent(1, 1) + clDOP.GetComponent(2, 2)) / dMeasStd_;
   // Compute correlation coefficients
   stPosError.dENCorr = clDOP.GetComponent(0, 1) / (stPosError.dEastStd * stPosError.dNorthStd);
   stPosError.dEUCorr = clDOP.GetComponent(0, 2) / (stPosError.dEastStd * stPosError.dHeightStd);
   stPosError.dEClkCorr = clDOP.GetComponent(0, 3) / (stPosError.dEastStd * sqrt(clDOP.GetComponent(3, 3)));
   stPosError.dNUCorr = clDOP.GetComponent(1, 2) / (stPosError.dNorthStd * stPosError.dHeightStd);
   stPosError.dNClkCorr = clDOP.GetComponent(1, 3) / (stPosError.dNorthStd * sqrt(clDOP.GetComponent(3, 3)));
   stPosError.dUClkCorr = clDOP.GetComponent(2, 3) / (stPosError.dHeightStd * sqrt(clDOP.GetComponent(3, 3)));

	pvPosError_->push_back(stPosError);
}

//------------------------------------------------------------------------------------------------
void GlobalCovarianceTest(CMatrix* pclDeltaPseudoRange_, CMatrix* pclObsCovMatrix_, 
                          ReliablityTest* pstRelTest_, Distributions* pclProbDist_,
                          DOUBLE dAlpha_)
{
   INT iDegreesOfFreedom = pclDeltaPseudoRange_->GetNumRows() - DATATYPEDEFINITIONS_NUMSTATES;
   CMatrix clVar("Variance", 1, 1);
   clVar = pclDeltaPseudoRange_->Transpose() * pclObsCovMatrix_->NumericInverse() * (*pclDeltaPseudoRange_);
   DOUBLE dHighBound = pclProbDist_->GetChiSquareAbscise((ULONG)iDegreesOfFreedom, dAlpha_/2);
   DOUBLE dLowBound = pclProbDist_->GetChiSquareAbscise((ULONG)iDegreesOfFreedom, 1.0 - dAlpha_ / 2);
   pstRelTest_->dGlobVar = clVar.GetComponent(0, 0);
   // Set limits for global variance test
   pstRelTest_->dHighBound = dHighBound;
   pstRelTest_->dLowBound = dLowBound;
   (pstRelTest_->dGlobVar > dHighBound) ? pstRelTest_->bGlobTest = FALSE : (pstRelTest_->dGlobVar < dLowBound) ? pstRelTest_->bGlobTest = FALSE : pstRelTest_->bGlobTest = TRUE;
}

//------------------------------------------------------------------------------------------------
void ReliabilityTest(EpochInfo* pstEpochData_, CMatrix* pclObsCovMatrixR_,
                     CMatrix* pclResCovMatrixCr_, CMatrix* pclDesignMatrixH_,
                     CMatrix* pclCoordinates_, BOOLEANO bECEFCoord_, CMatrix* pclErrorStateCovP_,
                     CMatrix* pclDeltaPseudoRange_, Distributions* pclProbDist_,
                     DOUBLE dAlpha_, DOUBLE dBeta_, ReliablityTest* pstRelTest_)
{
   // Non-centrality parameter
   DOUBLE dDelta = pclProbDist_->GetNormalAbscise(1.0 - dAlpha_ / 2.0) + pclProbDist_->GetNormalAbscise(1 - dBeta_);
   ReliabilityInfo stRelInfo;
   // Set current epoch time
   pstRelTest_->dGPSTime = pstEpochData_->dGPSTime;
   // Set current epoch number of satellites
   pstRelTest_->cNumSV = pstEpochData_->cNumSV;
   // Compute P * H^T * R^-1
   CMatrix clTempMatrix("Temp", DATATYPEDEFINITIONS_NUMSTATES, pclObsCovMatrixR_->GetNumRows());
   clTempMatrix = (*pclErrorStateCovP_) * pclDesignMatrixH_->Transpose() * pclObsCovMatrixR_->NumericInverse();
   // Compute Rotation matrices from ECEF to LLF
   CMatrix clMyCoordinates("Coordinates", DATATYPEDEFINITIONS_NUMSTATES - 1, 1);
   FillCoordinates(pclCoordinates_, &clMyCoordinates);
   CMatrix clMyCoordinates2(clMyCoordinates);
   
   // Transform ECEF to Curvilinear coordinates
   if (bECEFCoord_)
   {
      ECEFToGeoid(&clMyCoordinates, &clMyCoordinates2);
   }
   CMatrix clTempVector(clMyCoordinates2.AugmentDim(clMyCoordinates2, 1, 0));
   CMatrix clRle(ECEFToLLF(&clTempVector));
   // External Reliability vectors
   CMatrix clECEFExtern("ECEF External Reliability", DATATYPEDEFINITIONS_NUMSTATES, 1);
   CMatrix clLLFExtern("LLF External Reliability", DATATYPEDEFINITIONS_NUMSTATES, 1);
   CMatrix clExternVector("External Reliability Vector", pclObsCovMatrixR_->GetNumRows(), 1);

   for (INT i = 0; i < pclObsCovMatrixR_->GetNumRows(); i++)
   {
      // Check observation covariance matrix dimensions equals the number of satellites
      if (i < (INT)pstEpochData_->cNumSV)
      {
         // Set satellite PRN number
         stRelInfo.cPRN = pstEpochData_->vstSVInfo.at(i).cPRN;
         stRelInfo.dAzimuth = pstEpochData_->vstSVInfo.at(i).dAzimuth;
         stRelInfo.dElevation = pstEpochData_->vstSVInfo.at(i).dElevation;
      }
      else
      {
         stRelInfo.cPRN = 33;
         stRelInfo.dAzimuth = 0;
         stRelInfo.dElevation = 0;
      }
      // Set Internal reliability
      stRelInfo.dRelIntern = dDelta * pclObsCovMatrixR_->GetComponent(i, i) / sqrt(pclResCovMatrixCr_->GetComponent(i, i));
      // Set pseudo-range residual
      stRelInfo.dResidual = pclDeltaPseudoRange_->GetComponent(i, 0);
      // Set vector to minimum detectable blunder
      clExternVector.SetComponent(i, 0, stRelInfo.dRelIntern);
      // Compute ECEF external reliability and then rotate to LLF
      if (bECEFCoord_)
      {
         clECEFExtern = clTempMatrix * clExternVector;
         clLLFExtern = clRle * clECEFExtern;
         // Clock external reliability
         stRelInfo.dRelExternClk = clECEFExtern.GetComponent(3, 0);
      }
      else
      {
         clLLFExtern = clTempMatrix * clExternVector;
         clECEFExtern = clRle.Transpose() * clLLFExtern;
         // Clock external reliability
         stRelInfo.dRelExternClk = clLLFExtern.GetComponent(3, 0);
      }
      // Local level frame external reliability
      stRelInfo.dRelExternE = clLLFExtern.GetComponent(0, 0);
      stRelInfo.dRelExternN = clLLFExtern.GetComponent(1, 0);
      stRelInfo.dRelExternU = clLLFExtern.GetComponent(2, 0);
      // ECEF external reliability
      stRelInfo.dRelExternX = clECEFExtern.GetComponent(0, 0);
      stRelInfo.dRelExternY = clECEFExtern.GetComponent(1, 0);
      stRelInfo.dRelExternZ = clECEFExtern.GetComponent(2, 0);
      // Set vector back to zero
      clExternVector.SetComponent(i, 0, 0);
      pstRelTest_->vReliability.push_back(stRelInfo);      
   }
}