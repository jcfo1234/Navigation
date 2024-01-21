//----------------------------------------------------------------------
// geoidls.cpp
// Source file of geoidls.h
// Source code for curvilinear least squares computations
//----------------------------------------------------------------------
#include "..\inc\geoidls.h"

typedef CMatrix(*RfuncGeoid)(EpochInfo*, DOUBLE);

//------------------------------------------------------------------------------------------------
CMatrix ComputeDesignLLFMatrixH(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionGeoid_)
{
   // Verify the point of expansion is a 4x1 vector
   // <latitude,longitude,height,c*deltat>
   if (pclPointExpansionGeoid_->GetNumCols() != 1 && pclPointExpansionGeoid_->GetNumRows() != 4)
      throw 7;

   // Hi = CTS
   CMatrix clAzimuthElevationC("Azimuth Elevation C", 1, 3);
   CMatrix clPseudorangeGeoidT("T Matrix", 3, 3);
   CMatrix clGeoidToLlfS(ComputeGeoidToLlfS(pclPointExpansionGeoid_));
   CMatrix clDesignMatrixH("Design Matrix", (INT)pstEpochData_->cNumSV, 4);
   CMatrix clRowH("Row of Design Matrix", 1, 3);

   // Build the design matrix H
   for (INT i = 0; i < clDesignMatrixH.GetNumRows(); i++)
   {
      clAzimuthElevationC = ComputeAzimuthElevC(&pstEpochData_->vstSVInfo.at(i));
      clPseudorangeGeoidT = ComputePseudorangeGeoidT(&pstEpochData_->vstSVInfo.at(i), pclPointExpansionGeoid_);
      clRowH = clAzimuthElevationC * clPseudorangeGeoidT * clGeoidToLlfS;
      // Set components of design matrix H
      clDesignMatrixH.SetComponent(i, 0, clRowH.GetComponent(0, 0));
      clDesignMatrixH.SetComponent(i, 1, clRowH.GetComponent(0, 1));
      clDesignMatrixH.SetComponent(i, 2, clRowH.GetComponent(0, 2));
      clDesignMatrixH.SetComponent(i, 3, 1);
   }

   return clDesignMatrixH;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeAzimuthElevC(SVInfo* pstSVData_)
{
   CMatrix C("Azimuth Elevation C", 1, 3);
   C.SetComponent(0, 0, cos(DegToRad(pstSVData_->dElevation)) * sin(DegToRad(pstSVData_->dAzimuth)));
   C.SetComponent(0, 1, cos(DegToRad(pstSVData_->dElevation)) * cos(DegToRad(pstSVData_->dAzimuth)));
   C.SetComponent(0, 2, sin(DegToRad(pstSVData_->dElevation)));

   return C;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputePseudorangeGeoidT(SVInfo* pstSVData_, CMatrix* pclPointExpansionGeoid_)
{
   // Verify the point of expansion is a 4x1 vector
   // <latitude,longitude,height,c*deltat>
   if (pclPointExpansionGeoid_->GetNumCols() != 1 && pclPointExpansionGeoid_->GetNumRows() != 4)
      throw 8;

   // Extract latitude from point of expansion to compute N and M
   CMatrix dLatSingleton("Latitude", 1, 1);
   dLatSingleton.SetComponent(0, 0, pclPointExpansionGeoid_->GetComponent(0, 0));
   CMatrix NSingleton = CalcN(&dLatSingleton);
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
   T.SetComponent(0, 1, -pstSVData_->dXPosition * cos(pclPointExpansionGeoid_->GetComponent(1, 0)) - pstSVData_->dYPosition * sin(pclPointExpansionGeoid_->GetComponent(1, 0)));
   // T(0,2)
   T.SetComponent(0, 2, 0);
   // T(1,0)
   T.SetComponent(1, 0, -pstSVData_->dXPosition * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dYPosition * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) * sin(pclPointExpansionGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dZPosition * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) +
                         pow(dEcc, 2) * dNdLatitude * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) +
                         pow(dEcc, 2) * NSingleton.GetComponent(0, 0) * (pow(cos(pclPointExpansionGeoid_->GetComponent(0, 0)), 2) - pow(sin(pclPointExpansionGeoid_->GetComponent(0, 0)), 2)));
   // T(1,1)
   T.SetComponent(1, 1, pstSVData_->dXPosition * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * sin(pclPointExpansionGeoid_->GetComponent(1, 0)) -
                        pstSVData_->dYPosition * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(1, 0)));
   // T(1,2)
   T.SetComponent(1, 2, 0);
   // T(2,0)
   T.SetComponent(2, 0, -pstSVData_->dXPosition * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(1, 0)) -
                         pstSVData_->dYPosition * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * sin(pclPointExpansionGeoid_->GetComponent(1, 0)) +
                         pstSVData_->dZPosition * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) -
                         dNdLatitude * (1 - pow(dEcc, 2) * pow(sin(pclPointExpansionGeoid_->GetComponent(0, 0)), 2)) +
                         2 * pow(dEcc, 2) * NSingleton.GetComponent(0, 0) * sin(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(0, 0)));
   // T(2,1)
   T.SetComponent(2, 1, -pstSVData_->dXPosition * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) * sin(pclPointExpansionGeoid_->GetComponent(1, 0)) +
                         pstSVData_->dYPosition * cos(pclPointExpansionGeoid_->GetComponent(0, 0)) * cos(pclPointExpansionGeoid_->GetComponent(1, 0)));
   // T(2,2)
   T.SetComponent(2, 2, -1);

   return T;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidToLlfS(CMatrix* pclPointExpansionGeoid_)
{
   // Verify the point of expansion is a 4x1 vector
   // <latitude,longitude,height,c*deltat>
   if (pclPointExpansionGeoid_->GetNumCols() != 1 && pclPointExpansionGeoid_->GetNumRows() != 4)
      throw 9;

   // Extract latitude from point of expansion to compute N and M
   CMatrix dLatSingleton("Latitude", 1, 1);
   dLatSingleton.SetComponent(0, 0, pclPointExpansionGeoid_->GetComponent(0, 0));
   DOUBLE dN = CalcN(&dLatSingleton).GetComponent(0, 0);
   DOUBLE dM = CalcM(&dLatSingleton).GetComponent(0, 0);
   CMatrix S("S Matrix", 3, 3);
   // Set components of matrix
   S.SetComponent(0, 0, 0);
   S.SetComponent(0, 1, 1 / (dM + pclPointExpansionGeoid_->GetComponent(2, 0)));
   S.SetComponent(0, 2, 0);
   S.SetComponent(1, 0, 1 / ((dN + pclPointExpansionGeoid_->GetComponent(2, 0)) * cos(pclPointExpansionGeoid_->GetComponent(0, 0))));
   S.SetComponent(1, 1, 0);
   S.SetComponent(2, 0, 0);
   S.SetComponent(2, 1, 0);
   S.SetComponent(2, 2, 1);

   return S;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void UpdatePointExpansion(CMatrix* pclGeoidCoordinates_, CMatrix* pclLLFCorrections_)
{
   CMatrix clLatSingleton("Latitude", 1, 1);
   clLatSingleton.SetComponent(0, 0, pclGeoidCoordinates_->GetComponent(0, 0));
   CMatrix clMSingleton(CalcM(&clLatSingleton));
   CMatrix clNSingleton(CalcN(&clLatSingleton));
   DOUBLE dLatCorrection = pclLLFCorrections_->GetComponent(1, 0) / (clMSingleton.GetComponent(0, 0) + pclGeoidCoordinates_->GetComponent(2, 0));
   DOUBLE dLngCorrection = pclLLFCorrections_->GetComponent(0, 0) / ((clNSingleton.GetComponent(0, 0) + pclGeoidCoordinates_->GetComponent(2, 0)) * cos(pclGeoidCoordinates_->GetComponent(0, 0)));
   // Updated vector
   pclGeoidCoordinates_->SetComponent(0, 0, pclGeoidCoordinates_->GetComponent(0, 0) + dLatCorrection);
   pclGeoidCoordinates_->SetComponent(1, 0, pclGeoidCoordinates_->GetComponent(1, 0) + dLngCorrection);
   pclGeoidCoordinates_->SetComponent(2, 0, pclGeoidCoordinates_->GetComponent(2, 0) + pclLLFCorrections_->GetComponent(2, 0));
   pclGeoidCoordinates_->SetComponent(3, 0, pclGeoidCoordinates_->GetComponent(3, 0) + pclLLFCorrections_->GetComponent(3, 0));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidPosition(EpochInfo* pstEpochData_, Distributions* pclProbDist_,
                             CMatrix* pclGeoidCoordinates_, DOUBLE dStd_, string stKey_,
                             DOUBLE dAlpha_, DOUBLE dBeta_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bEnableBlunder_)
{
   map<string, RfuncGeoid> funcMap;
   funcMap["CONSTANT"] = &ConstObsCovariancMatrixR;
   funcMap["ELEVATION"] = &ElevObsCovariancMatrixR;
   // Find the function to compute the observation covariance matrix
   if (funcMap.find(stKey_) == funcMap.end())
      throw 11;

   // Perform blunder detection
   BOOLEANO bBlunderPass = FALSE;
   INT iNumBlunders = 1;
   CMatrix clGeoidCoordinates(*pclGeoidCoordinates_);
   // Compute ECEF coordinates for pseudo-range residuals
   CMatrix clEcefCoordinates_("ECEF", 4, 1);
   GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
   // Compute pseudo-range residuals
   CMatrix clDeltaPseudo(CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_));
   CMatrix clDesignH(ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_));
   CMatrix clPseudoCovarianceR(funcMap.at(stKey_)(pstEpochData_, dStd_));
   CMatrix clErrorCovarianceP(ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH));
   CMatrix clResidualCovarianceCr(ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR));
   CMatrix clErrorLLF(ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo));
   // Norm of error state vector
   DOUBLE dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
                             pow(clErrorLLF.GetComponent(1, 0), 2) + 
                             pow(clErrorLLF.GetComponent(2, 0), 2) +
                             pow(clErrorLLF.GetComponent(3, 0), 2));

   // Re-compute coordinates until blunder detectrion test passes
   while (!bBlunderPass)
   {
      // Iterate until solution is satisfactory
      while (dErrorNorm2 > clErrorCovarianceP.MatrixNorm2() / LLF_TOLERANCE_FACTOR)
      {
         // Update point of expansion
         UpdatePointExpansion(pclGeoidCoordinates_, &clErrorLLF);
         GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
         clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_);
         clDesignH = ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_);
         clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
         clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH);
         clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR);
         clErrorLLF = ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo);
         dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
                            pow(clErrorLLF.GetComponent(1, 0), 2) +
                            pow(clErrorLLF.GetComponent(2, 0), 2) +
                            pow(clErrorLLF.GetComponent(3, 0), 2));
      }
      // Blunder detection enabled
      if (bEnableBlunder_)
      {
         // Running out of measurements
         if (pstEpochData_->vstSVInfo.size() <= 5)
         {
            break;
         }
         // Do reliability testing before blunder detection
         if (iNumBlunders == 1)
         {
            ReliablityTest stRelTest;
            GlobalCovarianceTest(&clDeltaPseudo, &clPseudoCovarianceR, &stRelTest, pclProbDist_, dAlpha_);
            ReliabilityTest(pstEpochData_, &clPseudoCovarianceR, &clResidualCovarianceCr, &clDesignH, pclGeoidCoordinates_, FALSE, &clErrorCovarianceP, &clDeltaPseudo, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
            pvRelTest_->push_back(stRelTest);
         }
         bBlunderPass = BlunderDetection(pstEpochData_, pclProbDist_, &clDeltaPseudo, &clResidualCovarianceCr, dAlpha_, iNumBlunders);
         // Update vectors if blunder detection test failed
         if (!bBlunderPass)
         {
            *pclGeoidCoordinates_ = clGeoidCoordinates;
            GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
            clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_);
            clDesignH = ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_);
            clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
            // Compute error state corrections and covarainces
            clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH);
            clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR);
            clErrorLLF = ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo);
            dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
               pow(clErrorLLF.GetComponent(1, 0), 2) +
               pow(clErrorLLF.GetComponent(2, 0), 2) +
               pow(clErrorLLF.GetComponent(3, 0), 2));
         }
         iNumBlunders = iNumBlunders + 1;
      }
      // Blunder detection disabled
      else
      {
         ReliablityTest stRelTest;
         GlobalCovarianceTest(&clDeltaPseudo, &clPseudoCovarianceR, &stRelTest, pclProbDist_, dAlpha_);
         ReliabilityTest(pstEpochData_, &clPseudoCovarianceR, &clResidualCovarianceCr, &clDesignH, pclGeoidCoordinates_, FALSE, &clErrorCovarianceP, &clDeltaPseudo, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
         pvRelTest_->push_back(stRelTest);
         bBlunderPass = TRUE;
      }
   }

   return clErrorCovarianceP;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidPositionConstrained(EpochInfo* pstEpochData_, Distributions* pclProbDist_,
                                        CMatrix* pclGeoidCoordinates_, DOUBLE dStd_, string stKey_,
                                        DOUBLE dHgtConstrained_, DOUBLE dAlpha_, DOUBLE dBeta_,
                                        DOUBLE dHgtConstrainedStd_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bEnableBlunder_)
{
   map<string, RfuncGeoid> funcMap;
   funcMap["CONSTANT"] = &ConstObsCovariancMatrixR;
   funcMap["ELEVATION"] = &ElevObsCovariancMatrixR;
   // Find the function to compute the observation covariance matrix
   if (funcMap.find(stKey_) == funcMap.end())
      throw 11;

   // Perform blunder detection
   BOOLEANO bBlunderPass = FALSE;
   INT iNumBlunders = 1;
   CHAR cNumSV = pstEpochData_->cNumSV;
   CMatrix clGeoidCoordinates(*pclGeoidCoordinates_);
   BOOLEANO bHgtConstrainedBlunder = FALSE;
   // Compute ECEF coordinates for pseudo-range residuals
   CMatrix clEcefCoordinates_("ECEF", 4, 1);
   GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
   // Compute pseudo-range residuals
   CMatrix clDeltaPseudo(CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_));
   CMatrix clDesignH(ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_));
   CMatrix clPseudoCovarianceR(funcMap.at(stKey_)(pstEpochData_, dStd_));
   
   // Augment dimensions of matrices with the extra height constraint
   CMatrix clDeltaPseudoAug(clDeltaPseudo.AugmentDim(clDeltaPseudo, 1, 0));
   CMatrix clDesignHAug(clDesignH.AugmentDim(clDesignH, 1, 0));
   CMatrix clPseudoCovarianceRAug(clPseudoCovarianceR.AugmentDim(clPseudoCovarianceR, 1, 1));
   // Set constrained measurements differentials and covariances
   clDeltaPseudoAug.SetComponent((INT)pstEpochData_->cNumSV, 0, dHgtConstrained_ - pclGeoidCoordinates_->GetComponent(2, 0));
   clDesignHAug.SetComponent((INT)pstEpochData_->cNumSV, 2, 1);
   clPseudoCovarianceRAug.SetComponent((INT)pstEpochData_->cNumSV, (INT)pstEpochData_->cNumSV, pow(dHgtConstrainedStd_, 2));

   // Compute error state corrections and covariances
   CMatrix clErrorCovarianceP(ErrorStateCovarianceP(&clPseudoCovarianceRAug, &clDesignHAug));
   CMatrix clResidualCovarianceCr(ComputeResidualCovarianceMatrix(&clDesignHAug, &clPseudoCovarianceRAug));
   CMatrix clErrorLLF(ErrorCorrectionVector(&clPseudoCovarianceRAug, &clDesignHAug, &clDeltaPseudoAug));

   // Norm of error state vector
   DOUBLE dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
                             pow(clErrorLLF.GetComponent(1, 0), 2) +
                             pow(clErrorLLF.GetComponent(2, 0), 2) +
                             pow(clErrorLLF.GetComponent(3, 0), 2));

   // Re-compute coordinates until blunder detectrion test passes
   while (!bBlunderPass)
   {
      // Iterate until solution is satisfactory
      while (dErrorNorm2 > clErrorCovarianceP.MatrixNorm2() / LLF_TOLERANCE_FACTOR)
      {
         // Update point of expansion
         UpdatePointExpansion(pclGeoidCoordinates_, &clErrorLLF);
         GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
         // Pseudorange residuals
         clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_);
         clDesignH = ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_);
         clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
         // Height constrained observation was not the blunder
         if (!bHgtConstrainedBlunder)
         {
            // Augment vector with height constraint
            clDeltaPseudoAug = clDeltaPseudo.AugmentDim(clDeltaPseudo, 1, 0);
            clDesignHAug = clDesignH.AugmentDim(clDesignH, 1, 0);
            clPseudoCovarianceRAug = clPseudoCovarianceR.AugmentDim(clPseudoCovarianceR, 1, 1);
            // Set constrained measurements differentials and covariances
            clDeltaPseudoAug.SetComponent((INT)pstEpochData_->cNumSV, 0, dHgtConstrained_ - pclGeoidCoordinates_->GetComponent(2, 0));
            clDesignHAug.SetComponent((INT)pstEpochData_->cNumSV, 2, 1);
            clPseudoCovarianceRAug.SetComponent((INT)pstEpochData_->cNumSV, (INT)pstEpochData_->cNumSV, pow(dHgtConstrainedStd_, 2));
         }
         // Heigth constrained observation was the blunder
         else
         {
            clDeltaPseudoAug = clDeltaPseudo;
            clDesignHAug = clDesignH;
            clPseudoCovarianceRAug = clPseudoCovarianceR;
         }
         // Compute error state corrections and covarainces
         clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceRAug, &clDesignHAug);
         clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignHAug, &clPseudoCovarianceRAug);
         clErrorLLF = ErrorCorrectionVector(&clPseudoCovarianceRAug, &clDesignHAug, &clDeltaPseudoAug);
         dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
                            pow(clErrorLLF.GetComponent(1, 0), 2) +
                            pow(clErrorLLF.GetComponent(2, 0), 2) +
                            pow(clErrorLLF.GetComponent(3, 0), 2));
      }
      if (bEnableBlunder_)
      {
         // Running out of measurements
         if (pstEpochData_->vstSVInfo.size() <= 5)
         {
            break;
         }
         // Do reliability testing before blunder detection
         if (iNumBlunders == 1)
         {
            ReliablityTest stRelTest;
            GlobalCovarianceTest(&clDeltaPseudoAug, &clPseudoCovarianceRAug, &stRelTest, pclProbDist_, dAlpha_);
            ReliabilityTest(pstEpochData_, &clPseudoCovarianceRAug, &clResidualCovarianceCr, &clDesignHAug, pclGeoidCoordinates_, FALSE, &clErrorCovarianceP, &clDeltaPseudoAug, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
            pvRelTest_->push_back(stRelTest);
         }
         cNumSV = pstEpochData_->cNumSV;
         bBlunderPass = BlunderDetection(pstEpochData_, pclProbDist_, &clDeltaPseudoAug, &clResidualCovarianceCr, dAlpha_, iNumBlunders);
         // Update vectors if blunder detection test failed
         if (!bBlunderPass)
         {
            *pclGeoidCoordinates_ = clGeoidCoordinates;
            GeoidToECEF(pclGeoidCoordinates_, &clEcefCoordinates_);
            clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, &clEcefCoordinates_);
            clDesignH = ComputeDesignLLFMatrixH(pstEpochData_, pclGeoidCoordinates_);
            clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
            // If it was a pseudo-range which failed the blunder detection test
            if (pstEpochData_->cNumSV < cNumSV)
            {
               clDeltaPseudoAug = clDeltaPseudo.AugmentDim(clDeltaPseudo, 1, 0);
               clDesignHAug = clDesignH.AugmentDim(clDesignH, 1, 0);
               clPseudoCovarianceRAug = clPseudoCovarianceR.AugmentDim(clPseudoCovarianceR, 1, 1);
               // Set constrained measurements differentials and covariances
               clDeltaPseudoAug.SetComponent((INT)pstEpochData_->cNumSV, 0, dHgtConstrained_ - pclGeoidCoordinates_->GetComponent(2, 0));
               clDesignHAug.SetComponent((INT)pstEpochData_->cNumSV, 2, 1);
               clPseudoCovarianceRAug.SetComponent((INT)pstEpochData_->cNumSV, (INT)pstEpochData_->cNumSV, pow(dHgtConstrainedStd_, 2));
            }
            // Blunder test failure was the height observations
            else
            {
               clDeltaPseudoAug = clDeltaPseudo;
               clDesignHAug = clDesignH;
               clPseudoCovarianceRAug = clPseudoCovarianceR;
               bHgtConstrainedBlunder = TRUE;
            }
            // Compute error state corrections and covarainces
            clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceRAug, &clDesignHAug);
            clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignHAug, &clPseudoCovarianceRAug);
            clErrorLLF = ErrorCorrectionVector(&clPseudoCovarianceRAug, &clDesignHAug, &clDeltaPseudoAug);
            dErrorNorm2 = sqrt(pow(clErrorLLF.GetComponent(0, 0), 2) +
               pow(clErrorLLF.GetComponent(1, 0), 2) +
               pow(clErrorLLF.GetComponent(2, 0), 2) +
               pow(clErrorLLF.GetComponent(3, 0), 2));
         }
         iNumBlunders = iNumBlunders + 1;
      }
      else
      {
         ReliablityTest stRelTest;
         GlobalCovarianceTest(&clDeltaPseudoAug, &clPseudoCovarianceRAug, &stRelTest, pclProbDist_, dAlpha_);
         ReliabilityTest(pstEpochData_, &clPseudoCovarianceRAug, &clResidualCovarianceCr, &clDesignHAug, pclGeoidCoordinates_, FALSE, &clErrorCovarianceP, &clDeltaPseudoAug, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
         pvRelTest_->push_back(stRelTest);
         bBlunderPass = TRUE;
      }
   }

   return clErrorCovarianceP;
}