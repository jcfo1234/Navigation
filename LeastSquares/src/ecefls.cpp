//----------------------------------------------------------------------
// ecefls.cpp
// Source file of ecefls.h
// Source code for ecef least squares computations
//----------------------------------------------------------------------
#include "..\inc\ecefls.h"

typedef CMatrix (*Rfunc)(EpochInfo*, DOUBLE);

//------------------------------------------------------------------------------------------------
CMatrix ComputeDesignECEFMatrixH(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionEcef_)
{
   // Verify the point of expansion is a 4x1 vector
   // <X,Y,Z,c*deltat>
   if (pclPointExpansionEcef_->GetNumCols() != 1 && pclPointExpansionEcef_->GetNumRows() != 4)
      throw 6;

   CMatrix clGeometryH("Design Matrix H", (INT)pstEpochData_->cNumSV, 4);
   DOUBLE dRange = 0;

   // Compute each row of the H matrix
   for (INT i = 0; i < clGeometryH.GetNumRows(); i++)
   {
      dRange = sqrt(pow(pstEpochData_->vstSVInfo.at(i).dXPosition - pclPointExpansionEcef_->GetComponent(0, 0), 2) +
                    pow(pstEpochData_->vstSVInfo.at(i).dYPosition - pclPointExpansionEcef_->GetComponent(1, 0), 2) +
                    pow(pstEpochData_->vstSVInfo.at(i).dZPosition - pclPointExpansionEcef_->GetComponent(2, 0), 2));
      // Set row components of H matrix
      clGeometryH.SetComponent(i, 0, (pclPointExpansionEcef_->GetComponent(0, 0) - pstEpochData_->vstSVInfo.at(i).dXPosition) / dRange);
      clGeometryH.SetComponent(i, 1, (pclPointExpansionEcef_->GetComponent(1, 0) - pstEpochData_->vstSVInfo.at(i).dYPosition) / dRange);
      clGeometryH.SetComponent(i, 2, (pclPointExpansionEcef_->GetComponent(2, 0) - pstEpochData_->vstSVInfo.at(i).dZPosition) / dRange);
      clGeometryH.SetComponent(i, 3, 1);
   }

   return clGeometryH;
}

//------------------------------------------------------------------------------------------------
CMatrix ComputeECEFPosition(EpochInfo* pstEpochData_, Distributions* pclProbDist_, 
                            CMatrix* pclEcefCoordinates_, DOUBLE dStd_, string stKey_,
                            DOUBLE dAlpha_, DOUBLE dBeta_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bBlunderEnable_)
{   
   map<string, Rfunc> funcMap;
   funcMap["CONSTANT"] = &ConstObsCovariancMatrixR;
   funcMap["ELEVATION"] = &ElevObsCovariancMatrixR;

   // Find the function to compute the observation covariance matrix
   if (funcMap.find(stKey_) == funcMap.end())
      throw 10;

   // Perform blunder detection
   BOOLEANO bBlunderPass = FALSE;
   INT iNumBlunders = 1;
   CMatrix clEcefCoordinates(*pclEcefCoordinates_);
   // Compute pseudo-range residuals
   CMatrix clDeltaPseudo(CalculatePseudorangeResiduals(pstEpochData_, pclEcefCoordinates_));
   CMatrix clDesignH(ComputeDesignECEFMatrixH(pstEpochData_, pclEcefCoordinates_));
   CMatrix clPseudoCovarianceR(funcMap.at(stKey_)(pstEpochData_, dStd_));
   CMatrix clErrorCovarianceP(ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH));
   CMatrix clResidualCovarianceCr(ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR));
   CMatrix clErrorEcef(ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo));

   DOUBLE dErrorNorm2 = sqrt(pow(clErrorEcef.GetComponent(0, 0), 2) +
                             pow(clErrorEcef.GetComponent(1, 0), 2) + 
                             pow(clErrorEcef.GetComponent(2, 0), 2) +
                             pow(clErrorEcef.GetComponent(3, 0), 2));

   // Re-compute coordinates until blunder detectrion test passes
   while (!bBlunderPass && pstEpochData_->vstSVInfo.size() > 5)
   {
      // Iterate until solution is satisfactory
      while (dErrorNorm2 > clErrorCovarianceP.MatrixNorm2() / ECEF_TOLERANCE_FACTOR)
      {
         // Update point of expansion
         *pclEcefCoordinates_ = *pclEcefCoordinates_ + clErrorEcef;
         // Least squares computations
         clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, pclEcefCoordinates_);
         clDesignH = ComputeDesignECEFMatrixH(pstEpochData_, pclEcefCoordinates_);
         clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
         clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH);
         clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR);
         clErrorEcef = ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo);
         dErrorNorm2 = sqrt(pow(clErrorEcef.GetComponent(0, 0), 2) +
                            pow(clErrorEcef.GetComponent(1, 0), 2) +
                            pow(clErrorEcef.GetComponent(2, 0), 2) +
                            pow(clErrorEcef.GetComponent(3, 0), 2));
      }
      // Do Blunder detection
      if (bBlunderEnable_)
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
            ReliabilityTest(pstEpochData_, &clPseudoCovarianceR, &clResidualCovarianceCr, &clDesignH, pclEcefCoordinates_, TRUE, &clErrorCovarianceP, &clDeltaPseudo, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
            pvRelTest_->push_back(stRelTest);
         }
         bBlunderPass = BlunderDetection(pstEpochData_, pclProbDist_, &clDeltaPseudo, &clResidualCovarianceCr, dAlpha_, iNumBlunders);
         // Blunder detection test fails
         if (!bBlunderPass)
         {
            *pclEcefCoordinates_ = clEcefCoordinates;
            clDeltaPseudo = CalculatePseudorangeResiduals(pstEpochData_, pclEcefCoordinates_);
            clDesignH = ComputeDesignECEFMatrixH(pstEpochData_, pclEcefCoordinates_);
            clPseudoCovarianceR = funcMap.at(stKey_)(pstEpochData_, dStd_);
            clErrorCovarianceP = ErrorStateCovarianceP(&clPseudoCovarianceR, &clDesignH);
            clResidualCovarianceCr = ComputeResidualCovarianceMatrix(&clDesignH, &clPseudoCovarianceR);
            clErrorEcef = ErrorCorrectionVector(&clPseudoCovarianceR, &clDesignH, &clDeltaPseudo);
            dErrorNorm2 = sqrt(pow(clErrorEcef.GetComponent(0, 0), 2) +
               pow(clErrorEcef.GetComponent(1, 0), 2) +
               pow(clErrorEcef.GetComponent(2, 0), 2) +
               pow(clErrorEcef.GetComponent(3, 0), 2));
         }
         iNumBlunders = iNumBlunders + 1;
      }
      // Blunder detection disabled
      else
      {
         ReliablityTest stRelTest;
         GlobalCovarianceTest(&clDeltaPseudo, &clPseudoCovarianceR, &stRelTest, pclProbDist_, dAlpha_);
         ReliabilityTest(pstEpochData_, &clPseudoCovarianceR, &clResidualCovarianceCr, &clDesignH, pclEcefCoordinates_, TRUE, &clErrorCovarianceP, &clDeltaPseudo, pclProbDist_, dAlpha_, dBeta_, &stRelTest);
         pvRelTest_->push_back(stRelTest);
         bBlunderPass = TRUE;
      }
   }

   return clErrorCovarianceP;
}
