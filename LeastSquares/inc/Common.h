//----------------------------------------------------------------------
// Common.h
// Module common to ECEF and Curvilinear computations
//----------------------------------------------------------------------

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <math.h>
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "..\..\Distributions\inc\Distributions.h"
using namespace std;

//------------------------------------------------------------------------------------------------
// Calculate pseudo-range residual vector
// deltaP = P - Phat
// pstEpochData_ (In): Pointer to structure with satellite geometry per epoch and measurements
// pclPointExpansionEcef_ (In): Point of expansion to estimate satellite geometry
//------------------------------------------------------------------------------------------------
CMatrix CalculatePseudorangeResiduals(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionEcef_);

//------------------------------------------------------------------------------------------------
// Guess the measurement covariance matrix based on elevation
// R = dStd_^2 * diag(1/sin^2(elev_i))
// pstEpochData_ (In): Pointer to structure with satellite geometry per epoch and measurements
// dStd_ (In): Standard deviation of pseudo-range measurements at zenith
//------------------------------------------------------------------------------------------------
CMatrix ElevObsCovariancMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_);

//------------------------------------------------------------------------------------------------
// Guess the constant measurement covariance matrix
// R = diag(dStd_^2)
// iNumMeasurements (In): Number of measurements per epoch
// dStd_ (In): Standard deviation of pseudo-range measurements at zenith
//------------------------------------------------------------------------------------------------
CMatrix ConstObsCovariancMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_);

//------------------------------------------------------------------------------------------------
// Compute the error state covariance matrix
// P = (H^T*R^-1*H)^-1
// pclObsCovarianceMatrixR_ (In): Pointer to measurements covariance matrix R
// pclDesignMatrixH_ (In): Pointer to design matrix H
//------------------------------------------------------------------------------------------------
CMatrix ErrorStateCovarianceP(CMatrix* pclObsCovarianceMatrixR_, CMatrix* pclDesignMatrixH_);

//------------------------------------------------------------------------------------------------
// Compute the error correction vector
// dx = (H^T*R^(-1)*H)^-1*H^T*R^(-1)*dP
// pclObsCovarianceMatrixR_ (In): Pointer to measurements covariance matrix R
// pclDesignMatrixH_ (In): Pointer to design matrix H
// pclMeasurementResidual_ (In): Pointer to measurement residual vector
//------------------------------------------------------------------------------------------------
CMatrix ErrorCorrectionVector(CMatrix* pclObsCovarianceMatrixR_, CMatrix* pclDesignMatrixH_, CMatrix* pclMeasurementResidual_);

//------------------------------------------------------------------------------------------------
// Assign Geoid Point of expansion
// dLat_ (In): Latitude
// dLng_ (In): Longitude
// dHgt_ (In): Height
// dTShift_ (In): Time shift difference
// pclPointExpansionGeoid_ (Out): Geoid point of expansion vector for least squares
//------------------------------------------------------------------------------------------------
void SetPointOfExpansionGeoid(DOUBLE dLat_, DOUBLE dLng_, DOUBLE dHgt_, DOUBLE dTshift_, CMatrix* pclPointExpansionGeoid_);

//------------------------------------------------------------------------------------------------
// Assign ECEF Point of expansion
// dLat_ (In): Latitude
// dLng_ (In): Longitude
// dHgt_ (In): Height
// dTShift_ (In): Time shift difference
// pclPointExpansionEcef_ (Out): Ecef point of expansion vector for least squares
//------------------------------------------------------------------------------------------------
void SetPointOfExpansionEcef(DOUBLE dLat_, DOUBLE dLng_, DOUBLE dHgt_, DOUBLE dTshift_, CMatrix* pclPointExpansionEcef_);

//------------------------------------------------------------------------------------------------
// Compute Local Level Frame Dilution of Precision
// pclTransformationMatrix_ (In): Transformation matrix to Local Level Frame
// pclErrorCovMatrix_ (In): Error covariance matrix
//------------------------------------------------------------------------------------------------
CMatrix ComputeDOP(CMatrix* pclTransformationMatrix_, CMatrix* pclErrorCovMatrix_);

//------------------------------------------------------------------------------------------------
// Compute the residual covariance matrix
// pclDesignMatrixH_ (In): Design Geometry Matrix
// pclObsCovarianceMatrixR_ (In): Observation covariance matrix
//------------------------------------------------------------------------------------------------
CMatrix ComputeResidualCovarianceMatrix(CMatrix* pclDesignMatrixH_, CMatrix* pclObsCovarianceMatrixR_);

//------------------------------------------------------------------------------------------------
// Compute the residual covariance matrix
// pstEpochData_ (Out): Satellite set per epoch information structure
// pclProbDist_ (In): Class with Normal and Chi-squared probability distributions
// pclDesignMatrixH_ (In): Observation Geometry matrix H
// pclObsCovarianceMatrixR_ (In): Observation covariance matrix
// pclPointExpansionEcef_ (In): ECEF point of Expansion
// dAlpha_ (In): Confidence level
// dBeta_ (In): Power of the test
//------------------------------------------------------------------------------------------------
BOOLEANO BlunderDetection(EpochInfo* pstEpochData_, Distributions* pclProbDist_, CMatrix* pclPseudoRes_, CMatrix* pclResCovMatrixCr_, DOUBLE dAlpha_, INT iNumBlunders_);

//------------------------------------------------------------------------------------------------
// Fill coordinate vector
// pclCoord4D_ (In): Vector with coordinates and time shift
// pclCoord3D_ (Out): Vector with coordinates only
//------------------------------------------------------------------------------------------------
void FillCoordinates(CMatrix* pclCoord4D_, CMatrix* pclCoord3D_);

//------------------------------------------------------------------------------------------------
// Fill error in trajectory
// pvPosError_ (Out): Trajectory error
// pclGeoidCoordinates_ (In): Current 3D curvilinear coordinates
// pstTrajEpoch_ (In): Current 3D curvilinear truth coordinates
// pclErrCovP_ (In): Error state covariance matrix
// pclTransformation_ (In): Coordinate transformation matrix
//------------------------------------------------------------------------------------------------
void FillPositionErrorVector(vector<ErrorInfo>* pvPosError_, CMatrix* pclGeoidCoordinates_, 
                             TrajectoryInfo* pstTrajEpoch_, CMatrix* pclErrCovP_, 
                             CMatrix* pclTransformation_, CHAR cNumSV_, DOUBLE dMeasStd_);

//------------------------------------------------------------------------------------------------
// Compute global covariance test
// pclDeltaPseudoRange_ (In): Measurement residuals
// pclObsCovMatrix_ (In): Observation covariance matrix
// pstRelTest_ (Out): Reliability test structure
// pclProbDist_ (In): Normal and Xi-squared probability distributions
// dAlpha_ (In): Significance level
//------------------------------------------------------------------------------------------------
void GlobalCovarianceTest(CMatrix* pclDeltaPseudoRange_, CMatrix* pclObsCovMatrix_, 
                          ReliablityTest* pstRelTest_, Distributions* pclProbDist_,
                          DOUBLE dAlpha_);

//------------------------------------------------------------------------------------------------
// Compute internal and external reliability
// pstEpochData_ (In): Structure with satellite geometry
// pclObsCovMatrixR_ (In): Observation covariance matrix
// pclResCovMatrixCr_ (In): Residual covaraince matrix
// pclDesignMatrixH_ (In): Design Matrix
// pclCoordinates_ (In): Vector with coordinates
// bECEFCoord_ (In): Flag indicating if coordinates are given in ECEF frame
// pclDeltaPseudoRange_ (In): Pseudo-range residuals
// pclProbDist_ (In): Normal and Xi-squared probability distributions
// dAlpha_ (In): Significance level
// dBeta_ (In): Power of the test
// pstRelTest_ (Out): Vector with internal and external reliability values
//------------------------------------------------------------------------------------------------
void ReliabilityTest(EpochInfo* pstEpochData_, CMatrix* pclObsCovMatrixR_,
                     CMatrix* pclResCovMatrixCr_, CMatrix* pclDesignMatrixH_, 
                     CMatrix* pclCoordinates_, BOOLEANO bECEFCoord_, CMatrix* pclErrorStateCovP_,
                     CMatrix* pclDeltaPseudoRange_, Distributions* pclProbDist_,
                     DOUBLE dAlpha_, DOUBLE dBeta_, ReliablityTest* pstRelTest_);
#endif