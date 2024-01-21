//----------------------------------------------------------------------
// geoid.h
// Module that computes position in curvilinear coordinates
//----------------------------------------------------------------------

#ifndef GEOIDLS_H
#define GEOIDLS_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <math.h>
#include <map>
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "Common.h"
using namespace std;

#define LLF_TOLERANCE_FACTOR ((DOUBLE) 1000)

//------------------------------------------------------------------------------------------------
// Compute design observation matrix
// delta_measurements = H * delta_states
// pstEpochData_ (In): Pointer to structure with satellite geometry per epoch and measurements
// pclPointExpansionGeoid_ (In): Point of expansion to estimate satellite geometry
//------------------------------------------------------------------------------------------------
CMatrix ComputeDesignLLFMatrixH(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionGeoid_);

//------------------------------------------------------------------------------------------------
// Compute matrix of azimuth and elevation angles
// C
// pstSVData_ (In): Pointer to structure with satellite information
//------------------------------------------------------------------------------------------------
CMatrix ComputeAzimuthElevC(SVInfo* pstSVData_);

//------------------------------------------------------------------------------------------------
// Compute matrix of range errors in Local Level Frame in terms of curvilinear coordinates
// T
// pstSVData_ (In): Pointer to structure with satellite information
// pclPointExpansionGeoid_ (In): Point of expansion to estimate satellite geometry
//------------------------------------------------------------------------------------------------
CMatrix ComputePseudorangeGeoidT(SVInfo* pstSVData_, CMatrix* pclPointExpansionGeoid_);

//------------------------------------------------------------------------------------------------
// Compute transformation matrix from Geoid to Local level frame coordinates
// S
// pclPointExpansionGeoid_ (In): Point of expansion to estimate satellite geometry
//------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidToLlfS(CMatrix* pclPointExpansionGeoid_);

//------------------------------------------------------------------------------------------------
// Compute Geoid position
// pstEpochData_ (Out): Pointer to structure with satellite geometry per epoch and measurements
// pclProbDist_ (In) : Class with normal and Chi-squared probability distributions
// pclGeoidCoordinates_ (In): Point of expansion to estimate receiver position
// dStd_ (In): Zenith elevation standard deviation
// stKey_ (In): Constant or Elevation based observation covaraince matrix selector
// dAlpha_ (In): Significance
// dBeta_ (In): Power of the test.s
//------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidPosition(EpochInfo* pstEpochData_, Distributions* pclProbDist_, 
                             CMatrix* pclGeoidCoordinates_, DOUBLE dStd_, string stKey_,
                             DOUBLE dAlpha_, DOUBLE dBeta_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bEnableBlunder_);

//------------------------------------------------------------------------------------------------
// Compute Geoid position with height constrained
// pstEpochData_ (Out): Pointer to structure with satellite geometry per epoch and measurements
// pclProbDist_ (In) : Class with normal and Chi-squared probability distributions
// pclGeoidCoordinates_ (In): Point of expansion to estimate receiver position
// dStd_ (In): Zenith elevation standard deviation
// stKey_ (In): Constant or Elevation based observation covaraince matrix selector
// dHgtConstrained_ (In): Value of height constrained
// dAlpha_ (In): Significance
// dBeta_ (In): Power of the test.s
//------------------------------------------------------------------------------------------------
CMatrix ComputeGeoidPositionConstrained(EpochInfo* pstEpochData_, Distributions* pclProbDist_,
                                        CMatrix* pclGeoidCoordinates_, DOUBLE dStd_, string stKey_,
                                        DOUBLE dHgtConstrained_, DOUBLE dAlpha_, DOUBLE dBeta_,
                                        DOUBLE dHgtConstrainedStd_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bEnableBlunder_);

//------------------------------------------------------------------------------------------------
// Update point of expansion
// pclGeoidCoordinates_ (In): Point of expansion to estimate receiver position
// pclLLFCorrections_ (In): Local level frame correction vector
//------------------------------------------------------------------------------------------------
void UpdatePointExpansion(CMatrix* pclGeoidCoordinates_, CMatrix* pclLLFCorrections_);
#endif