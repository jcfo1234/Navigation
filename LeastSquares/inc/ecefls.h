//----------------------------------------------------------------------
// ecefls.h
// Module that computes position in ECEF coordinates
//----------------------------------------------------------------------

#ifndef ECEFLS_H
#define ECEFLS_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <map>
#include <math.h>
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "Common.h"
using namespace std;

#define ECEF_TOLERANCE_FACTOR ((DOUBLE) 1000)

//------------------------------------------------------------------------------------------------
// Compute design observation matrix
// delta_measurements = H * delta_states
// pstEpochData_ (In): Pointer to structure with satellite geometry per epoch and measurements
// pclPointExpansionEcef_ (In): Point of expansion to estimate satellite geometry
//------------------------------------------------------------------------------------------------
CMatrix ComputeDesignECEFMatrixH(EpochInfo* pstEpochData_, CMatrix* pclPointExpansionEcef_);

//------------------------------------------------------------------------------------------------
// Compute ECEF position
// pstEpochData_ (In): Pointer to structure with satellite geometry per epoch and measurements
// pclProbDist_ (In): Class with Normal and Chi-squared probability distributions
// pclEcefCoordinates_ (In): Point of expansion to estimate receiver position
// dStd_ (In): Constant zenith standard deviation
// stKey_ (In): Select from constant or elevation observation covariance matrix
// dAlpha_ (In): Significance level
// dBeta_ (In): Power of the test
// pvRelTest_ (Out): Internal/External reliability test vector
//------------------------------------------------------------------------------------------------
CMatrix ComputeECEFPosition(EpochInfo* pstEpochData_, Distributions* pclProbDist_, 
                            CMatrix* pclEcefCoordinates_, DOUBLE dStd_, string stKey_, 
                            DOUBLE dAlpha_, DOUBLE dBeta_, vector<ReliablityTest>* pvRelTest_, BOOLEANO bBlunderEnable_);

#endif