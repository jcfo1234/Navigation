//----------------------------------------------------------------------
// Kalman.h
// Module that computes position with Kalman filter
//----------------------------------------------------------------------

#ifndef KALMAN_H
#define KALMAN_H

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <tchar.h>
#include <algorithm>
#include <map>
#include <math.h>
#include "..\..\FileCommon\inc\FileCommon.h"
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "..\..\Distributions\inc\Distributions.h"
#include "PRW.h"
#include "VRW.h"

//----------------------------------------------------------------------------
// Kalman solution: Compute Kalman solutions
// bTheVRW_: Enable processing using Velocity Random Walk model
// bThePRW_: Enable processing using Position Random Walk model
// bUsePseudorate_: Enable using Doppler measurements
// pdTheInitLatitude: Pointer to the initial latitude
// pdTheInitLongitude: Pointer to the initial longitude
// pdTheInitHeight: Pointer to the initial height
// pdTheInitVelEast: Pointer to the initial velocity east
// pdTheInitVelNorth: Pointer to the initial velocity north
// pdTheInitVelUp: Pointer to the initial velocity vertical
// pdTheInitClkBias: Pointer to the initial clock bias
// pdTheInitClkDrift: Pointer to the initial clock drift
// pdTheTruthLatitude_: Pointer to the truth latitude (static only)
// pdTheTruthLongitude_: Pointer to the truth longitude (static only)
// pdTheTruthHeight_: Pointer to the truth height (static only)
// psTheTruthFilePath_: Pointer to the truth trajectory file
// psTheFullFilePath_: Pointer to the data file
// pdThePSRStd_: Pointer to the measurement pseudorange standard deviation
// pdThePSRRateStd_: Pointer to the measurement Doppler standard deviation
// pdThePRWSpectralQ_: Position random walk states (East, North, Up) spectral noise density
// pdTheVRWSpectralQ_: Velocity random walk states (VEast, VNorth, VUp) spectral noise density
// pdTheClkBiasQ_: Clock bias state spectral noise density
// pdTheClkDriftQ_: Clock drift state spectral noise density
// pdTheInitPosStd_: Initial position states (East, North, Up) standard deviation
// pdTheInitVelStd_: Initial velocity states (VEast, VNorth, VUp) standard deviation
//----------------------------------------------------------------------------
void KalmanSolution(BOOLEANO bTheVRW_, BOOLEANO bThePRW_, BOOLEANO bUsePseudorate_,
                    DOUBLE* pdTheInitLatitude, DOUBLE* pdTheInitLongitude, DOUBLE* pdTheInitHeight,
                    DOUBLE* pdTheInitVelEast, DOUBLE* pdTheInitVelNorth, DOUBLE* pdTheInitVelUp,
                    DOUBLE* pdTheInitClkBias, DOUBLE* pdTheInitClkDrift, DOUBLE* pdTheTruthLatitude_,
                    DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_, string* psTheTruthFilePath_,
                    string* psTheFullFilePath_, DOUBLE* pdThePSRStd_, DOUBLE* pdThePSRRateStd_,
                    DOUBLE* pdThePRWSpectralQ_, DOUBLE* pdTheVRWSpectralQ_, DOUBLE* pdTheClkBiasQ_,
                    DOUBLE* pdTheClkDriftQ_, DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_);

#endif