//----------------------------------------------------------------------
// KPolar.h
// Circular motion modeling for Kalman Filter
//----------------------------------------------------------------------

#ifndef KPOLAR_H
#define KPOLAR_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <map>
#include <math.h>
#include "..\..\FileCommon\inc\FileCommon.h"
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "Kmeasurements.h"
using namespace std;

#define KPOLAR_STATE_DIMENSION ((ULONG) 8)
#define KPOLAR_NUM_NOISE_SOURCES ((ULONG) 5)


//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix KPolarStateSpaceMatrixF();

//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix KPolarNoiseShapeMatrixG();

//-------------------------------------------------------------------
// KPolarJacobian
// pclPolar_: Radius, polar angle [rad], angular velocity [rad/sec]
//-------------------------------------------------------------------
CMatrix KPolarJacobian(CMatrix* pclPolar_);

//-------------------------------------------------------------------
// KPolarTransformH
// pclPolar_: polar state (radius, polar angle, angular velocity)
//-------------------------------------------------------------------
CMatrix KPolarTransformH(CMatrix* pclPolar_);

//-------------------------------------------------------------------
// KPolarInitializeState
// Initializes the circular motion model system state
//-------------------------------------------------------------------
CMatrix KPolarInitializeState(DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_,
                              DOUBLE* pdTheHeight_, DOUBLE* pdTheVelEast_,
                              DOUBLE* pdTheVelNorth_, DOUBLE* pdTheVelUp_,
                              DOUBLE* pdTheClkBias_, DOUBLE* pdTheClkDrift_,
                              DOUBLE* pdTheRadius_);

//-------------------------------------------------------------------
// KPolarInitializeStateCovarianceP
// Initialize the initial state covaraince matrix
//-------------------------------------------------------------------
CMatrix KPolarInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_, CMatrix* pclPolar_);

//-------------------------------------------------------------------
// KPolarSpectralDensityQ
// dRWSpectralQ_ (In): Spectral densities of polar coordinates
// dClkBiasQ_ (In): Spectral density of the Clock bias noise
// dClkDriftQ_ (In): Spectral density of the Clock drift noise
//-------------------------------------------------------------------
CMatrix KPolarSpectralDensityQ(DOUBLE* pdRWSpectralQ_, DOUBLE* pdClkBiasQ_, DOUBLE* pdClkDriftQ_, CMatrix* pclPolar_);

//-------------------------------------------------------------------
// KPolarPredictionLoop
// pclF_: System state space matrix
// pclG_: Noise shaping matrix
// pclQ_: Spectral density matrix
// pclSysState_: System state
// pclSysCovP_: State covariance matrix
// dStepSize_: Prediction rate
//-------------------------------------------------------------------
void KPolarPredictionLoop(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
                          CMatrix* pclSysState_, CMatrix* pclSysCovP_, DOUBLE dStepSize_);

//-------------------------------------------------------------------
// KPolarUpdateLoop
// pclSysState_: Apriori and then aposteriori system state
//               (latitude, longitude, height, clock bias, clock drift)
//               latitude and longitude must be in radians
// pclSysCovP_: Apriori and then aposteriori state covariance matrix
// pstEpochInfo_: Satellite geometry per epoch
// dStdPSR_: Pseudorange standard deviation at zenith
// dStdPSRRate_: Pseudorange rate standard deviation at zenith
// bUsePseudorate_: Enable using pseudorate measurements
//-------------------------------------------------------------------
void KPolarUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_,
                      EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_, 
                      DOUBLE dStdPSRRate_, BOOLEANO bUsePseudorate_ = TRUE);

//----------------------------------------------------------------------------
// Kalman solution: Compute Kalman solution in polar coordinates
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
void KPolarSolution(DOUBLE* pdTheInitLatitude, DOUBLE* pdTheInitLongitude, DOUBLE* pdTheInitHeight,
                    DOUBLE* pdTheInitVelEast, DOUBLE* pdTheInitVelNorth, DOUBLE* pdTheInitVelUp,
                    DOUBLE* pdTheInitClkBias, DOUBLE* pdTheInitClkDrift, DOUBLE* pdTheTruthLatitude_,
                    DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_, string* psTheTruthFilePath_,
                    string* psTheFullFilePath_, DOUBLE* pdThePSRStd_, DOUBLE* pdThePSRRateStd_,
                    DOUBLE* pdThePRWSpectralQ_, DOUBLE* pdTheVRWSpectralQ_, DOUBLE* pdTheClkBiasQ_,
                    DOUBLE* pdTheClkDriftQ_, DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_,
                    BOOLEANO bTheUseDoppler, DOUBLE* pdTheInitRadius_); 

//----------------------------------------------------------------------------
void KPolarFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_,
                             CMatrix* pclStateCovP_, vector<VRWErrorInfo>* pvVRWKalmanError_, 
                             CHAR cNumSV_, CMatrix* pclPolar_);

#endif