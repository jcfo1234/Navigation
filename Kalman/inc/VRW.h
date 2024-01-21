//----------------------------------------------------------------------
// VRW.h
// Velocity Random Walk dynamic model
//----------------------------------------------------------------------

#ifndef VRW_H
#define VRW_H

#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <algorithm>
#include <map>
#include <math.h>
#include "..\..\DataTypeDefinitions.h"
#include "..\..\Matrix\inc\Matrix.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"
#include "Kmeasurements.h"
using namespace std;

#define VRW_STATE_DIMENSION ((ULONG) 8)
#define VRW_NUM_NOISE_SOURCES ((ULONG) 5)

//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix VRWStateSpaceMatrixF();

//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix VRWNoiseShapeMatrixG();

//-------------------------------------------------------------------
// VRWInitializeState
// Initializes the Velocity Random Walk model system state
//-------------------------------------------------------------------
CMatrix VRWInitializeState(DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_,
                           DOUBLE* pdTheHeight_, DOUBLE* pdTheVelEast_,
                           DOUBLE* pdTheVelNorth_, DOUBLE* pdTheVelUp_, 
                           DOUBLE* pdTheClkBias_, DOUBLE* pdTheClkDrift_);

//-------------------------------------------------------------------
// VRWInitializeStateCovarianceP
// Initialize the initial state covaraince matrix
//-------------------------------------------------------------------
CMatrix VRWInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_, DOUBLE* pdTheInitVelStd_);

//-------------------------------------------------------------------
// VRWUpdateLoop
// pclSysState_: Apriori and then aposteriori system state
//               (latitude, longitude, height, clock bias, clock drift)
//               latitude and longitude must be in radians
// pclSysCovP_: Apriori and then aposteriori state covariance matrix
// pstEpochInfo_: Satellite geometry per epoch
// dStdPSR_: Pseudorange standard deviation at zenith
// dStdPSRRate_: Pseudorange rate standard deviation at zenith
// bUsePseudorate_: Enable using pseudorate measurements
//-------------------------------------------------------------------
void VRWUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_,
                   EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_, DOUBLE dStdPSRRate_, BOOLEANO bUsePseudorate_=FALSE);

//-------------------------------------------------------------------
// VRWFillPositionError
// pstTruthTraj_: Truth trajectory
// pstKalmanTraj_: Kalman estimated trajectory
// pclStateCovP_: Aposteriori system state covariance matrix
// pstPRWKalmanError_: Trajectory error and correlation coefficients
// cNumSV_: Number of satellites per epoch
//-------------------------------------------------------------------
void VRWFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_,
                          CMatrix* pclStateCovP_, vector<VRWErrorInfo>* pvVRWKalmanError_, CHAR cNumSV_);

#endif
