//----------------------------------------------------------------------
// PRW.h
// Position Random Walk dynamic model
//----------------------------------------------------------------------

#ifndef PRW_H
#define PRW_H

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

#define PRW_STATE_DIMENSION ((ULONG) 5)

//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix PRWStateSpaceMatrixF();

//-------------------------------------------------------------------
// Compute state space matrix
//-------------------------------------------------------------------
CMatrix PRWNoiseShapeMatrixG();

//-------------------------------------------------------------------
// PRWInitializeState
// Initializes the Position Random Walk model system state
//-------------------------------------------------------------------
CMatrix PRWInitializeState(DOUBLE* pdTheLatitude, DOUBLE* pdTheLongitude,
                           DOUBLE* pdTheHeight, DOUBLE* pdTheClkBias_, 
                           DOUBLE* pdTheClkDrift_);

//-------------------------------------------------------------------
// PRWInitializeStateCovarianceP
// Initialize the initial state covaraince matrix
//-------------------------------------------------------------------
CMatrix PRWInitializeStateCovarianceP(DOUBLE* pdTheInitPosStd_);

//-------------------------------------------------------------------
// PRWUpdateLoop
// pclSysState_: Apriori and then aposteriori system state
//               (latitude, longitude, height, clock bias, clock drift)
//               latitude and longitude must be in radians
// pclSysCovP_: Apriori and then aposteriori state covariance matrix
// pstEpochInfo_: Satellite geometry per epoch
// dStdPSR_: Pseudorange standard deviation at zenith
//-------------------------------------------------------------------
void PRWUpdateLoop(CMatrix* pclSysState_, CMatrix* pclSysCovP_, 
                   EpochInfo* pstEpochInfo_, DOUBLE dStdPSR_);

//-------------------------------------------------------------------
// PRWFillPositionError
// pstTruthTraj_: Truth trajectory
// pstKalmanTraj_: Kalman estimated trajectory
// pclStateCovP_: Aposteriori system state covariance matrix
// pstPRWKalmanError_: Trajectory error and correlation coefficients
// cNumSV_: Number of satellites per epoch
//-------------------------------------------------------------------
void PRWFillPositionError(TrajectoryInfo* pstTruthTraj_, TrajectoryInfo* pstKalmanTraj_, 
                          CMatrix* pclStateCovP_, vector<PRWErrorInfo>* pvPRWKalmanError_, CHAR cNumSV_);
#endif
