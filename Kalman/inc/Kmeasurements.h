//----------------------------------------------------------------------
// KCommon.h
// Module common to Position and Velocity Random Walk models
//----------------------------------------------------------------------

#ifndef KCOMMON_H
#define KCOMMON_H

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

//-------------------------------------------------------------------
// KGeometryMatrixH:
// Assembles geometry measurements matrix H
// pclPseudorangeH_ (In): Pseudorange Matrix H_rho
// pclPseudorateH_ (In): Pseudorange rate Matrix H_rho^dot
// iNumExtraDim_ (In): Number of extra dimensions for VRW with pseudoranges only
//-------------------------------------------------------------------
CMatrix KGeometryMatrixH(CMatrix* pclPseudorangeH_, CMatrix* pclPseudorateH_=NULL, INT iNumExtraDim_=0);

//-------------------------------------------------------------------
// KPseudorangeMatrixH:
// Compute Geometry design matrix H for pseudo-ranges
// pstEpochData_ (In): Satellite information data per epoch
// pclSysStateGeoid_ (In): Current epoch latitude (radians), longitude (radians) and height (m)
//-------------------------------------------------------------------
CMatrix KPseudorangeMatrixH(EpochInfo* pstEpochData_, CMatrix* pclSysStateGeoid_);

//-------------------------------------------------------------------
// KPseudorateMatrixH:
// Compute Geometry design matrix H for pseudo-range rates
// pstEpochData_ (In): Satellite information data per epoch
// pclSysStateGeoid_ (In): Current epoch state
//-------------------------------------------------------------------
CMatrix KPseudorateMatrixH(EpochInfo* pstEpochData_, CMatrix* pclSysStateGeoid_);

//-------------------------------------------------------------------
// Compute Azimuth and Elevation dependent matrix C
// pstSVData_ (In): Current epoch, Current Satellite information data
//-------------------------------------------------------------------
CMatrix KAzimuthElevationC(SVInfo* pstSVData_);

//-------------------------------------------------------------------
// KPseudorangeGeoidT
// Compute Geoid Matrix based on current latitude, longitude and 
// height estimates
// pstSVData_ (In): Current epoch, Current Satellite information data
// pclSysStateGeoid_ (In): Latitude (radians), longitude (radians) and height (m)
//-------------------------------------------------------------------
CMatrix KPseudorangeGeoidT(SVInfo* pstSVData_, CMatrix* pclSysStateGeoid_);

//-------------------------------------------------------------------
// KPseudorateGeoidTS
// Compute Geoid Velocity Matrix based on current latitude, longitude 
// and height
// pstSVData_ (In): Current epoch, Current Satellite information data
// pclSysStateGeoid_ (In): Latitude (radians), longitude (radians), height (m)
//-------------------------------------------------------------------
CMatrix KPseudorateGeoidTS(SVInfo* pstSVData_, CMatrix* pclSysStateGeoid_);

//------------------------------------------------------------------------------------------------
// Compute transformation matrix from Geoid to Local level frame coordinates
// S
// pclSysStateGeoid_ (In): (Latitude (radians), Longitude (radians), Heigth (m)) states
//------------------------------------------------------------------------------------------------
CMatrix KGeoidToLLFS(CMatrix* pclSysStateGeoid_);

//------------------------------------------------------------------------------------------------
// KResiudals
// Compute measurement residuals for Kalman filter
// pclPseudoRangeResidual_ (In): Vector of pseudo-range residuals
// pclPseudoRateResidual_ (In): Vector of pseudorange-rate residuals
//------------------------------------------------------------------------------------------------
CMatrix KResiudals(CMatrix* pclPseudoRangeResidual_, CMatrix* pclPseudoRateResidual_=NULL);

//------------------------------------------------------------------------------------------------
// KPseudorangeResiduals
// Computes residuals between pseudorange measurements and estimated pseudorange
// pstEpochData_ (In): Satellite geometry per epoch
// pclSysStateEcef_ (In): Estimated ECEF coordinates
//------------------------------------------------------------------------------------------------
CMatrix KPseudorangeResiduals(EpochInfo* pstEpochData_, CMatrix* pclSysStateEcef_);

//------------------------------------------------------------------------------------------------
// KPseudorateResiduals
// Computes residuals between pseudorange rate measurements and estimated pseudorange rate
// pstEpochData_ (In): Satellite geometry per epoch
// pclSysStateEcef_ (In): Estimated ECEF coordinates and velocities
//------------------------------------------------------------------------------------------------
CMatrix KPseudorateResiduals(EpochInfo* pstEpochData_, CMatrix* pclSysStateEcef_);

//------------------------------------------------------------------------------------------------
// KObsCovMatrixR
// Computes observation covariance matrix
// pclObsCovPseudorangeR_ (In): Pseudorange observation covariance matrix
// pclObsCovPseudorateR_ (In): Pseudorange rate observation covariance matrix
//------------------------------------------------------------------------------------------------
CMatrix KObsCovMatrixR(CMatrix* pclObsCovPseudorangeR_, CMatrix* pclObsCovPseudorateR_=NULL);

//------------------------------------------------------------------------------------------------
// KObsCovSubMatrixR
// Computes observation covariance matrix from pseudorange or pseudorange rate observations
// pstEpochData_ (In): Satellite geometry per epoch
// dStd_ (In): Pseudorange standard deviation at zenith
//------------------------------------------------------------------------------------------------
CMatrix KObsCovSubMatrixR(EpochInfo* pstEpochData_, DOUBLE dStd_);

//-------------------------------------------------------------------
// KPredictionLoop
// pclF_: System state space matrix
// pclG_: Noise shaping matrix
// pclQ_: Spectral density matrix
// pclSysState_: System state
// pclSysCovP_: State covariance matrix
// dStepSize_: Prediction rate
//-------------------------------------------------------------------
void KPredictionLoop(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
                     CMatrix* pclSysState_, CMatrix* pclSysCovP_, DOUBLE dStepSize_);

//-------------------------------------------------------------------
// KProcessNoiseQ
// pclF_: System state space matrix
// pclG_: Noise shaping matrix
// pclQ_: Spectral density matrix
// dStepSize_: Prediction rate
// ulNumSubIntervals_: Number of sub-intervals in the prediction rate
//-------------------------------------------------------------------
CMatrix KProcessNoiseQ(CMatrix* pclF_, CMatrix* pclG_, CMatrix* pclQ_,
                       DOUBLE dStepSize_, ULONG ulNumSubIntervals_ = 10);

//-------------------------------------------------------------------
// KSpectralDensityQ
// dRWSpectralQ_ (In): Spectral densities of position or velocity RW
// dClkBiasQ_ (In): Spectral density of the Clock bias noise
// dClkDriftQ_ (In): Spectral density of the Clock drift noise
//-------------------------------------------------------------------
CMatrix KSpectralDensityQ(DOUBLE* pdRWSpectralQ_, DOUBLE* pdClkBiasQ_, DOUBLE* pdClkDriftQ_);

#endif