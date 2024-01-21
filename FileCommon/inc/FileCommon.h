//----------------------------------------------------------------------
// FileCommon.h
// Module used to read binary files given a predefined format
//----------------------------------------------------------------------

#ifndef FILECOMMON_H
#define FILECOMMON_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <Windows.h>
#include "..\..\DataTypeDefinitions.h"
#include "..\..\CoordTrans\inc\CoordTrans.h"

using namespace std;

#define FILECOMMON_GPSTIME_TOLERANCE ((DOUBLE) 0.0001)

//----------------------------------------------------------------------------
// PopulateSVEpochData: Reads binary formatted file and populates Epoch vector
// sFilePath_ (In): Full path of binary input file
// vEpochInfo_ (Out): Vector with information of Satellite data epoch by epoch
//----------------------------------------------------------------------------
void PopulateSVEpochData(string sFilePath_, vector<EpochInfo>* pvEpochInfo_);

//----------------------------------------------------------------------------
// SendSVEpochDatatoText: Reads binary formatted file and outputs to text file
// sInFilePath_ (In): Full path of binary input file
// sOutFilePath_ (Out): Full path of output text file
//----------------------------------------------------------------------------
void SendSVEpochDatatoText(string sInFilePath_, string sOutFilePath_);

//----------------------------------------------------------------------------
void TransformReliabilitySV(vector<ReliablityTest>* pvRelTest_, vector<ReliabilitySVInfo>* pvRelSVTest_, CHAR cPRN_);

//----------------------------------------------------------------------------
void SetRootDirectory(string sInFilePath_, string* psOutRootDir_);

//----------------------------------------------------------------------------
void PrintPlots(string sInFilePath_, vector<ErrorInfo>* pvPosError_, vector<ReliablityTest>* pvRelTest_, string sMethod_);

//----------------------------------------------------------------------------
// ParseLine: Separate text file line into meaningful fields
//----------------------------------------------------------------------------
void ParseLine(string* psLine_, TrajectoryInfo* pstTrajInfo_);

//----------------------------------------------------------------------------
// AlignEpochs: Align truth and data file epochs
//----------------------------------------------------------------------------
void AlignEpochs(ifstream& clTruthFile_, vector<EpochInfo>::iterator& itEpoch_,
                 TrajectoryInfo* pstTrajInfo_, vector<EpochInfo>* pvEpochData_);

//------------------------------------------------------------------------------------------------
// Fill trajectory values
// pstTrajInfo_ (Out): Current epoch trajectory values
// dLatitude_ (In): Current latitude
// dLongitude_ (In): Current longitude
// dHeight_ (In): Current Height
// dGPSTime_ (In): Current epoch
// dNorthVel_ (In): Current north velocity
// dEastVel_ (In): Current east velocity
// dHeightVel_ (In): Current height velocity
//------------------------------------------------------------------------------------------------
void FillTrajectory(TrajectoryInfo* pstTrajInfo_,
                    DOUBLE dLatitude_, DOUBLE dLongitude_, DOUBLE dHeight_, DOUBLE dGPSTime_,
                    DOUBLE dNorthVel_ = 0, DOUBLE dEastVel_ = 0, DOUBLE dHeightVel_ = 0);

//------------------------------------------------------------------------------------------------
void KalmanPosErrorsPlots(vector<TrajectoryInfo>* pvTruthTraj_, vector<TrajectoryInfo>* pvKalmanTraj_,
                          string sInFilePath_, string sMethod_);

//------------------------------------------------------------------------------------------------
void KalmanVelErrorsPlots(vector<TrajectoryInfo>* pvTruthTraj_, vector<TrajectoryInfo>* pvKalmanTraj_,
                          string sInFilePath_, string sMethod_);

//------------------------------------------------------------------------------------------------
void KalmanTrajectory(vector<TrajectoryInfo>* pvTraj_, string sInFilePath_, string sMethod_, string sFileName_);

//------------------------------------------------------------------------------------------------
void KalmanCrossCorrPlots(void* pvErrorInfo_, string sInFilePath_, string sMethod_);

#endif