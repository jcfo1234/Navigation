//----------------------------------------------------------------------
// LeastSquares.h
// Module that computes position with Least Squares
//----------------------------------------------------------------------

#ifndef LEASTSQUARES_H
#define LEASTSQUARES_H

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
#include "ecefls.h"
#include "geoidls.h"

//----------------------------------------------------------------------------
// LeastSquaresSolution: Compute least square solutions
//----------------------------------------------------------------------------
void LeastSquaresSolution(BOOLEANO bTheECEFLS_, BOOLEANO bTheGeoidLS_, BOOLEANO bTheGeoidConstrainedLS_,
                          DOUBLE* pdTheLatitude_, DOUBLE* pdTheLongitude_, DOUBLE* pdTheHeight_, DOUBLE* pdTheTShift_,
                          DOUBLE* pdTheTruthLatitude_, DOUBLE* pdTheTruthLongitude_, DOUBLE* pdTheTruthHeight_,
                          string* psTheTruthFullFilePath_, Distributions* pclTheProbabilityDistributions_, string* psTheObsCovKey_,
                          string* psTheFullFilePath_, DOUBLE* pdTheMeasStd_, DOUBLE* pdTheConfidenceLevel_, DOUBLE* pdTheTestPower_,
                          DOUBLE* pdTheHeightConstrained_, DOUBLE* pdTheHeightObsStd_, BOOLEANO bTheBlunderEnable_);

#endif