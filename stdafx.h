// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"
#include "DataTypeDefinitions.h"
#include "FileCommon\inc\FileCommon.h"
#include "Matrix\inc\Matrix.h"
#include "CoordTrans\inc\CoordTrans.h"
#include "LeastSquares\inc\LeastSquares.h"
#include "Kalman\inc\Kalman.h"
#include "Kalman\inc\KPolar.h"
#include "Distributions\inc\Distributions.h"

#include <stdio.h>
#include <tchar.h>
#include <string>

// Number of defined options
#define ESTIMATION_NUM_OPTIONS ((ULONG) 31)
// Options position indices
#define ESTIMATION_CHANGELATITUDE_INDEX ((ULONG) 0)
#define ESTIMATION_CHANGELONGITUDE_INDEX ((ULONG) 1)
#define ESTIMATION_CHANGEHEIGHT_INDEX ((ULONG) 2)
#define ESTIMATION_CHANGETIMESHIFT_INDEX ((ULONG) 3)
#define ESTIMATION_CHANGEFILE_INDEX ((ULONG) 4)
#define ESTIMATION_ECEFENABLE_INDEX ((ULONG) 5)
#define ESTIMATION_GEOIDENABLE_INDEX ((ULONG) 6)
#define ESTIMATION_GEOIDCONSTRAINEDENABLE_INDEX ((ULONG) 7)
#define ESTIMATION_COVCONST_INDEX ((ULONG) 8)
#define ESTIMATION_COVELEV_INDEX ((ULONG) 9)
#define ESTIMATION_CONFIDENCELEVEL_INDEX ((ULONG) 10)
#define ESTIMATION_TESTPOWER_INDEX ((ULONG) 11)
#define ESTIMATION_HEIGHTCNST_INDEX ((ULONG) 12)
#define ESTIMATION_TRUTHFILE_INDEX ((ULONG) 13)
#define ESTIMATION_CHANGEVEAST_INDEX ((ULONG) 14)
#define ESTIMATION_CHANGEVNORTH_INDEX ((ULONG) 15)
#define ESTIMATION_CHANGEVVERT_INDEX ((ULONG) 16)
#define ESTIMATION_CHANGECLOCKDRIFT_INDEX ((ULONG) 17)
#define ESTIMATION_COVVELEV_INDEX ((ULONG) 18)
#define ESTIMATION_SPECTRALPRWQ_INDEX ((ULONG) 19)
#define ESTIMATION_SPECTRALVRWQ_INDEX ((ULONG) 20)
#define ESTIMATION_SPECTRALTSHQ_INDEX ((ULONG) 21)
#define ESTIMATION_SPECTRALTDRIFTQ_INDEX ((ULONG) 22)
#define ESTIMATION_INITPOSSTD_INDEX ((ULONG) 23)
#define ESTIMATION_INITVELSTD_INDEX ((ULONG) 24)
#define ESTIMATION_VRWENABLE_INDEX ((ULONG) 25)
#define ESTIMATION_DOPPLERENABLE_INDEX ((ULONG) 26)
#define ESTIMATION_PRWENABLE_INDEX ((ULONG) 27)
#define ESTIMATION_BLUNDERENABLE_INDEX ((ULONG) 28)
#define ESTIMATION_POLARENABLE_INDEX ((ULONG) 29)
#define ESTIMATION_POLARRADIUS_INDEX ((ULONG) 30)