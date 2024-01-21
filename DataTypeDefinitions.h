#ifndef DATATYPEDEFINITIONS_H
#define DATATYPEDEFINITIONS_H

#include <vector>

// Data type definitions
typedef double DOUBLE;
typedef float FLOAT;
typedef int INT;
typedef unsigned int UINT;
typedef unsigned long ULONG;
typedef unsigned long long ULONGLONG;
typedef long LONG;
typedef long long LONGLONG;
typedef short SHORT;
typedef unsigned short USHORT;
typedef unsigned char UCHAR;
typedef char CHAR;
typedef bool BOOLEANO;

// Constant definitions
// Speed of light in vacuum, m/s
#define SPEED_OF_LIGHT ((DOUBLE) 299792458.0)
#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif
#define DATATYPEDEFINITIONS_NUMSTATES ((INT) 4)
#define DATATYPEDEFINITIONS_NUMPRN ((INT) 33)

// Structure definitions
typedef struct SVInfo
{
   CHAR   cPRN;       // Satellite PRN number
   DOUBLE dPSR;       // Satellite pseudorange
   DOUBLE dPSRRate;   // Satellite pseudorange rate
   DOUBLE dAzimuth;   // Satellite azimuth
   DOUBLE dElevation; // Satellite elevation
   DOUBLE dXPosition; // Satellite X cartesian position
   DOUBLE dYPosition; // Satellite Y cartesian position
   DOUBLE dZPosition; // Satellite Z cartesian position
   DOUBLE dXVelocity; // Satellite X cartesian velocity
   DOUBLE dYVelocity; // Satellite Y cartesian velocity
   DOUBLE dZVelocity; // Satellite Z cartesian velocity
} SVInfo;

typedef struct EpochInfo
{
   DOUBLE dGPSTime;               // GPS Time
   CHAR   cNumSV;                 // Number of Satellites
   std::vector<SVInfo> vstSVInfo; // Vector with satellite data
} EpochInfo;

typedef struct ErrorInfo
{
   CHAR   cNumSV;        // Number of satellites
   DOUBLE dGPSTime;      // GPS Time
   DOUBLE dEastError;    // Truth east error
   DOUBLE dNorthError;   // Truth north error
   DOUBLE dHeightError;  // Truth Height error
   DOUBLE dEastStd;      // Estimated east error standard deviation
   DOUBLE dNorthStd;     // Estimated north error standard deviation
   DOUBLE dHeightStd;    // Estimated height error standard deviation
   DOUBLE dHDop;         // Horizontal dilution of precision
   DOUBLE dVDop;         // Vertical dilution of precision
   DOUBLE dPDop;         // 3-D dilution of precision
   DOUBLE dENCorr;       // E-N Correlation coefficient
   DOUBLE dEUCorr;       // E-U Correlation coefficient
   DOUBLE dEClkCorr;     // E-Clock Correlation coefficient
   DOUBLE dNUCorr;       // N-U correlation coefficient
   DOUBLE dNClkCorr;     // N-Clock Correlation coefficient
   DOUBLE dUClkCorr;     // U-Clock Correlation coefficient
} ErrorInfo;

typedef struct PRWErrorInfo
{
   CHAR   cNumSV;           // Number of satellites
   DOUBLE dGPSTime;         // GPS Time
   DOUBLE dEastError;       // Truth east error
   DOUBLE dNorthError;      // Truth north error
   DOUBLE dHeightError;     // Truth Height error
   DOUBLE dEastStd;         // Estimated east error standard deviation
   DOUBLE dNorthStd;        // Estimated north error standard deviation
   DOUBLE dHeightStd;       // Estimated height error standard deviation
   DOUBLE dENCorr;          // E-N Correlation coefficient
   DOUBLE dEUCorr;          // E-U Correlation coefficient
   DOUBLE dEClkCorr;        // E-Clock Correlation coefficient
   DOUBLE dEClkDriftCorr;   // E-Clock drift correlation coefficient
   DOUBLE dNUCorr;          // N-U correlation coefficient
   DOUBLE dNClkCorr;        // N-Clock Correlation coefficient
   DOUBLE dNClkDriftCorr;   // N-Clock drift correlation coefficient
   DOUBLE dUClkCorr;        // U-Clock Correlation coefficient
   DOUBLE dUClkDriftCorr;   // U-Clock drift correlation coefficient
   DOUBLE dClkClkDriftCorr; // Clock bias-Clock drift correlation coefficient
} PRWErrorInfo;

typedef struct VRWErrorInfo
{
   CHAR   cNumSV;           // Number of satellites
   DOUBLE dGPSTime;         // GPS Time
   DOUBLE dEastError;       // Truth east error
   DOUBLE dNorthError;      // Truth north error
   DOUBLE dHeightError;     // Truth Height error
   DOUBLE dVelEastError;    // Truth east velocity error
   DOUBLE dVelNorthError;   // Truth north velocity error
   DOUBLE dVelUpError;      // Truth vertical velocity error
   DOUBLE dEastStd;         // Estimated east error standard deviation
   DOUBLE dNorthStd;        // Estimated north error standard deviation
   DOUBLE dHeightStd;       // Estimated height error standard deviation
   DOUBLE dVelEastStd;      // Estimated velocity east error standard deviation
   DOUBLE dVelNorthStd;     // Estimated velocity north error standard deviation
   DOUBLE dVelUpStd;        // Estimated velocity up error standard deviation
   DOUBLE dENCorr;          // E-N Correlation coefficient
   DOUBLE dEUCorr;          // E-U Correlation coefficient
   DOUBLE dEVECorr;         // E-Vel E Correlation coefficient
   DOUBLE dEVNCorr;         // E-Vel N Correlation coefficient
   DOUBLE dEVUCorr;         // E-Vel U Correlation coefficient
   DOUBLE dEClkCorr;        // E-Clock Correlation coefficient
   DOUBLE dEClkDriftCorr;   // E-Clock drift correlation coefficient
   DOUBLE dNUCorr;          // N-U correlation coefficient
   DOUBLE dNVECorr;         // N-Vel E Correlation coefficient
   DOUBLE dNVNCorr;         // N-Vel N Correlation coefficient
   DOUBLE dNVUCorr;         // N-Vel U Correlation coefficient
   DOUBLE dNClkCorr;        // N-Clock Correlation coefficient
   DOUBLE dNClkDriftCorr;   // N-Clock drift correlation coefficient
   DOUBLE dUVECorr;         // U-Vel E Correlation coefficient
   DOUBLE dUVNCorr;         // U-Vel N Correlation coefficient
   DOUBLE dUVUCorr;         // U-Vel U Correlation coefficient
   DOUBLE dUClkCorr;        // U-Clock Correlation coefficient
   DOUBLE dUClkDriftCorr;   // U-Clock drift correlation coefficient
   DOUBLE dVEVNCorr;        // Vel E-Vel N Correlation coefficient
   DOUBLE dVEVUCorr;        // Vel E-Vel U Correlation coefficient
   DOUBLE dVEClkCorr;       // Vel E-Clock bias correlation coefficient
   DOUBLE dVEClkDriftCorr;  // Vel E-Clock drift correlation coefficient
   DOUBLE dVNVUCorr;        // Vel N-Vel U Correlation coefficient
   DOUBLE dVNClkCorr;       // Vel N-Clock bias correlation coefficient
   DOUBLE dVNClkDriftCorr;  // Vel N-Clock drift correlation coefficient
   DOUBLE dVUClkCorr;       // Vel U-Clock bias correlation coefficient
   DOUBLE dVUClkDriftCorr;  // Vel U-Clock drift correlation coefficient
   DOUBLE dClkClkDriftCorr; // Clock bias-Clock drift correlation coefficient
} VRWErrorInfo;

typedef struct TrajectoryInfo
{
	DOUBLE dGPSTime;     // Epoch time
	DOUBLE dLatitude;    // Truth latitude
	DOUBLE dLongitude;   // Truth longitude
	DOUBLE dHeight;      // Truth height
   DOUBLE dNorthVel;    // Truth North Velocity
   DOUBLE dEastVel;     // Truth East Velocity
   DOUBLE dHeightVel;   // Truth Vertical Velocity
} TrajectoryInfo;

typedef struct ReliabilityInfo
{
   CHAR cPRN;            // Satellite PRN
   DOUBLE dElevation;    // Elevation angle
   DOUBLE dAzimuth;      // Azimuth angle
   DOUBLE dRelIntern;    // Internal Reliability
   DOUBLE dRelExternX;   // External Reliability ECEF X
   DOUBLE dRelExternY;   // External Reliability ECEF Y
   DOUBLE dRelExternZ;   // External Reliability ECEF Z
   DOUBLE dRelExternE;   // External Reliability East
   DOUBLE dRelExternN;   // External Reliability North
   DOUBLE dRelExternU;   // External Reliability Up
   DOUBLE dRelExternClk; // External Reliability Clock
   DOUBLE dResidual;
} ReliabilityInfo;

typedef struct ReliabilitySVInfo
{
   DOUBLE dGPSTime;                   // Current epoch time
   ReliabilityInfo stReliabilityInfo; // Satellite reliability info
} ReliabilitySVInfo;

typedef struct ReliablityTest
{
   DOUBLE dGPSTime;                           // Epoch time
   CHAR cNumSV;                               // Number of satellites in epoch
   DOUBLE dGlobVar;                           // Global variance test value
   DOUBLE dHighBound;                         // Chi-square high significance limit
   DOUBLE dLowBound;                          // Chi-square low significance limit
   BOOLEANO bGlobTest;                        // Global variance pass/fail
   std::vector<ReliabilityInfo> vReliability; // Vector of internal and external reliabilities
} ReliablityTest;

#endif