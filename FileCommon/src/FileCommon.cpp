//----------------------------------------------------------------------
// FileCommon.cpp
// Source file of FileCommon.h
//----------------------------------------------------------------------
#include "..\inc\FileCommon.h"

//----------------------------------------------------------------------------
void ParseLine(string* psLine_, TrajectoryInfo* pstTrajInfo_)
{
   istringstream ss(*psLine_);
   string sField;
   // Flags indicating field was found
   BOOL bGPSTimeFound = FALSE;
   BOOL bLatitudeFound = FALSE;
   BOOL bLongitudeFound = FALSE;
   BOOL bHeightFound = FALSE;
   BOOL bNorthVelFound = FALSE;
   BOOL bEastVelFound = FALSE;
   BOOL bHeightVelFound = FALSE;

   // Read line and separate fields
   while (getline(ss, sField, ' '))
   {
      // Remove unwanted spaces
      if (!sField.empty())
      {
         // GPS Time Field
         if (!bGPSTimeFound)
         {
            pstTrajInfo_->dGPSTime = stod(sField.data());
            bGPSTimeFound = TRUE;
            continue;
         }
         // Latitude Field
         if (!bLatitudeFound)
         {
            pstTrajInfo_->dLatitude = DegToRad(stod(sField.data()));
            bLatitudeFound = TRUE;
            continue;
         }
         // Longitude Field
         if (!bLongitudeFound)
         {
            pstTrajInfo_->dLongitude = DegToRad(stod(sField.data()));
            bLongitudeFound = TRUE;
            continue;
         }
         // Height Field
         if (!bHeightFound)
         {
            pstTrajInfo_->dHeight = stod(sField.data());
            bHeightFound = TRUE;
            continue;
         }
         // North Velocity Field
         if (!bNorthVelFound)
         {
            pstTrajInfo_->dNorthVel = stod(sField.data());
            bNorthVelFound = TRUE;
            continue;
         }
         // East Velocity Field
         if (!bEastVelFound)
         {
            pstTrajInfo_->dEastVel = stod(sField.data());
            bEastVelFound = TRUE;
            continue;
         }
         // Up velocity Field
         if (!bHeightVelFound)
         {
            pstTrajInfo_->dHeightVel = stod(sField.data());
            bHeightVelFound = TRUE;
            continue;
         }
      }
   }
}

//----------------------------------------------------------------------------
void AlignEpochs(ifstream& clTruthFile_, vector<EpochInfo>::iterator& itEpoch_,
                 TrajectoryInfo* pstTrajInfo_, vector<EpochInfo>* pvEpochData_)
{
   // Truth Data variables
   string sFileRead("");
   string sWord("");

   // Read first line of file
   getline(clTruthFile_, sFileRead);
   ParseLine(&sFileRead, pstTrajInfo_);
   while (abs(itEpoch_->dGPSTime - pstTrajInfo_->dGPSTime) > FILECOMMON_GPSTIME_TOLERANCE && !clTruthFile_.eof() && itEpoch_ != pvEpochData_->end())
   {
      // Align truth epoch to data file epoch
      if ((itEpoch_->dGPSTime - pstTrajInfo_->dGPSTime) > FILECOMMON_GPSTIME_TOLERANCE)
      {
         clTruthFile_ >> sFileRead;
         ParseLine(&sFileRead, pstTrajInfo_);
      }
      // Align data file epoch to truth epoch
      else if ((pstTrajInfo_->dGPSTime - itEpoch_->dGPSTime) > FILECOMMON_GPSTIME_TOLERANCE)
      {
         itEpoch_++;
      }
   }
}

//----------------------------------------------------------------------------
void PopulateSVEpochData(string sFilePath_, vector<EpochInfo>* pvEpochInfo_)
{
   ifstream clBinFile;
   SVInfo stSatelliteData;
   EpochInfo stEpochData;

   clBinFile.open(sFilePath_, ios::binary);

   //# Format of data is:
   //  GPS Time (DOUBLE)
   //  Number of Satellites (UCHAR)
   //  Satellite ID PRN (UCHAR)
   //  Pseudorange (DOUBLE)
   //  Pseudorange Rate (DOUBLE)
   //  Satellite Azimuth (DOUBLE)
   //  Satellite Elevation (DOUBLE)
   //  Satellite X Position (DOUBLE)
   //  Satellite Y Position (DOUBLE)
   //  Satellite Z Position (DOUBLE)
   //  Satellite X Velocity (DOUBLE)
   //  Satellite Y Velocity (DOUBLE)
   //  Satellite Z Velocity (DOUBLE)
   if (clBinFile.is_open())
   {
      // Beginning of file
      clBinFile.seekg(0, clBinFile.beg);
      while (clBinFile)
      {
         stEpochData.vstSVInfo.clear();
         // Extract GPS Time (Epoch)
         clBinFile.read((CHAR*)&stEpochData.dGPSTime, sizeof(DOUBLE));
         // Extract Number of Satellites (Epoch)
         clBinFile.read(&stEpochData.cNumSV, sizeof(CHAR));
         // Extract Satellite data Info (SV Data)
         for (INT iIndex = 0; iIndex < (INT)stEpochData.cNumSV; iIndex++)
         {
            // Extract PRN Number (SV Data)
            clBinFile.read(&stSatelliteData.cPRN, sizeof(CHAR));
            // Extract Pseudorange (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dPSR, sizeof(DOUBLE));
            // Extract Pseudorange rate (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dPSRRate, sizeof(DOUBLE));
            // Extract Azimuth (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dAzimuth, sizeof(DOUBLE));
            // Extract Elevation (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dElevation, sizeof(DOUBLE));
            // Extract Satellite X Position (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dXPosition, sizeof(DOUBLE));
            // Extract Satellite Y Position (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dYPosition, sizeof(DOUBLE));
            // Extract Satellite Z Position (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dZPosition, sizeof(DOUBLE));
            // Extract Satellite X Velocity (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dXVelocity, sizeof(DOUBLE));
            // Extract Satellite Y Velocity (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dYVelocity, sizeof(DOUBLE));
            // Extract Satellite Z Velocity (SV Data)
            clBinFile.read((CHAR*)&stSatelliteData.dZVelocity, sizeof(DOUBLE));
            // Verify read operation had not reached end of file
            if (clBinFile)
               stEpochData.vstSVInfo.push_back(stSatelliteData);
            else
               break;
         }
         // Verify read operation had not reached end of file
         if (clBinFile)
            pvEpochInfo_->push_back(stEpochData);
         cout << fixed << setprecision(3) << stEpochData.dGPSTime << endl;
      }
   }

   clBinFile.close();
}

//----------------------------------------------------------------------------
void SendSVEpochDatatoText(string sInFilePath_, string sOutFilePath_)
{
   ifstream clBinFile;
   ofstream clTextFile;
   CHAR cBuffer;
   CHAR cNumSV;
   DOUBLE dBuffer;
   DOUBLE dGPSTime;

   clBinFile.open(sInFilePath_, ios::binary);
   SetRootDirectory(sInFilePath_, &sOutFilePath_);
   sOutFilePath_.append("\\Data.txt");
   clTextFile.open(sOutFilePath_);

   if (clBinFile.is_open())
   {
      // Beginning of file
      clBinFile.seekg(0, clBinFile.beg);
      while (clBinFile)
      {
         // Extract GPS Time (Epoch)
         clBinFile.read((CHAR*)&dGPSTime, sizeof(DOUBLE));
         clTextFile << fixed << setprecision(3) << dGPSTime << ",";
         // Extract Number of Satellites (Epoch)
         clBinFile.read(&cNumSV, sizeof(CHAR));
         clTextFile << (INT)cNumSV << ",";
         // Extract Satellite data Info (SV Data)
         for (INT iIndex = 0; iIndex < (INT)cNumSV; iIndex++)
         {
            // Extract PRN Number (SV Data)
            clBinFile.read(&cBuffer, sizeof(CHAR));
            clTextFile << (INT)cBuffer << ",";
            // Extract Pseudorange (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Pseudorange rate (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Azimuth (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Elevation (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite X Position (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite Y Position (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite Z Position (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite X Velocity (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite Y Velocity (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Extract Satellite Z Velocity (SV Data)
            clBinFile.read((CHAR*)&dBuffer, sizeof(DOUBLE));
            clTextFile << fixed << setprecision(3) << dBuffer << ",";
            // Verify read operation had not reached end of file
            if (!clBinFile)
               break;
         }

         if (clBinFile)
            clTextFile << endl;
         else
            break;
         cout << fixed << setprecision(3) << dGPSTime << endl;
      }
   }
}

//----------------------------------------------------------------------------
void TransformReliabilitySV(vector<ReliablityTest>* pvRelTest_, vector<ReliabilitySVInfo>* pvRelSVTest_, CHAR cPRN_)
{
   ReliabilitySVInfo stSVInfoTimeSeries;
   // Do a PRN Time series of reliability measurements
   for (vector<ReliablityTest>::iterator itRelTest = pvRelTest_->begin(); itRelTest != pvRelTest_->end(); itRelTest++)
   {
      stSVInfoTimeSeries.dGPSTime = itRelTest->dGPSTime;
      for (vector<ReliabilityInfo>::iterator itRelSV = itRelTest->vReliability.begin(); itRelSV < itRelTest->vReliability.end(); itRelSV++)
      {
         // Match requested PRN
         if ((INT)itRelSV->cPRN == (INT)cPRN_)
         {
            stSVInfoTimeSeries.stReliabilityInfo.dRelIntern = itRelSV->dRelIntern;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternX = itRelSV->dRelExternX;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternY = itRelSV->dRelExternY;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternZ = itRelSV->dRelExternZ;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternE = itRelSV->dRelExternE;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternN = itRelSV->dRelExternN;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternU = itRelSV->dRelExternU;
            stSVInfoTimeSeries.stReliabilityInfo.dRelExternClk = itRelSV->dRelExternClk;
            stSVInfoTimeSeries.stReliabilityInfo.dResidual = itRelSV->dResidual;
            stSVInfoTimeSeries.stReliabilityInfo.cPRN = cPRN_;
            stSVInfoTimeSeries.stReliabilityInfo.dAzimuth = itRelSV->dAzimuth;
            stSVInfoTimeSeries.stReliabilityInfo.dElevation = itRelSV->dElevation;
            pvRelSVTest_->push_back(stSVInfoTimeSeries);
            break;
         }
      }
   }
}

//----------------------------------------------------------------------------
void SetRootDirectory(string sInFilePath_, string* psOutRootDir_)
{
   string sTemp;
   INT iPos = 0;
   psOutRootDir_->clear();
   while (sInFilePath_.find("\\", iPos) != string::npos)
   {
      iPos = sInFilePath_.find("\\", iPos) + 1;
   }
   sTemp = sInFilePath_.substr(0, iPos - 1);
   psOutRootDir_->append(sTemp);
}

//------------------------------------------------------------------------------------------------
void FillTrajectory(TrajectoryInfo* pstTrajInfo_,
                    DOUBLE dLatitude_, DOUBLE dLongitude_, DOUBLE dHeight_, DOUBLE dGPSTime_,
                    DOUBLE dNorthVel_, DOUBLE dEastVel_, DOUBLE dHeightVel_)
{
   pstTrajInfo_->dGPSTime = dGPSTime_;
   pstTrajInfo_->dLatitude = DegToRad(dLatitude_);
   pstTrajInfo_->dLongitude = DegToRad(dLongitude_);
   pstTrajInfo_->dHeight = dHeight_;
   pstTrajInfo_->dNorthVel = dNorthVel_;
   pstTrajInfo_->dEastVel = dEastVel_;
   pstTrajInfo_->dHeightVel = dHeightVel_;
}

//----------------------------------------------------------------------------
void PrintPlots(string sInFilePath_, vector<ErrorInfo>* pvPosError_, vector<ReliablityTest>* pvRelTest_, string sMethod_)
{
   string sRootDirectory;
   string sOutFile;
   string sOutFile1;
   ofstream clTextFile;
   ofstream clTextFile1;
   vector<ReliabilitySVInfo> avRelPRNTimeSeries[33];

   // Set the root directory where all result files are going to be dumped
   SetRootDirectory(sInFilePath_, &sRootDirectory);
   sRootDirectory.append("\\").append(sMethod_);
   CreateDirectory(sRootDirectory.c_str(), NULL);

   // Print global variance test results
   sOutFile = sRootDirectory;
   sOutFile.append("\\GlobalVariance").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Global Variance Test" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=Variance (m^2)" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=4" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // Global Variance Plot
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=Global Variance" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Low variance bound plot
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=Low Variance Bound" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // High variance bound plot
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=High Variance Bound" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Test variance result
   clTextFile << "\n[series4]" << endl;
   clTextFile << "name=Test Variance Result" << endl;
   clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;

   for (vector<ReliablityTest>::iterator itRel = pvRelTest_->begin(); itRel != pvRelTest_->end(); itRel++)
   {
      clTextFile << fixed << setprecision(3) << itRel->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << itRel->dGlobVar << ",";
      clTextFile << fixed << setprecision(6) << itRel->dLowBound << ",";
      clTextFile << fixed << setprecision(6) << itRel->dHighBound << ",";
      clTextFile << fixed << setprecision(6) << itRel->bGlobTest << "," << endl;
   }
   clTextFile.close();

   // Reliability per PRN files
   for (INT i = 1; i <= DATATYPEDEFINITIONS_NUMPRN; i++)
   {
      avRelPRNTimeSeries[i - 1].clear();
      TransformReliabilitySV(pvRelTest_, &avRelPRNTimeSeries[i - 1], i);
      // If PRN was tracked at all
      if (!avRelPRNTimeSeries[i - 1].empty())
      {
         sOutFile = sRootDirectory;
         sOutFile1 = sRootDirectory;
         sOutFile.append("\\ReliabilityPRN").append(to_string(i)).append(".GIT");
         sOutFile1.append("\\PRNAzimuthElevation").append(to_string(i)).append(".GIT");
         clTextFile.open(sOutFile);
         clTextFile1.open(sOutFile1);
         clTextFile << "[header]" << endl;
         clTextFile << "title=Reliability Measurements for PRN " << to_string(i) << endl;
         clTextFile << "xlabel=Time(s)" << endl;
         clTextFile << "ylabel=Reliability Test Results (m)" << endl;
         clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=9" << endl << "miny=-10" << endl << "maxy=10" << endl;
         // Internal Reliability Plot
         clTextFile << "\n[series1]" << endl;
         clTextFile << "name=Internal Reliability" << endl;
         clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability ECEF X Plot
         clTextFile << "\n[series2]" << endl;
         clTextFile << "name=External Reliability ECEF X" << endl;
         clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability ECEF Y Plot
         clTextFile << "\n[series3]" << endl;
         clTextFile << "name=External Reliability ECEF Y" << endl;
         clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability ECEF Z Plot
         clTextFile << "\n[series4]" << endl;
         clTextFile << "name=External Reliability ECEF Z" << endl;
         clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability LLF E Plot
         clTextFile << "\n[series5]" << endl;
         clTextFile << "name=External Reliability LLF E" << endl;
         clTextFile << "column=5" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability LLF N Plot
         clTextFile << "\n[series6]" << endl;
         clTextFile << "name=External Reliability LLF N" << endl;
         clTextFile << "column=6" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability LLF U Plot
         clTextFile << "\n[series7]" << endl;
         clTextFile << "name=External Reliability LLF U" << endl;
         clTextFile << "column=7" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // External Reliability Clock Plot
         clTextFile << "\n[series8]" << endl;
         clTextFile << "name=External Reliability LLF Clk" << endl;
         clTextFile << "column=8" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // Pseudo-range residual
         clTextFile << "\n[series9]" << endl;
         clTextFile << "name=Pseudorange Residual" << endl;
         clTextFile << "column=9" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         clTextFile << "\n[data]" << endl;
         // Azimuth and Elevation angles
         clTextFile1 << "[header]" << endl;
         clTextFile1 << "title=Azimuth and Elevation for PRN " << to_string(i) << endl;
         clTextFile1 << "xlabel=Time(s)" << endl;
         clTextFile1 << "ylabel=Degrees" << endl;
         clTextFile1 << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=2" << endl << "miny=-10" << endl << "maxy=10" << endl;
         // Azimuth series
         clTextFile1 << "\n[series1]" << endl;
         clTextFile1 << "name=Azimuth Angle" << endl;
         clTextFile1 << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         // Elevation series
         clTextFile1 << "\n[series2]" << endl;
         clTextFile1 << "name=Elevation Angle" << endl;
         clTextFile1 << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
         clTextFile1 << "\n[data]" << endl;

         for (vector<ReliabilitySVInfo>::iterator itRelPRN_ = avRelPRNTimeSeries[i - 1].begin(); itRelPRN_ != avRelPRNTimeSeries[i - 1].end(); itRelPRN_++)
         {
            // Reliability files
            clTextFile << fixed << setprecision(3) << itRelPRN_->dGPSTime << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelIntern << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternX << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternY << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternZ << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternE << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternN << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternU << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dRelExternClk << ",";
            clTextFile << fixed << setprecision(6) << itRelPRN_->stReliabilityInfo.dResidual << "," << endl;
            // Azimuth and elevation angle files
            clTextFile1 << fixed << setprecision(3) << itRelPRN_->dGPSTime << ",";
            clTextFile1 << fixed << setprecision(3) << itRelPRN_->stReliabilityInfo.dAzimuth << ",";
            clTextFile1 << fixed << setprecision(3) << itRelPRN_->stReliabilityInfo.dElevation << "," << endl;
         }
         clTextFile.close();
         clTextFile1.close();
      }
   }

   // Position errors and standard deviations
   sOutFile = sRootDirectory;
   sOutFile.append("\\PositionErrors").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Position Error Stats" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=Position Errors Stats (m)" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=7" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // East position error
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=East-Error" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North position error
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=North-Error" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up position error
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=Up-Error" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East standard deviation
   clTextFile << "\n[series4]" << endl;
   clTextFile << "name=East-StdDev" << endl;
   clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North standard deviation
   clTextFile << "\n[series5]" << endl;
   clTextFile << "name=North-StdDev" << endl;
   clTextFile << "column=5" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up standard deviation
   clTextFile << "\n[series6]" << endl;
   clTextFile << "name=Up-StdDev" << endl;
   clTextFile << "column=6" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Number of satellites
   clTextFile << "\n[series7]" << endl;
   clTextFile << "name=Num-SV" << endl;
   clTextFile << "column=7" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   for (vector<ErrorInfo>::iterator itPosError = pvPosError_->begin(); itPosError != pvPosError_->end(); itPosError++)
   {
      clTextFile << fixed << setprecision(3) << itPosError->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dEastError << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dNorthError << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dHeightError << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dEastStd << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dNorthStd << ",";
      clTextFile << fixed << setprecision(6) << itPosError->dHeightStd << ",";
      clTextFile << to_string(itPosError->cNumSV) << endl;
   }
   clTextFile.close();

   // Dilution of precision parameters
   sOutFile = sRootDirectory;
   sOutFile.append("\\DOP").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Dilution of Precision" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=DOP" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=4" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // Horizontal DOP
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=H-DOP" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Vertical DOP
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=V-DOP" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Position DOP
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=P-DOP" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Number of satellites
   clTextFile << "\n[series4]" << endl;
   clTextFile << "name=Num-SV" << endl;
   clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   for (vector<ErrorInfo>::iterator itDOP = pvPosError_->begin(); itDOP != pvPosError_->end(); itDOP++)
   {
      clTextFile << fixed << setprecision(3) << itDOP->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << itDOP->dHDop << ",";
      clTextFile << fixed << setprecision(6) << itDOP->dVDop << ",";
      clTextFile << fixed << setprecision(6) << itDOP->dPDop << ",";
      clTextFile << to_string(itDOP->cNumSV) << endl;
   }
   clTextFile.close();

   // Cross-correlation parameters
   sOutFile = sRootDirectory;
   sOutFile.append("\\CrossCorr").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Cross Correlations" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=XCORR" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=7" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // East-North cross-correlation
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=EN-XCORR" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East-Up cross-correlation
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=EU-XCORR" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East-Clock cross-correlation
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=ECLK-XCORR" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North-Up cross-correlation
   clTextFile << "\n[series4]" << endl;
   clTextFile << "name=NU-XCORR" << endl;
   clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North-Clock cross-correlation
   clTextFile << "\n[series5]" << endl;
   clTextFile << "name=NCLK-XCORR" << endl;
   clTextFile << "column=5" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up-Clock cross-correlation
   clTextFile << "\n[series6]" << endl;
   clTextFile << "name=UCLK-XCORR" << endl;
   clTextFile << "column=6" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Number of satellites
   clTextFile << "\n[series7]" << endl;
   clTextFile << "name=Num-SV" << endl;
   clTextFile << "column=7" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   for (vector<ErrorInfo>::iterator itXcorr = pvPosError_->begin(); itXcorr != pvPosError_->end(); itXcorr++)
   {
      clTextFile << fixed << setprecision(3) << itXcorr->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dENCorr << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dEUCorr << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dEClkCorr << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dNUCorr << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dNClkCorr << ",";
      clTextFile << fixed << setprecision(6) << itXcorr->dUClkCorr << ",";
      clTextFile << to_string(itXcorr->cNumSV) << endl;
   }
   clTextFile.close();

}

//------------------------------------------------------------------------------------------------
void KalmanPosErrorsPlots(vector<TrajectoryInfo>* pvTruthTraj_, vector<TrajectoryInfo>* pvKalmanTraj_,
                          string sInFilePath_, string sMethod_)
{
   string sRootDirectory;
   string sOutFile;
   ofstream clTextFile;

   // Set the root directory where all result files are going to be dumped
   SetRootDirectory(sInFilePath_, &sRootDirectory);
   sRootDirectory.append("\\").append(sMethod_);
   CreateDirectory(sRootDirectory.c_str(), NULL);

   // Position errors
   sOutFile = sRootDirectory;
   sOutFile.append("\\PositionErrors").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Position Error Stats" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=Position Errors Stats (m)" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=3" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // East position error
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=East-Error" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North position error
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=North-Error" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up position error
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=Up-Error" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   vector<TrajectoryInfo>::iterator pIterKalman;
   vector<TrajectoryInfo>::iterator pIterTruth;
   // Compute Local level frame errors
   DOUBLE dLatError = 0;
   DOUBLE dLngError = 0;
   DOUBLE dEastError = 0;
   DOUBLE dNorthError = 0;
   CMatrix clLatSingleton("Latitude", 1, 1);
   CMatrix clN("N", 1, 1);
   CMatrix clM("M", 1, 1);
   // Iterate over the two trajectories to find the errors in position
   for (pIterTruth = pvTruthTraj_->begin(), pIterKalman = pvKalmanTraj_->begin();
        pIterTruth != pvTruthTraj_->end(), pIterKalman != pvKalmanTraj_->end(); 
        pIterTruth++, pIterKalman++)
   {
      // Compute East and North errors
      clLatSingleton.SetComponent(0, 0, pIterTruth->dLatitude);
      clN = CalcN(&clLatSingleton);
      clM = CalcM(&clLatSingleton);
      dLatError = pIterTruth->dLatitude - pIterKalman->dLatitude;
      dLngError = pIterTruth->dLongitude - pIterKalman->dLongitude;
      dEastError = dLngError * (clN.GetComponent(0, 0) + pIterTruth->dHeight) * cos(pIterTruth->dLatitude);
      dNorthError = dLatError * (clM.GetComponent(0, 0) + pIterTruth->dHeight);
      // Output results to PositionErrors.GIT
      clTextFile << fixed << setprecision(3) << pIterTruth->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << dEastError << ",";
      clTextFile << fixed << setprecision(6) << dNorthError << ",";
      clTextFile << fixed << setprecision(6) << pIterTruth->dHeight - pIterKalman->dHeight << endl;
   }
   clTextFile.close();
}

//------------------------------------------------------------------------------------------------
void KalmanVelErrorsPlots(vector<TrajectoryInfo>* pvTruthTraj_, 
                          vector<TrajectoryInfo>* pvKalmanTraj_,
                          string sInFilePath_, string sMethod_)
{
   string sRootDirectory;
   string sOutFile;
   ofstream clTextFile;

   // Set the root directory where all result files are going to be dumped
   SetRootDirectory(sInFilePath_, &sRootDirectory);
   sRootDirectory.append("\\").append(sMethod_);
   CreateDirectory(sRootDirectory.c_str(), NULL);

   // Position errors
   sOutFile = sRootDirectory;
   sOutFile.append("\\VelocityErrors").append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Velocity Error Stats" << endl;
   clTextFile << "xlabel=Time(s)" << endl;
   clTextFile << "ylabel=Velocity Errors Stats (m/s)" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=3" << endl << "miny=-10" << endl << "maxy=10" << endl;
   // East position error
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=VEast-Error" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North position error
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=VNorth-Error" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up position error
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=VUp-Error" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   vector<TrajectoryInfo>::iterator pIterKalman;
   vector<TrajectoryInfo>::iterator pIterTruth;
   // Iterate over the two trajectories to find the errors in position
   for (pIterTruth = pvTruthTraj_->begin(), pIterKalman = pvKalmanTraj_->begin();
        pIterTruth != pvTruthTraj_->end(), pIterKalman != pvKalmanTraj_->end();
        pIterTruth++, pIterKalman++)
   {
      // Output results to VelocityErrors.GIT
      clTextFile << fixed << setprecision(3) << pIterTruth->dGPSTime << ",";
      clTextFile << fixed << setprecision(6) << pIterTruth->dEastVel - pIterKalman->dEastVel << ",";
      clTextFile << fixed << setprecision(6) << pIterTruth->dNorthVel - pIterKalman->dNorthVel << ",";
      clTextFile << fixed << setprecision(6) << pIterTruth->dHeightVel - pIterKalman->dHeightVel << endl;
   }
   clTextFile.close();
}

//------------------------------------------------------------------------------------------------
void KalmanTrajectory(vector<TrajectoryInfo>* pvTraj_, string sInFilePath_, string sMethod_, string sFileName_)
{
   string sRootDirectory;
   string sOutFile;
   ofstream clTextFile;

   // Set the root directory where all result files are going to be dumped
   SetRootDirectory(sInFilePath_, &sRootDirectory);
   sRootDirectory.append("\\").append(sMethod_);
   CreateDirectory(sRootDirectory.c_str(), NULL);

   sOutFile = sRootDirectory;
   sOutFile.append("\\").append(sFileName_).append(".GIT");
   clTextFile.open(sOutFile);
   clTextFile << "[header]" << endl;
   clTextFile << "title=Kalman Truth Latitude vs Longitude" << endl;
   clTextFile << "xlabel=Longitude (Degrees)" << endl;
   clTextFile << "ylabel=Latitude (Degrees)" << endl;
   clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=1" << endl << "miny=-10" << endl << "maxy=10" << endl;
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=Trajectory" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   clTextFile << "\n[data]" << endl;
   for (vector<TrajectoryInfo>::iterator pIterTruth = pvTraj_->begin(); pIterTruth != pvTraj_->end(); pIterTruth++)
   {
      clTextFile << fixed << setprecision(6) << RadToDeg(pIterTruth->dLongitude) << ",";
      clTextFile << fixed << setprecision(6) << RadToDeg(pIterTruth->dLatitude) << endl;
   }

   clTextFile.close();
}

//------------------------------------------------------------------------------------------------
void KalmanCrossCorrPlots(void* pvErrorInfo_, string sInFilePath_, string sMethod_)
{
   string sRootDirectory;
   string sOutFile;
   ofstream clTextFile;

   // Set the root directory where all result files are going to be dumped
   SetRootDirectory(sInFilePath_, &sRootDirectory);
   sRootDirectory.append("\\").append(sMethod_);
   CreateDirectory(sRootDirectory.c_str(), NULL);

   if (sMethod_.compare("PRW") == 0)
   {
      sOutFile = sRootDirectory;
      sOutFile.append("\\PRWCrossCorrelation").append(".GIT");
      clTextFile.open(sOutFile);
      clTextFile << "[header]" << endl;
      clTextFile << "title=PRW Cross Correlation" << endl;
      clTextFile << "xlabel=Time(s)" << endl;
      clTextFile << "ylabel=Cross-Correlation" << endl;
      clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=11" << endl << "miny=-10" << endl << "maxy=10" << endl;
   }
   else if (sMethod_.compare("VRW") == 0 || sMethod_.compare("POLAR") == 0)
   {
      sOutFile = sRootDirectory;
      sOutFile.append("\\VRWCrossCorrelation").append(".GIT");
      clTextFile.open(sOutFile);
      clTextFile << "[header]" << endl;
      clTextFile << "title=VRW Cross Correlation" << endl;
      clTextFile << "xlabel=Time(s)" << endl;
      clTextFile << "ylabel=Cross-Correlation" << endl;
      clTextFile << "grid=both" << endl << "pointsize=5" << endl << "xcolumn=0" << endl << "series=29" << endl << "miny=-10" << endl << "maxy=10" << endl;
   }
   // East-North Cross-Correlation
   clTextFile << "\n[series1]" << endl;
   clTextFile << "name=East-North XCorr" << endl;
   clTextFile << "column=1" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East-Up Cross-Correlation
   clTextFile << "\n[series2]" << endl;
   clTextFile << "name=East-Up XCorr" << endl;
   clTextFile << "column=2" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East-Clock bias Cross-Correlation
   clTextFile << "\n[series3]" << endl;
   clTextFile << "name=East-ClkBias XCorr" << endl;
   clTextFile << "column=3" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // East-Clock drift Cross-Correlation
   clTextFile << "\n[series4]" << endl;
   clTextFile << "name=East-ClkDrift XCorr" << endl;
   clTextFile << "column=4" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North-Up Cross-Correlation
   clTextFile << "\n[series5]" << endl;
   clTextFile << "name=North-Up XCorr" << endl;
   clTextFile << "column=5" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North-Clock bias Cross-Correlation
   clTextFile << "\n[series6]" << endl;
   clTextFile << "name=North-ClkBias XCorr" << endl;
   clTextFile << "column=6" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // North-Clock drift Cross-Correlation
   clTextFile << "\n[series7]" << endl;
   clTextFile << "name=North-ClkDrift XCorr" << endl;
   clTextFile << "column=7" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up-Clock bias Cross-Correlation
   clTextFile << "\n[series8]" << endl;
   clTextFile << "name=Up-ClkBiast XCorr" << endl;
   clTextFile << "column=8" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Up-Clock drift Cross-Correlation
   clTextFile << "\n[series9]" << endl;
   clTextFile << "name=Up-ClkDrift XCorr" << endl;
   clTextFile << "column=9" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Clock Bias-Clock Drift Cross-Correlation
   clTextFile << "\n[series10]" << endl;
   clTextFile << "name=ClkBias-ClkDrift XCorr" << endl;
   clTextFile << "column=10" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   // Number of Satellites
   clTextFile << "\n[series11]" << endl;
   clTextFile << "name=Num SV" << endl;
   clTextFile << "column=11" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
   if (sMethod_.compare("PRW") == 0)
   {
      clTextFile << "\n[data]" << endl;
      for (vector<PRWErrorInfo>::iterator pIterPRW = ((vector<PRWErrorInfo>*)pvErrorInfo_)->begin();
           pIterPRW != ((vector<PRWErrorInfo>*)pvErrorInfo_)->end(); pIterPRW++)
      {
         clTextFile << fixed << setprecision(3) << pIterPRW->dGPSTime << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dENCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dEUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dEClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dEClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dNUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dNClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dNClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dUClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dUClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterPRW->dClkClkDriftCorr << ",";
         clTextFile << (INT)pIterPRW->cNumSV << endl;
      }
   }
   else if (sMethod_.compare("VRW") == 0 || sMethod_.compare("POLAR") == 0)
   {
      // East-Velocity east Cross-Correlation
      clTextFile << "\n[series12]" << endl;
      clTextFile << "name=East-VEast XCorr" << endl;
      clTextFile << "column=12" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // East-Velocity north Cross-Correlation
      clTextFile << "\n[series13]" << endl;
      clTextFile << "name=East-VNorth XCorr" << endl;
      clTextFile << "column=13" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // East-Velocity up Cross-Correlation
      clTextFile << "\n[series14]" << endl;
      clTextFile << "name=East-VUp XCorr" << endl;
      clTextFile << "column=14" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // North-Velocity east Cross-Correlation
      clTextFile << "\n[series15]" << endl;
      clTextFile << "name=North-VEast XCorr" << endl;
      clTextFile << "column=15" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // North-Velocity north Cross-Correlation
      clTextFile << "\n[series16]" << endl;
      clTextFile << "name=North-VNorth XCorr" << endl;
      clTextFile << "column=16" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // North-Velocity Up Cross-Correlation
      clTextFile << "\n[series17]" << endl;
      clTextFile << "name=North-VUp XCorr" << endl;
      clTextFile << "column=17" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Up-Velocity east Cross-Correlation
      clTextFile << "\n[series18]" << endl;
      clTextFile << "name=Up-VEast XCorr" << endl;
      clTextFile << "column=18" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Up-Velocity north Cross-Correlation
      clTextFile << "\n[series19]" << endl;
      clTextFile << "name=Up-VNorth XCorr" << endl;
      clTextFile << "column=19" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Up-Velocity up Cross-Correlation
      clTextFile << "\n[series20]" << endl;
      clTextFile << "name=Up-VUp XCorr" << endl;
      clTextFile << "column=20" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity east-Velocity north Cross-Correlation
      clTextFile << "\n[series21]" << endl;
      clTextFile << "name=VEast-VNorth XCorr" << endl;
      clTextFile << "column=21" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity east-Velocity up Cross-Correlation
      clTextFile << "\n[series22]" << endl;
      clTextFile << "name=VEast-VUp XCorr" << endl;
      clTextFile << "column=22" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity east-Clock bias Cross-Correlation
      clTextFile << "\n[series23]" << endl;
      clTextFile << "name=VEast-ClkBias XCorr" << endl;
      clTextFile << "column=23" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity east-Clock drift Cross-Correlation
      clTextFile << "\n[series24]" << endl;
      clTextFile << "name=VEast-ClkDrift XCorr" << endl;
      clTextFile << "column=24" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity north-Velocity up Cross-Correlation
      clTextFile << "\n[series25]" << endl;
      clTextFile << "name=VNorth-VUp XCorr" << endl;
      clTextFile << "column=25" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity north-Clock bias Cross-Correlation
      clTextFile << "\n[series26]" << endl;
      clTextFile << "name=VNorth-ClkBias XCorr" << endl;
      clTextFile << "column=26" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity north-Clock drift Cross-Correlation
      clTextFile << "\n[series27]" << endl;
      clTextFile << "name=VNorth-ClkDrift XCorr" << endl;
      clTextFile << "column=27" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity Up-Clock bias Cross-Correlation
      clTextFile << "\n[series28]" << endl;
      clTextFile << "name=VUp-ClkBias XCorr" << endl;
      clTextFile << "column=28" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      // Velocity Up-Clock drift Cross-Correlation
      clTextFile << "\n[series29]" << endl;
      clTextFile << "name=VUp-ClkDrift XCorr" << endl;
      clTextFile << "column=29" << endl << "point=dot" << endl << "connect=1" << endl << "connectdist=-1" << endl << "linewidth=1" << endl;
      clTextFile << "\n[data]" << endl;
      for (vector<VRWErrorInfo>::iterator pIterVRW = ((vector<VRWErrorInfo>*)pvErrorInfo_)->begin();
           pIterVRW != ((vector<VRWErrorInfo>*)pvErrorInfo_)->end(); pIterVRW++)
      {
         clTextFile << fixed << setprecision(3) << pIterVRW->dGPSTime << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dENCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dUClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dUClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dClkClkDriftCorr << ",";
         clTextFile << (INT)pIterVRW->cNumSV << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEVECorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEVNCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dEVUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNVECorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNVNCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dNVUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dUVECorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dUVNCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dUVUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVEVNCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVEVUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVEClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVEClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVNVUCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVNClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVNClkDriftCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVUClkCorr << ",";
         clTextFile << fixed << setprecision(3) << pIterVRW->dVUClkDriftCorr << endl;
      }
   }

   clTextFile.close();
}
