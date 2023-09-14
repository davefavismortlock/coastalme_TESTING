/*!
 *
 * \class CSimulation
 * \brief This class runs CoastalME simulations
 * \details TODO This is a more detailed description of the CSimulation class
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file simulation.h
 * \brief Contains CSimulation definitions
 *
 */

#ifndef SIMULATION_H
#define SIMULATION_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <ctime>
using std::localtime;
using std::time;
using std::time_t;

#include <fstream>
using std::ofstream;

#include <string>
using std::string;

#include <utility>
using std::pair;

#include <gdal_priv.h>

#include "line.h"
#include "i_line.h"

#include "inc/cshore.h"

int const
   NRNG = 2,
   SAVEMAX = 100000;

class CGeomRasterGrid; // Forward declarations
class CRWCoast;
class CGeomProfile;
class CGeomCoastPolygon;
class CRWCliff;
class CSedInputEvent;

class CSimulation
{
private:
   bool
      m_bHaveFineSediment,
      m_bHaveSandSediment,
      m_bHaveCoarseSediment,
      m_bRasterGISSaveAll,
      m_bVectorGISSaveAll,
      m_bBasementElevSave,
      m_bSedimentTopSurfSave,
      m_bTopSurfSave,
      m_bSliceSave,
      m_bSeaDepthSave,
      m_bAvgSeaDepthSave,
      m_bWaveHeightSave,
      m_bAvgWaveHeightSave,
      m_bWaveAngleSave,
      m_bAvgWaveAngleSave,
      m_bWaveAngleAndHeightSave,
      m_bAvgWaveAngleAndHeightSave,
      m_bDeepWaterWaveAngleAndHeightSave,
      m_bWaveEnergySinceCollapseSave,
      m_bMeanWaveEnergySave,
      m_bBreakingWaveHeightSave,
      m_bBeachProtectionSave,
      m_bPotentialPlatformErosionSave,
      m_bActualPlatformErosionSave,
      m_bTotalPotentialPlatformErosionSave,
      m_bTotalActualPlatformErosionSave,
      m_bPotentialBeachErosionSave,
      m_bActualBeachErosionSave,
      m_bTotalPotentialBeachErosionSave,
      m_bTotalActualBeachErosionSave,
      m_bBeachDepositionSave,
      m_bTotalBeachDepositionSave,
      m_bLandformSave,
      m_bLocalSlopeSave,
      m_bInterventionClassSave,
      m_bInterventionHeightSave,
      m_bSuspSedSave,
      m_bAvgSuspSedSave,
      m_bFineUnconsSedSave,
      m_bSandUnconsSedSave,
      m_bCoarseUnconsSedSave,
      m_bFineConsSedSave,
      m_bSandConsSedSave,
      m_bCoarseConsSedSave,
      m_bRasterCoastlineSave,
      m_bRasterNormalSave,
      m_bDistWeightSave,
      m_bActiveZoneSave,
      m_bCliffCollapseSave,
      m_bTotCliffCollapseSave,
      m_bCliffCollapseDepositionSave,
      m_bTotCliffCollapseDepositionSave,
      m_bRasterPolygonSave,
      m_bPotentialPlatformErosionMaskSave,
      m_bSeaMaskSave,
      m_bBeachMaskSave,
      m_bShadowZoneCodesSave,
      m_bDeepWaterWaveAngleSave,
      m_bDeepWaterWaveHeightSave,
      m_bDeepWaterWavePeriodSave,
      m_bPolygonUnconsSedUpOrDownDrift,
      m_bPolygonUnconssedGainOrLoss,
      m_bSaveRegular,
      m_bCoastSave,
      m_bNormalsSave,
      m_bInvalidNormalsSave,
      m_bCoastCurvatureSave,
      m_bPolygonNodeSave,
      m_bPolygonBoundarySave,
      m_bCliffNotchSave,
      m_bShadowBoundarySave,
      m_bShadowDowndriftBoundarySave,
      m_bSeaAreaTSSave,
      m_bStillWaterLevelTSSave,
      m_bActualPlatformErosionTSSave,
      m_bCliffCollapseErosionTSSave,
      m_bCliffCollapseDepositionTSSave,
      m_bCliffCollapseNetTSSave,
      m_bBeachErosionTSSave,
      m_bBeachDepositionTSSave,
      m_bBeachSedimentChangeNetTSSave,
      m_bSuspSedTSSave,
      m_bFloodSetupSurgeTSSave,
      m_bFloodSetupSurgeRunupTSSave,
      m_bSaveGISThisIter,
      m_bOutputProfileData,
      m_bOutputParallelProfileData,
      m_bOutputLookUpData,
      m_bOmitSearchNorthEdge,
      m_bOmitSearchSouthEdge,
      m_bOmitSearchWestEdge,
      m_bOmitSearchEastEdge,
      m_bErodeShorePlatformAlternateDirection,
      m_bDoCoastPlatformErosion,
      m_bDoCliffCollapse,
      m_bDoBeachSedimentTransport,
      m_bGDALCanCreate,
      m_bGDALCanWriteFloat,
      m_bGDALCanWriteInt32,
      m_bScaleRasterOutput,
      m_bWorldFile,
      m_bSingleDeepWaterWaveValues,
      m_bHaveWaveStationData,
      m_bSedimentInput,
      m_bSedimentInputAtPoint,
      m_bSedimentInputAtCoast,
      m_bSedimentInputAlongLine,
      m_bSedimentInputEventSave,
      m_bSedimentInputThisIter,
      m_bDoFlood,
      m_bWaveSetupSave,
      m_bStormSurgeSave,
      m_bRunUpSave,
      m_bSetupSurgeFloodMaskSave,
      m_bSetupSurgeRunupFloodMaskSave,
      m_bRasterWaveFloodLineSave,
      m_bVectorWaveFloodLineSave,
      m_bFloodLocation,
      m_bFloodSWLSetupLine,
      m_bFloodSWLSetupSurgeLine,
      m_bFloodSWLSetupSurgeRunupLine,
      m_bGISSaveDigitsSequential;

   char **m_papszGDALRasterOptions;
   char **m_papszGDALVectorOptions;

   int
      m_nXGridMax,
      m_nYGridMax,
      m_nLayers,
      m_nCoastSmooth,
      m_nCoastSmoothWindow,
      m_nSavGolCoastPoly,
      m_nProfileSmoothWindow,
      m_nCoastNormalAvgSpacing,           // In cells
      m_nCoastCurvatureInterval,          // A length, measured in coastline points
      m_nNaturalCapeNormals,
      m_nGISMaxSaveDigits,
      m_nGISSave,
      m_nUSave,
      m_nThisSave,
      m_nCoastMax,
      m_nCoastMin,
      m_nNThisIterCliffCollapse,
      m_nNTotCliffCollapse,
      m_nGlobalPolygonID,                 // There are m_nGlobalPolygonID + 1 polygons at any time (all coasts)
      m_nUnconsSedimentHandlingAtGridEdges,
      m_nBeachErosionDepositionEquation,
      m_nMissingValue,
      m_nXMinBoundingBox,
      m_nXMaxBoundingBox,
      m_nYMinBoundingBox,
      m_nYMaxBoundingBox,
      m_nWavePropagationModel,
      m_nSimStartSec,
      m_nSimStartMin,
      m_nSimStartHour,
      m_nSimStartDay,
      m_nSimStartMonth,
      m_nSimStartYear,
      m_nDeepWaterWaveDataNTimeSteps,
      m_nLogFileDetail,
      m_nRunUpEquation;

   GDALDataType
      m_GDALWriteIntDataType,
      m_GDALWriteFloatDataType;

   long
      m_lGDALMaxCanWrite,
      m_lGDALMinCanWrite;

   unsigned long
      m_ulIter,
      m_ulTotTimestep,
      m_ulRandSeed[NRNG],
      m_ulNumCells,
      m_ulThisIterNumSeaCells,
      m_ulThisIterNumCoastCells,
      m_ulThisIterNumPotentialPlatformErosionCells,
      m_ulThisIterNumActualPlatformErosionCells,
      m_ulThisIterNumPotentialBeachErosionCells,
      m_ulThisIterNumActualBeachErosionCells,
      m_ulThisIterNumBeachDepositionCells,
      m_ulTotPotentialPlatformErosionOnProfiles,
      m_ulTotPotentialPlatformErosionBetweenProfiles,
      m_ulMissingValueBasementCells;

   double
      m_dDurationUnitsMult,
      m_dNorthWestXExtCRS,
      m_dNorthWestYExtCRS,
      m_dSouthEastXExtCRS,
      m_dSouthEastYExtCRS,
      m_dExtCRSGridArea,
      m_dCellSide,                                              // Length of cell side (in external CRS units)
      m_dCellArea,                                              // Area of cell  (in external CRS units)
      m_dCellDiagonal,                                          // Length of cell's diagonal (in external CRS units)
      m_dInvCellSide,                                           // Inverse of m_dCellSide
      m_dInvCellDiagonal,                                       // Inverse of m_dCellDiagonal
      m_dSimDuration,                                           // Duration of simulation, in hours
      m_dTimeStep,
      m_dSimElapsed,                                            // Time simulated so far, in hours
      m_dRegularSaveTime,
      m_dRegularSaveInterval,
      m_dUSaveTime[SAVEMAX],
      m_dClkLast,                                               // Last value returned by clock()
      m_dCPUClock,                                              // Total elapsed CPU time
      m_dGeoTransform[6],
      m_dSeaWaterDensity,
      m_dOrigSWL,
      m_dFinalSWL,
      m_dDeltaSWLPerTimestep,
      m_dThisIterSWL,
      m_dThisIterDiffTotWaterLevel,
      m_dThisIterDiffWaveSetupWaterLevel,
      m_dThisIterDiffWaveSetupSurgeWaterLevel,
      m_dThisIterDiffWaveSetupSurgeRunupWaterLevel,
      m_dAccumulatedSeaLevelChange,
      m_dMinSWL,
      m_dMaxSWL,
      m_dBreakingWaveHeight,
      m_dC_0,                                                   // Deep water wave speed (m/s)
      m_dL_0,                                                   // Deep water wave length (m)
      m_dWaveDepthRatioForWaveCalcs,
      m_dBreakingWaveHeightDepthRatio,
      m_dAllCellsDeepWaterWaveHeight,
      m_dAllCellsDeepWaterWaveAngle,
      m_dAllCellsDeepWaterWavePeriod,
      m_dMaxUserInputWaveHeight,
      m_dMaxUserInputWavePeriod,                                // Used to constrain depth of closure
      m_dR,
      m_dD50Fine,
      m_dD50Sand,
      m_dD50Coarse,
      m_dBeachSedimentDensity,
      m_dBeachSedimentPorosity,
      m_dFineErodibility,
      m_dSandErodibility,
      m_dCoarseErodibility,
      m_dFineErodibilityNormalized,
      m_dSandErodibilityNormalized,
      m_dCoarseErodibilityNormalized,
      m_dKLS,
      m_dKamphuis,
      m_dG,
      m_dInmersedToBulkVolumetric,
      m_dDepthOfClosure,
      m_dCoastNormalAvgSpacing,                                 // In m
      m_dCoastNormalLength,
      m_dThisIterTotSeaDepth,
      m_dThisIterPotentialPlatformErosion,
      m_dThisIterActualPlatformErosionFineCons,
      m_dThisIterActualPlatformErosionSandCons,
      m_dThisIterActualPlatformErosionCoarseCons,
      m_dThisIterPotentialBeachErosion,
      m_dThisIterBeachErosionFine,
      m_dThisIterBeachErosionSand,
      m_dThisIterBeachErosionCoarse,
      m_dThisIterBeachDepositionSand,
      m_dThisIterBeachDepositionCoarse,
      m_dThisIterFineSedimentToSuspension,
      m_dThisIterPotentialSedLostBeachErosion,
      m_dThisIterLeftGridUnconsFine,
      m_dThisIterLeftGridUnconsSand,
      m_dThisIterLeftGridUnconsCoarse,
      m_dThisIterCliffCollapseFineErodedDuringDeposition,
      m_dThisIterCliffCollapseSandErodedDuringDeposition,
      m_dThisIterCliffCollapseCoarseErodedDuringDeposition,
      m_dThisIterDepositionSandDiff,
      m_dThisIterDepositionCoarseDiff,
      m_dDepthOverDBMax,                                        // Used in erosion potential look-up function
      m_dTotPotentialPlatformErosionOnProfiles,
      m_dTotPotentialPlatformErosionBetweenProfiles,
      m_dProfileMaxSlope,
      m_dMaxBeachElevAboveSWL,
      m_dCliffErosionResistance,
      m_dNotchDepthAtCollapse,
      m_dNotchBaseBelowSWL,
      m_dCliffDepositionA,
      m_dCliffDepositionPlanviewWidth,
      m_dCliffTalusMinDepositionLength,
      m_dMinCliffTalusHeightFrac,
      m_dThisIterCliffCollapseErosionFineUncons,
      m_dThisIterCliffCollapseErosionSandUncons,
      m_dThisIterCliffCollapseErosionCoarseUncons,
      m_dThisIterCliffCollapseErosionFineCons,
      m_dThisIterCliffCollapseErosionSandCons,
      m_dThisIterCliffCollapseErosionCoarseCons,
      m_dThisIterUnconsSandCliffDeposition,
      m_dThisIterUnconsCoarseCliffDeposition,
      m_dCoastNormalRandSpaceFact,
      m_dDeanProfileStartAboveSWL,
      m_dMissingValue,
      m_dWaveDataWrapHours,
      m_dThisIterTopElevMax,
      m_dThisIterTopElevMin,
      m_dThisiterUnconsFineInput,
      m_dThisiterUnconsSandInput,
      m_dThisiterUnconsCoarseInput,
      m_dStartIterSuspFine,
      m_dStartIterUnconsSand,
      m_dStartIterUnconsCoarse,
      m_dStartIterConsFine,
      m_dStartIterConsSand,
      m_dStartIterConsCoarse;

   // These grand totals are all long doubles, the aim is to minimize rounding errors when many very small numbers are added to a single much larger number, see e.g. http://www.ddj.com/cpp/184403224
   long double
      m_ldGTotPotentialPlatformErosion,
      m_ldGTotFineActualPlatformErosion,
      m_ldGTotSandActualPlatformErosion,
      m_ldGTotCoarseActualPlatformErosion,
      m_ldGTotPotentialSedLostBeachErosion,
      m_ldGTotActualFineLostBeachErosion,
      m_ldGTotActualSandLostBeachErosion,
      m_ldGTotActualCoarseLostBeachErosion,
      m_ldGTotSandSedLostCliffCollapse,
      m_ldGTotCoarseSedLostCliffCollapse,
      m_ldGTotCliffCollapseFine,
      m_ldGTotCliffCollapseSand,
      m_ldGTotCliffCollapseCoarse,
      m_ldGTotCliffTalusFineToSuspension,
      m_ldGTotCliffTalusSandDeposition,
      m_ldGTotCliffTalusCoarseDeposition,
      m_ldGTotCliffCollapseFineErodedDuringDeposition,
      m_ldGTotCliffCollapseSandErodedDuringDeposition,
      m_ldGTotCliffCollapseCoarseErodedDuringDeposition,
      m_ldGTotPotentialBeachErosion,
      m_ldGTotActualFineBeachErosion,
      m_ldGTotActualSandBeachErosion,
      m_ldGTotActualCoarseBeachErosion,
      m_ldGTotSandBeachDeposition,
      m_ldGTotCoarseBeachDeposition,
      m_ldGTotSuspendedSediment,
      m_ldGTotSandDepositionDiff,
      m_ldGTotCoarseDepositionDiff,
      m_ldGTotFineSedimentInput,
      m_ldGTotSandSedimentInput,
      m_ldGTotCoarseSedimentInput;

   string
      m_strCMEDir,
      m_strCMEIni,
      m_strMailAddress,
      m_strDataPathName,
      m_strRasterGISOutFormat,
      m_strVectorGISOutFormat,
      m_strInitialBasementDEMFile,
      m_strInitialLandformFile,
      m_strInterventionClassFile,
      m_strInterventionHeightFile,
      m_strInitialSuspSedimentFile,
      m_strShapeFunctionFile,
      m_strTideDataFile,
      m_strLogFile,
      m_strOutPath,
      m_strOutFile,
      m_strPalFile,
      m_strGDALBasementDEMDriverCode, // Basement DEM (raster)
      m_strGDALBasementDEMDriverDesc,
      m_strGDALBasementDEMProjection,
      m_strGDALBasementDEMDataType,
      m_strGDALLDriverCode, // Initial landform class (raster)
      m_strGDALLDriverDesc,
      m_strGDALLProjection,
      m_strGDALLDataType,
      m_strGDALICDriverCode, // Initial intervention class (raster)
      m_strGDALICDriverDesc,
      m_strGDALICProjection,
      m_strGDALICDataType,
      m_strGDALIHDriverCode, // Initial intervention class (raster)
      m_strGDALIHDriverDesc,
      m_strGDALIHProjection,
      m_strGDALIHDataType,
      m_strGDALIWDriverCode, // Initial water depth (raster)
      m_strGDALIWDriverDesc,
      m_strGDALIWProjection,
      m_strGDALIWDataType,
      m_strGDALISSDriverCode, // Initial suspended sediment (raster)
      m_strGDALISSDriverDesc,
      m_strGDALISSProjection,
      m_strGDALISSDataType,
      m_strOGRDWWVDriverCode, // Initial deep water wave stations (vector)
      m_strOGRDWWVGeometry,
      m_strOGRDWWVDataType,
      m_strOGRSedInputDriverCode, // Sediment input event locations (vector)
      m_strOGRSedInputGeometry,
      m_strOGRSedInputDataType,
      m_strOGRFloodDriverCode, // Flood input locations (point or vector)
      m_strOGRFloodGeometry,
      m_strOGRFloodDataType,
      m_strGDALRasterOutputDriverLongname,
      m_strGDALRasterOutputDriverExtension,
      m_strOGRVectorOutputExtension,
      m_strRunName,
      m_strDurationUnits,
      m_strDeepWaterWaveStationsShapefile,
      m_strDeepWaterWavesTimeSeriesFile,
      m_strSedimentInputEventShapefile,
      m_strSedimentInputEventTimeSeriesFile,
      // m_strLevel,
      m_strFloodLocationShapefile;

   int m_nLevel;

   struct RandState
   {
      unsigned long s1, s2, s3;
   } m_ulRState[NRNG];

   time_t
      m_tSysStartTime,
      m_tSysEndTime;

   ofstream
      OutStream,
      SeaAreaTSStream,
      StillWaterLevelTSStream,
      ErosionTSStream,
      CliffCollapseErosionTSStream,
      CliffCollapseDepositionTSStream,
      CliffCollapseNetTSStream,
      BeachErosionTSStream,
      BeachDepositionTSStream,
      BeachSedimentChangeNetTSStream,
      SedLoadTSStream,
      FloodSetupSurgeTSStream,
      FloodSetupSurgeRunupTSStream;

   vector<bool>
      m_bConsChangedThisIter,
      m_bUnconsChangedThisIter;

   vector<int>
      m_VnProfileToSave,
      m_VnDeepWaterWaveStationID,  // ID for deep water wave station, this corresponds with the ID in the wave time series file
      m_VnSedimentInputLocationID, // ID for sediment input location, this corresponds with the ID in the sediment input time series file
      m_VnFloodLocationID,         // ID for flood location
      m_VnSavGolIndexCoast;        // Savitzky-Golay shift index for the coastline vector(s)

   vector<unsigned long>
      m_VulProfileTimestep,
      m_VlDeepWaterWaveValuesAtTimestep; // Calculate deep water wave values at these timesteps

   vector<double>
      m_VdSliceElev,
      m_VdErosionPotential,                   // For erosion potential lookup
      m_VdDepthOverDB,                        // For erosion potential lookup
      m_VdSavGolFCRWCoast,                    // Savitzky-Golay filter coefficients for the coastline vector(s)
      m_VdSavGolFCGeomProfile,                // Savitzky-Golay filter coefficients for the profile vectors
      m_VdTideData,                           // Tide data: one record per timestep, is the change (m) from still water level for that timestep
      m_VdDeepWaterWaveStationX,              // X co-ordinate (grid CRS) for deep water wave station
      m_VdDeepWaterWaveStationY,              // Y co-ordinate (grid CRS) for deep water wave station
      m_VdThisIterDeepWaterWaveStationHeight, // This iteration wave height at deep water wave station
      m_VdThisIterDeepWaterWaveStationAngle,  // This iteration wave orientation at deep water wave station
      m_VdThisIterDeepWaterWaveStationPeriod, // This iteration wave period at deep water wave station
      m_VdTSDeepWaterWaveStationHeight,       // Time series of wave heights at deep water wave station
      m_VdTSDeepWaterWaveStationAngle,        // Time series of wave orientation at deep water wave station
      m_VdTSDeepWaterWaveStationPeriod,       // Time series of wave period at deep water wave station
      m_VdSedimentInputLocationX,             // X co-ordinate (grid CRS) for sediment input event
      m_VdSedimentInputLocationY,             // X co-ordinate (grid CRS) for sediment input event
      m_VdFloodLocationX,                     // X co-ordinate (grid CRS) for total water level flooding
      m_VdFloodLocationY;                     // X co-ordinate (grid CRS) for total water level flooding

   vector<string>
      m_VstrInitialFineUnconsSedimentFile,
      m_VstrInitialSandUnconsSedimentFile,
      m_VstrInitialCoarseUnconsSedimentFile,
      m_VstrInitialFineConsSedimentFile,
      m_VstrInitialSandConsSedimentFile,
      m_VstrInitialCoarseConsSedimentFile,
      m_VstrGDALIUFDriverCode,
      m_VstrGDALIUFDriverDesc,
      m_VstrGDALIUFProjection,
      m_VstrGDALIUFDataType,
      m_VstrGDALIUSDriverCode,
      m_VstrGDALIUSDriverDesc,
      m_VstrGDALIUSProjection,
      m_VstrGDALIUSDataType,
      m_VstrGDALIUCDriverCode,
      m_VstrGDALIUCDriverDesc,
      m_VstrGDALIUCProjection,
      m_VstrGDALIUCDataType,
      m_VstrGDALICFDriverCode,
      m_VstrGDALICFDriverDesc,
      m_VstrGDALICFProjection,
      m_VstrGDALICFDataType,
      m_VstrGDALICSDriverCode,
      m_VstrGDALICSDriverDesc,
      m_VstrGDALICSProjection,
      m_VstrGDALICSDataType,
      m_VstrGDALICCDriverCode,
      m_VstrGDALICCDriverDesc,
      m_VstrGDALICCProjection,
      m_VstrGDALICCDataType;

   // Pointer to the raster grid object
   CGeomRasterGrid *m_pRasterGrid;

   // The coastline objects
   vector<CRWCoast> m_VCoast;
   // vector<CRWCoast> m_VFloodWaveSetup;
   vector<CRWCoast> m_VFloodWaveSetupSurge;
   vector<CRWCoast> m_VFloodWaveSetupSurgeRunup;

   // Pointers to coast polygon objects
   vector<CGeomCoastPolygon *> m_pVCoastPolygon;

   // Edge cells
   vector<CGeom2DIPoint> m_VEdgeCell;

   // And the grid edge that each edge cell belongs to
   vector<int> m_VEdgeCellEdge;

   // And the location to compute the total water level for flooding
   vector<int> m_VCellFloodLocation;

   // Sediment input events
   vector<CSedInputEvent *> m_pVSedInputEvent;

private:
   // Input and output routines
   static int nHandleCommandLineParams(int, char const* []);
   bool bReadIniFile(void);
   bool bReadRunDataFile(void);
   bool bOpenLogFile(void);
   bool bSetUpTSFiles(void);
   void WriteStartRunDetails(void);
   bool bWritePerTimestepResults(void);
   bool bWriteTSFiles(void);
   int nWriteEndRunDetails(void);
   int nReadShapeFunctionFile(void);
   int nReadWaveStationTimeSeriesFile(int const);
   int nReadSedimentInputEventTimeSeriesFile(void);
   int nReadTideDataFile(void);
   int nSaveProfile(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*);
   bool bWriteProfileData(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*) const;
   int nSaveParProfile(int const, int const, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*);
   bool bWriteParProfileData(int const, int const, int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>*const, vector<double> const*) const;
   void WriteLookUpData(void);

   // GIS input and output stuff
   int nReadRasterBasementDEM(void);
   int nReadRasterGISFile(int const, int const);
   int nReadVectorGISFile(int const);
   bool bWriteRasterGISFile(int const, string const*, int const = 0, double const = 0);
   bool bWriteVectorGISFile(int const, string const*);
   void GetRasterOutputMinMax(int const, double&, double&, int const, double const);
   void SetRasterFileCreationDefaults(void);
   int nInterpolateWavesToPolygonCells(vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*);

   // Initialization
   bool bCreateErosionPotentialLookUp(vector<double>*, vector<double>*, vector<double>*);

   // Top-level simulation routines
   static int nUpdateIntervention(void);
   int nCheckForSedimentInputEvent(void);
   int nCalcExternalForcing(void);
   int nInitGridAndCalcStillWaterLevel(void);
   int nLocateSeaAndCoasts(void);
   int nLocateFloodAndCoasts(void);
   static int nLocateAllEstuaries(void);
   int nAssignAllCoastalLandforms(void);
   int nAssignNonCoastlineLandforms(void);
   int nDoAllPropagateWaves(void);
   int nDoAllShorePlatFormErosion(void);
   int nDoAllWaveEnergyToCoastLandforms(void);
   int nDoCliffCollapse(int const, CRWCliff*, double&, double&, double&, double&, double&);
   int nDoCliffCollapseDeposition(int const, CRWCliff*, double const, double const, double const, double const);
   int nUpdateGrid(void);

   // Lower-level simulation routines
   void FindAllSeaCells(void);
   int FindAllInundatedCells(void);
   void FloodFillSea(int const, int const);
   void FloodFillLand(int const, int const); // Recursive Scanline 8-Way Floodfill Algorithm (floodFillScanline)
   int nTraceCoastLine(unsigned int const, int const, int const, vector<bool>*, vector<CGeom2DIPoint> const*);
   int nTraceAllCoasts(void);
   int nTraceFloodCoastLine(unsigned int const, int const, int const, vector<bool>*, vector<CGeom2DIPoint> const*);
   int nTraceAllFloodCoasts(void);
   void DoCoastCurvature(int const, int const);
   int nCreateAllProfilesAndCheckForIntersection(void);
   int nCreateAllProfiles(void);
   void CreateNaturalCapeNormalProfiles(int const, int&, int const, vector<bool>*, vector<pair<int, double>> const*);
   void CreateRestOfNormalProfiles(int const, int&, int const, double const, vector<bool>*, vector<pair<int, double>> const*);
   void CreateInterventionProfiles(int const, int& /*, int const*/);
   int nCreateProfile(int const, int const, int&);
   int nCreateGridEdgeProfile(bool const, int const, int&);
   int nPutAllProfilesOntoGrid(void);
   int nModifyAllIntersectingProfiles(void);
   static bool bCheckForIntersection(CGeomProfile *const, CGeomProfile *const, int&, int&, double&, double&, double&, double&);
   void MergeProfilesAtFinalLineSegments(int const, int const, int const, int const, int const, double const, double const, double const, double const);
   void TruncateOneProfileRetainOtherProfile(int const, int const, int const, double const, double const, int const, int const, bool const);
   int nInsertPointIntoProfilesIfNeededThenUpdate(int const, int const, double const, double const, int const, int const, int const, bool const);
   void TruncateProfileAndAppendNew(int const, int const, int const, vector<CGeom2DPoint> const*, vector<vector<pair<int, int>>> const*);
   void RasterizeProfile(int const, int const, vector<CGeom2DIPoint>*, vector<bool>*, bool&, bool&, bool&, bool&, bool&);
   static void CalcDeanProfile(vector<double>*, double const, double const, double const, bool const, int const, double const);
   static double dSubtractProfiles(vector<double> const*, vector<double> const*, vector<bool> const*);
   int nRasterizeCliffCollapseProfile(vector<CGeom2DPoint> const*, vector<CGeom2DIPoint>*) const;
   int nCalcPotentialPlatformErosionOnProfile(int const, int const);
   int nCalcPotentialPlatformErosionBetweenProfiles(int const, int const, int const);
   void ConstructParallelProfile(int const, int const, int const, int const, int const, vector<CGeom2DIPoint>*const, vector<CGeom2DIPoint>*, vector<CGeom2DPoint> *);
   double dCalcBeachProtectionFactor(int const, int const, double const);
   void FillInBeachProtectionHoles(void);
   void FillPotentialPlatformErosionHoles(void);
   void DoActualPlatformErosionOnCell(int const, int const);
   double dLookUpErosionPotential(double const) const;
   static CGeom2DPoint PtChooseEndPoint(int const, CGeom2DPoint const*, CGeom2DPoint const*, double const, double const, double const, double const);
   int nGetCoastNormalEndPoint(int const, int const, int const, CGeom2DPoint const*, double const, CGeom2DPoint *, CGeom2DIPoint*);
   int nLandformToGrid(int const, int const);
   int nCalcWavePropertiesOnProfile(int const, int const, int const, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<bool>*);
   int nGetThisProfileElevationVectorsForCShore(int const, int const, int const, vector<double>*, vector<double>*, vector<double>*);
   int nCreateCShoreInfile(int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, double const, double const, double const, double const, double const, double const, double const, double const, vector<double> const*, vector<double> const*, vector<double> const*);
   int nReadCShoreOutput(int const, string const*, int const, int const, vector<double> const*, vector<double>*);   
   static void InterpolateCShoreOutput(vector<double> const*, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double>*, vector<double>*, vector<double>*, vector<double>*);
   static double dCalcWaveAngleToCoastNormal(double const, double const, int const);
   void CalcCoastTangents(int const);
   void InterpolateWavePropsBetweenProfiles(int const, int const, int const);
   void InterpolateWaveHeightToCoastPoints(int const);
   void InterpolateWavePropertiesToCells(int const, int const, int const);
   void ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const, int const);
   static double dCalcCurvature(int const, CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void CalcD50AndFillWaveCalcHoles(void);
   int nDoAllShadowZones(void);
   static bool bOnOrOffShoreAndUpOrDownCoast(double const, double const, int const, bool&);
   static CGeom2DIPoint PtiFollowWaveAngle(CGeom2DIPoint const*, double const, double&);
   int nFindAllShadowZones(void);
   int nFloodFillShadowZone(int const, CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);
   int nDoShadowZoneAndDownDriftZone(int const, int const, int const, int const);
   void ProcessDownDriftCell(int const, int const, int const, double const, int const);
   void ProcessShadowZoneCell(int const, int const, int const, CGeom2DIPoint const*, int const, int const, int const);
   int nCreateAllPolygons(void);
   void RasterizePolygonJoiningLine(CGeom2DPoint const*, CGeom2DPoint const*);
   static bool bIsWithinPolygon(CGeom2DPoint const*, vector<CGeom2DPoint> const*);
   static CGeom2DPoint PtFindPointInPolygon(vector<CGeom2DPoint> const*, int const);
   void MarkPolygonCells(void);
   int nDoPolygonSharedBoundaries(void);
   void DoAllPotentialBeachErosion(void);
   int nDoAllActualBeachErosionAndDeposition(void);
   // int nEstimateBeachErosionOnPolygon(int const, int const, double const);
   int nEstimateErosionOnPolygon(int const, int const, double const, double&, double&, double&);
   int nEstimateUnconsErosionOnParallelProfile(/*int const, int const,*/ int const, int const, /* int const, */ int const, vector<CGeom2DIPoint> const*, vector<double> const*, double&, double&, double&, double&, double&);
   int nDoParallelProfileUnconsErosion( int const, int const,  int const, int const, int const,  int const,  int const, vector<CGeom2DIPoint> const*, vector<double> const*, double&, double&, double&);
   // void EstimateUnconsErosionOnCell(int const, int const, int const, double const, double&, double&, double&);
   void ErodeCellBeachSedimentSupplyLimited(int const, int const, int const, int const, double const, double&);
   int nEstimateMovementUnconsToAdjacentPolygons(int const, int const);
   int nDoUnconsErosionOnPolygon(int const, int const, int const, double const, double&);
   int nDoUnconsDepositionOnPolygon(int const, int const, int const, double, double&);
   void CalcDepthOfClosure(void);
   int nInterpolateAllDeepWaterWaveValues(void);
   int nSetAllCoastpointDeepWaterWaveValues(void);
   int nDoSedimentInputEvent(int const);
   void AllPolygonsUpdateStoredUncons(int const);

   // GIS utility routines
   int nMarkBoundingBoxEdgeCells(void);
   bool bCheckRasterGISOutputFormat(void);
   bool bCheckVectorGISOutputFormat(void);
   bool bSaveAllRasterGISFiles(void);
   bool bSaveAllVectorGISFiles(void);
   bool bIsWithinValidGrid(int const, int const) const;
   bool bIsWithinValidGrid(CGeom2DIPoint const*) const;
   double dGridCentroidXToExtCRSX(int const) const;
   double dGridCentroidYToExtCRSY(int const) const;
   double dGridXToExtCRSX(double const) const;
   double dGridYToExtCRSY(double const) const;
   double dExtCRSXToGridCentroidX(double const) const;
   double dExtCRSYToGridCentroidY(double const) const;
   CGeom2DIPoint PtiExtCRSToGrid(CGeom2DPoint const*) const;
   CGeom2DPoint PtGridCentroidToExt(CGeom2DIPoint const*) const;
   double dExtCRSXToGridX(double const) const;
   double dExtCRSYToGridY(double const) const;
   static double dGetDistanceBetween(CGeom2DPoint const*, CGeom2DPoint const*);
   static double dGetDistanceBetween(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static double dTriangleAreax2(CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void KeepWithinValidGrid(int, int, int&, int&) const;
   void KeepWithinValidGrid(CGeom2DIPoint const*, CGeom2DIPoint*) const;
   static double dKeepWithin360(double const);
   //    vector<CGeom2DPoint> VGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DPoint PtGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(CGeom2DIPoint const*, CGeom2DIPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(int const, int const, int const, int const, double const, int const);
   static CGeom2DPoint PtAverage(CGeom2DPoint const*, CGeom2DPoint const*);
   static CGeom2DPoint PtAverage(vector<CGeom2DPoint> *);
   static CGeom2DIPoint PtiAverage(CGeom2DIPoint const*, CGeom2DIPoint const*);
   static CGeom2DIPoint PtiAverage(vector<CGeom2DIPoint>*);
   static CGeom2DIPoint PtiWeightedAverage(CGeom2DIPoint const*, CGeom2DIPoint const*, double const);
   static CGeom2DIPoint PtiPolygonCentroid(vector<CGeom2DIPoint>*);
   static double dAngleSubtended(CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);
   static int nGetOppositeDirection(int const);
   static void GetSlopeAndInterceptFromPoints(CGeom2DIPoint const*, CGeom2DIPoint const*, double&, double&);
   CGeom2DIPoint PtiFindClosestCoastPoint(int const, int const);
   int nConvertMetresToNumCells(double const);

   // Utility routines
   static void AnnounceStart(void);
   void AnnounceLicence(void);
   void AnnounceReadBasementDEM(void) const;
   static void AnnounceAddLayers(void);
   static void AnnounceReadRasterFiles(void);
   static void AnnounceReadVectorFiles(void);
   void AnnounceReadLGIS(void) const;
   void AnnounceReadICGIS(void) const;
   void AnnounceReadIHGIS(void) const;
   static void AnnounceInitializing(void);
   void AnnounceReadInitialSuspSedGIS(void) const;
   void AnnounceReadInitialFineUnconsSedGIS(int const) const;
   void AnnounceReadInitialSandUnconsSedGIS(int const) const;
   void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
   void AnnounceReadInitialFineConsSedGIS(int const) const;
   void AnnounceReadInitialSandConsSedGIS(int const) const;
   void AnnounceReadInitialCoarseConsSedGIS(int const) const;
   void AnnounceReadDeepWaterWaveValuesGIS(void) const;
   void AnnounceReadSedimentEventInputValuesGIS(void) const;
   void AnnounceReadFloodLocationGIS(void) const;
   void AnnounceReadTideData(void) const;
   static void AnnounceReadSCAPEShapeFunctionFile(void);
   static void AnnounceAllocateMemory(void);
   static void AnnounceIsRunning(void);
   static void AnnounceSimEnd(void);
   void StartClock(void);
   bool bFindExeDir(char const*);
   bool bTimeToQuit(void);
   static int nDoTimeUnits(string const*);
   int nDoSimulationTimeMultiplier(string const*);
   static double dGetTimeMultiplier(string const*);
   static bool bParseDate(string const*, int&, int&, int&);
   static bool bParseTime(string const*, int&, int&, int&);
   void DoTimestepTotals(void);
   static string strGetBuild(void);
   static string strGetComputerName(void);
   void DoCPUClockReset(void);
   void CalcTime(double const);
   static string strDispTime(double const, bool const, bool const);
   static string strDispSimTime(double const);
   void AnnounceProgress(void);
   static string strGetErrorText(int const);
   string strListRasterFiles(void) const;
   string strListVectorFiles(void) const;
   string strListTSFiles(void) const;
   void CalcProcessStats(void);
   void CalcSavitzkyGolayCoeffs(void);
   CGeomLine LSmoothCoastSavitzkyGolay(CGeomLine *, int const, int const) const;
   CGeomLine LSmoothCoastRunningMean(CGeomLine *) const;
   vector<double> dVSmoothProfileSlope(vector<double>*);
   //    vector<double> dVCalCGeomProfileSlope(vector<CGeom2DPoint>*, vector<double>*);
   vector<double> dVSmoothProfileSavitzkyGolay(vector<double>*, vector<double>*);
   vector<double> dVSmoothProfileRunningMean(vector<double>*);
   static void CalcSavitzkyGolay(double[], int const, int const, int const, int const, int const);
   static string pstrChangeToBackslash(string const*);
   static string pstrChangeToForwardSlash(string const*);
   static string strTrim(string const*);
   static string strTrimLeft(string const*);
   static string strTrimRight(string const*);
   static string strToLower(string const*);
   //  static string strToUpper(string const*);
   static string strRemoveSubstr(string *, string const*);
   static vector<string> *VstrSplit(string const*, char const, vector<string> *);
   static vector<string> VstrSplit(string const*, char const);
   static double dCrossProduct(double const, double const, double const, double const, double const, double const);
   static double dGetMean(vector<double> const*);
   static double dGetStdDev(vector<double> const*);
   static void AppendEnsureNoGap(vector<CGeom2DIPoint>*, CGeom2DIPoint const*);
   static bool bIsNumeric(string const*);
   unsigned long ulConvertToTimestep(string const*);
   void WritePolygonShareTable(int const);
   void WritePolygonPreExistingSediment(int const);
   void WritePolygonShorePlatformErosion(int const);
   void WritePolygonCliffCollapseErosion(int const);
   void WritePolygonSedimentBeforeMovement(int const);
   void WritePolygonPotentialErosion(int const);
   void WritePolygonUnconsErosion(int const);
   void WritePolygonUnsortedSequence(int const, vector<vector<int> >&);
   void WritePolygonSortedSequence(int const, vector<vector<int> >&);
   void WritePolygonEstimatedMovement(int const, vector<vector<int> >&);
   void WritePolygonActualMovement(int const, vector<vector<int> >&);

   // Random number stuff
   static unsigned long ulGetTausworthe(unsigned long const, unsigned long const, unsigned long const, unsigned long const, unsigned long const);
   void InitRand0(unsigned long const);
   void InitRand1(unsigned long const);
   unsigned long ulGetRand0(void);
   unsigned long ulGetRand1(void);
   static unsigned long ulGetLCG(unsigned long const); // Used by all generators
   double dGetRand0d1(void);
   //    int nGetRand0To(int const);
   int nGetRand1To(int const);
   //    double dGetRand0GaussPos(double const, double const);
   double dGetRand0Gaussian(void);
   //    double dGetCGaussianPDF(double const);
   void Rand1Shuffle(int *, int);
#ifdef RANDCHECK
   void CheckRand(void) const;
#endif

public:
   ofstream LogStream;

   CSimulation(void);
   ~CSimulation(void);

   //! Returns the NODATA value
   double dGetMissingValue(void) const;

   //! Returns this timestep's still water level
   double dGetThisIterSWL(void) const;

   //! Returns this timestep's total water level
   double dGetThisIterTotWaterLevel(void) const;

   //! Returns the vertical tolerance for beach cells to be included in smoothing
   double dGetMaxBeachElevAboveSWL(void) const;

   //! Returns the cell size
   //    double dGetCellSide(void) const;

   //! Returns the size of the grid in the X direction
   int nGetGridXMax(void) const;

   //! Returns the size of the grid in the Y direction
   int nGetGridYMax(void) const;

   //! Returns the global d50 value for fine sediment
   double dGetD50Fine(void) const;

   //! Returns the global d50 value for sand sediment
   double dGetD50Sand(void) const;

   //! Returns the global d50 value for coarse sediment
   double dGetD50Coarse(void) const;

   //! Runs the simulation
   int nDoSimulation(int, char const* []);

   //! Carries out end-of-simulation tidying (error messages etc.)
   void DoSimulationEnd(int const);
};
#endif // SIMULATION_H
