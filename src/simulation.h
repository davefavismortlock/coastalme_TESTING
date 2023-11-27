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

#include <stack>
using std::stack;

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
   //! Does this simulation consider fine-sized sediment?
   bool m_bHaveFineSediment;

   //! Does this simulation consider sand-sized sediment?
   bool m_bHaveSandSediment;

   //! Does this simulation consider coarse-sized sediment?
   bool m_bHaveCoarseSediment;

   //! Save all raster GIS files?
   bool m_bRasterGISSaveAll;

   //! Save all vector GIS files?
   bool m_bVectorGISSaveAll;

   //! Save basement raster DEMs?
   bool m_bBasementElevSave;

   //! Save sediment top surface raster DEMs?
   bool m_bSedimentTopSurfSave;

   //! Save fop surface (sediment and sea) raster DEMs?
   bool m_bTopSurfSave;

   //! Save slices?
   bool m_bSliceSave;

   //! Save sea depth raster GIS files?
   bool m_bSeaDepthSave;

   //! Save average sea depth raster GIS files?
   bool m_bAvgSeaDepthSave;

   //! Save wave height raster GIS files?
   bool m_bWaveHeightSave;

   //! Save wave height raster GIS files?
   bool m_bAvgWaveHeightSave;

   //! Save wave angle raster GIS files?
   bool m_bWaveAngleSave;

   //! Save average wave angle raster GIS files?
   bool m_bAvgWaveAngleSave;

   //! Save wave angle and wave height raster GIS files?
   bool m_bWaveAngleAndHeightSave;

   //! Save average wave angle and average wave height raster GIS files?
   bool m_bAvgWaveAngleAndHeightSave;

   //! Save deep water wave angle and wave height raster GIS files?
   bool m_bDeepWaterWaveAngleAndHeightSave;

   //! Save wave energy since cliff collapse raster GIS files?
   bool m_bWaveEnergySinceCollapseSave;

   //! Save mean wave energy raster GIS files?
   bool m_bMeanWaveEnergySave;

   //! Save breaking wave height raster GIS files?
   bool m_bBreakingWaveHeightSave;

   //! Save beach protection raster GIS files>
   bool m_bBeachProtectionSave;

   //! Save potential shore platform erosion raster GIS files?
   bool m_bPotentialPlatformErosionSave;

   //! Save actual (supply-limited) shore platform erosion raster GIS files?
   bool m_bActualPlatformErosionSave;

   //! Save total potential shore platform erosion raster GIS files?
   bool m_bTotalPotentialPlatformErosionSave;

   //! Save total actual (supply-limited) shore platform erosion raster GIS files?
   bool m_bTotalActualPlatformErosionSave;

   //! Save potential beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bPotentialBeachErosionSave;

   //! Save actual (supply-limited) beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bActualBeachErosionSave;

   //! Save total potential beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bTotalPotentialBeachErosionSave;

   //! Save total actual (supply-limited) beach (unconsolidated sediment) erosion raster GIS files?
   bool m_bTotalActualBeachErosionSave;

   //! Save beach (unconsolidated sediment) deposition raster GIS files?
   bool m_bBeachDepositionSave;

   //! Save total beach (unconsolidated sediment) deposition raster GIS files?
   bool m_bTotalBeachDepositionSave;

   //! Save coast landform raster GIS files?
   bool m_bLandformSave;

   //! Save local slope raster GIS files?
   bool m_bLocalSlopeSave;

   //! Save intervention class raster GIS files?
   bool m_bInterventionClassSave;

   //! Save intervention height raster GIS files?
   bool m_bInterventionHeightSave;

   //! Save suspended sediment raster GIS files?
   bool m_bSuspSedSave;

   //! Save average suspended sediment raster GIS files?
   bool m_bAvgSuspSedSave;

   //! Save fine unconsolidated sediment raster GIS files?
   bool m_bFineUnconsSedSave;

   //! Save sand unconsolidated sediment raster GIS files?
   bool m_bSandUnconsSedSave;

   //! Save coarse unconsolidated sediment raster GIS files?
   bool m_bCoarseUnconsSedSave;

   //! Save fine consolidated sediment raster GIS files?
   bool m_bFineConsSedSave;

   //! Save sand consolidated sediment raster GIS files?
   bool m_bSandConsSedSave;

   //! Save coarse consolidated sediment raster GIS files?
   bool m_bCoarseConsSedSave;

   //! Save raster coastline GIS files?
   bool m_bRasterCoastlineSave;

   //! Save raster coastline-normal GIS files?
   bool m_bRasterNormalSave;

   //! Save active zone raster GIS files?
   bool m_bActiveZoneSave;

   //! Save cliff collapse raster GIS files?
   bool m_bCliffCollapseSave;

   //! Save total cliff collapse raster GIS files?
   bool m_bTotCliffCollapseSave;

   //! Save cliff collapse deposition raster GIS files?
   bool m_bCliffCollapseDepositionSave;

   //! Save total cliff collapse deposition raster GIS files?
   bool m_bTotCliffCollapseDepositionSave;

   //! Save raster polygon raster GIS files?
   bool m_bRasterPolygonSave;

   //! Save potential platform erosion mask raster GIS files?
   bool m_bPotentialPlatformErosionMaskSave;

   //! Save sea mask raster GIS files?
   bool m_bSeaMaskSave;

   //! Save beach mask raster GIS files?
   bool m_bBeachMaskSave;

   //! Save wave shadow zones raster GIS files?
   bool m_bShadowZoneCodesSave;

   //! Save deep water wave angle raster GIS files?
   bool m_bDeepWaterWaveAngleSave;

   //! Save deep water wave height raster GIS files?
   bool m_bDeepWaterWaveHeightSave;

   //! Save deep water wave period raster GIS files?
   bool m_bDeepWaterWavePeriodSave;

   //! Save polygon unconsolidated sediment up- or down-drift raster GIS files?
   bool m_bPolygonUnconsSedUpOrDownDriftSave;

   //! Save polygon unconsolidated sediment gain or loss raster GIS files?
   bool m_bPolygonUnconsSedGainOrLossSave;

   //! Save GIS files at regular intervals?
   bool m_bSaveRegular;

   //! Save
   bool m_bCoastSave;

   //! Save coastline-normal vector GIS files?
   bool m_bNormalsSave;

   //! Save invalid coastline-normal vector GIS files?
   bool m_bInvalidNormalsSave;

   //! Save coastline-curvature vector GIS files?
   bool m_bCoastCurvatureSave;

   //! Save polygon node vector GIS files?
   bool m_bPolygonNodeSave;

   //! Save polygon boundary vector GIS files?
   bool m_bPolygonBoundarySave;

   //! Save cliff notch incision depth vector GIS files?
   bool m_bCliffNotchSave;

   //! Save wave shadow boundary vector GIS files?
   bool m_bShadowBoundarySave;

   //! Save wave shadow downdrift boundary vector GIS files?
   bool m_bShadowDowndriftBoundarySave;

   bool m_bSeaAreaTSSave;
   bool m_bStillWaterLevelTSSave;
   bool m_bActualPlatformErosionTSSave;
   bool m_bCliffCollapseErosionTSSave;
   bool m_bCliffCollapseDepositionTSSave;
   bool m_bCliffCollapseNetTSSave;
   bool m_bBeachErosionTSSave;
   bool m_bBeachDepositionTSSave;
   bool m_bBeachSedimentChangeNetTSSave;
   bool m_bSuspSedTSSave;
   bool m_bFloodSetupSurgeTSSave;
   bool m_bFloodSetupSurgeRunupTSSave;
   bool m_bSaveGISThisIter;
   bool m_bOutputProfileData;
   bool m_bOutputParallelProfileData;
   bool m_bOutputLookUpData;
   bool m_bOmitSearchNorthEdge;
   bool m_bOmitSearchSouthEdge;
   bool m_bOmitSearchWestEdge;
   bool m_bOmitSearchEastEdge;
   bool m_bErodeShorePlatformAlternateDirection;
   bool m_bDoCoastPlatformErosion;
   bool m_bDoCliffCollapse;
   bool m_bDoBeachSedimentTransport;
   bool m_bGDALCanCreate;
   bool m_bGDALCanWriteFloat;
   bool m_bGDALCanWriteInt32;
   bool m_bScaleRasterOutput;
   bool m_bWorldFile;
   bool m_bSingleDeepWaterWaveValues;
   bool m_bHaveWaveStationData;
   bool m_bSedimentInput;
   bool m_bSedimentInputAtPoint;
   bool m_bSedimentInputAtCoast;
   bool m_bSedimentInputAlongLine;
   bool m_bSedimentInputEventSave;
   bool m_bSedimentInputThisIter;
   bool m_bDoFlood;
   bool m_bWaveSetupSave;
   bool m_bStormSurgeSave;
   bool m_bRunUpSave;
   bool m_bSetupSurgeFloodMaskSave;
   bool m_bSetupSurgeRunupFloodMaskSave;
   bool m_bRasterWaveFloodLineSave;
   bool m_bVectorWaveFloodLineSave;
   bool m_bFloodLocation;
   bool m_bFloodSWLSetupLine;
   bool m_bFloodSWLSetupSurgeLine;
   bool m_bFloodSWLSetupSurgeRunupLine;
   bool m_bGISSaveDigitsSequential;

   char** m_papszGDALRasterOptions;
   char** m_papszGDALVectorOptions;

   int m_nXGridMax;
   int m_nYGridMax;
   int m_nLayers;
   int m_nCoastSmooth;
   int m_nCoastSmoothWindow;
   int m_nSavGolCoastPoly;
   int m_nProfileSmoothWindow;

   //! Average spacing between cost normals, measured in cells
   int m_nCoastNormalAvgSpacing;

   //! Coast curvature interval is a length, measured in coastline points
   int m_nCoastCurvatureInterval;
   int m_nNaturalCapeNormals;
   int m_nGISMaxSaveDigits;
   int m_nGISSave;
   int m_nUSave;
   int m_nThisSave;
   int m_nCoastMax;
   int m_nCoastMin;
   int m_nNThisIterCliffCollapse;
   int m_nNTotCliffCollapse;

   //! Global (all coasts) polygon ID. There are m_nGlobalPolygonID + 1 polygons at any time
   int m_nGlobalPolygonID;
   int m_nUnconsSedimentHandlingAtGridEdges;
   int m_nBeachErosionDepositionEquation;
   int m_nMissingValue;
   int m_nXMinBoundingBox;
   int m_nXMaxBoundingBox;
   int m_nYMinBoundingBox;
   int m_nYMaxBoundingBox;
   int m_nWavePropagationModel;
   int m_nSimStartSec;
   int m_nSimStartMin;
   int m_nSimStartHour;
   int m_nSimStartDay;
   int m_nSimStartMonth;
   int m_nSimStartYear;
   int m_nDeepWaterWaveDataNTimeSteps;
   int m_nLogFileDetail;
   int m_nRunUpEquation;
   int m_nLevel;

   GDALDataType m_GDALWriteIntDataType;
   GDALDataType m_GDALWriteFloatDataType;

   long m_lGDALMaxCanWrite;
   long m_lGDALMinCanWrite;

   unsigned long m_ulIter;
   unsigned long m_ulTotTimestep;
   unsigned long m_ulRandSeed[NRNG];
   unsigned long m_ulNumCells;
   unsigned long m_ulThisIterNumSeaCells;
   unsigned long m_ulThisIterNumCoastCells;
   unsigned long m_ulThisIterNumPotentialPlatformErosionCells;
   unsigned long m_ulThisIterNumActualPlatformErosionCells;
   unsigned long m_ulThisIterNumPotentialBeachErosionCells;
   unsigned long m_ulThisIterNumActualBeachErosionCells;
   unsigned long m_ulThisIterNumBeachDepositionCells;
   unsigned long m_ulTotPotentialPlatformErosionOnProfiles;
   unsigned long m_ulTotPotentialPlatformErosionBetweenProfiles;
   unsigned long m_ulMissingValueBasementCells;

   double m_dDurationUnitsMult;
   double m_dNorthWestXExtCRS;
   double m_dNorthWestYExtCRS;
   double m_dSouthEastXExtCRS;
   double m_dSouthEastYExtCRS;
   double m_dExtCRSGridArea;
   double m_dCellSide;                                              // Length of cell side (in external CRS units)
   double m_dCellArea;                                              // Area of cell  (in external CRS units)
   double m_dCellDiagonal;                                          // Length of cell's diagonal (in external CRS units)
   double m_dInvCellSide;                                           // Inverse of m_dCellSide
   double m_dInvCellDiagonal;                                       // Inverse of m_dCellDiagonal
   double m_dSimDuration;                                           // Duration of simulation, in hours
   double m_dTimeStep;
   double m_dSimElapsed;                                            // Time simulated so far, in hours
   double m_dRegularSaveTime;
   double m_dRegularSaveInterval;
   double m_dUSaveTime[SAVEMAX];
   double m_dClkLast;                                               // Last value returned by clock()
   double m_dCPUClock;                                              // Total elapsed CPU time
   double m_dGeoTransform[6];
   double m_dSeaWaterDensity;
   double m_dOrigSWL;
   double m_dFinalSWL;
   double m_dDeltaSWLPerTimestep;
   double m_dThisIterSWL;
   double m_dThisIterDiffTotWaterLevel;
   double m_dThisIterDiffWaveSetupWaterLevel;
   double m_dThisIterDiffWaveSetupSurgeWaterLevel;
   double m_dThisIterDiffWaveSetupSurgeRunupWaterLevel;
   double m_dAccumulatedSeaLevelChange;
   double m_dMinSWL;
   double m_dMaxSWL;
   double m_dBreakingWaveHeight;
   double m_dC_0;                                                   // Deep water wave speed (m/s)
   double m_dL_0;                                                   // Deep water wave length (m)
   double m_dWaveDepthRatioForWaveCalcs;
   double m_dBreakingWaveHeightDepthRatio;
   double m_dAllCellsDeepWaterWaveHeight;
   double m_dAllCellsDeepWaterWaveAngle;
   double m_dAllCellsDeepWaterWavePeriod;
   double m_dMaxUserInputWaveHeight;
   double m_dMaxUserInputWavePeriod;                                // Used to constrain depth of closure
   double m_dR;
   double m_dD50Fine;
   double m_dD50Sand;
   double m_dD50Coarse;
   double m_dBeachSedimentDensity;
   double m_dBeachSedimentPorosity;
   double m_dFineErodibility;
   double m_dSandErodibility;
   double m_dCoarseErodibility;
   double m_dFineErodibilityNormalized;
   double m_dSandErodibilityNormalized;
   double m_dCoarseErodibilityNormalized;
   double m_dKLS;
   double m_dKamphuis;
   double m_dG;
   double m_dInmersedToBulkVolumetric;
   double m_dDepthOfClosure;
   double m_dCoastNormalAvgSpacing;                                 // In m
   double m_dCoastNormalLength;
   double m_dThisIterTotSeaDepth;
   double m_dThisIterPotentialPlatformErosion;
   double m_dThisIterActualPlatformErosionFineCons;
   double m_dThisIterActualPlatformErosionSandCons;
   double m_dThisIterActualPlatformErosionCoarseCons;
   double m_dThisIterPotentialBeachErosion;
   double m_dThisIterBeachErosionFine;
   double m_dThisIterBeachErosionSand;
   double m_dThisIterBeachErosionCoarse;
   double m_dThisIterBeachDepositionSand;
   double m_dThisIterBeachDepositionCoarse;
   double m_dThisIterFineSedimentToSuspension;
   double m_dThisIterPotentialSedLostBeachErosion;
   double m_dThisIterLeftGridUnconsFine;
   double m_dThisIterLeftGridUnconsSand;
   double m_dThisIterLeftGridUnconsCoarse;
   double m_dThisIterCliffCollapseFineErodedDuringDeposition;
   double m_dThisIterCliffCollapseSandErodedDuringDeposition;
   double m_dThisIterCliffCollapseCoarseErodedDuringDeposition;
   double m_dThisIterDepositionSandDiff;
   double m_dThisIterDepositionCoarseDiff;
   double m_dDepthOverDBMax;                                        // Used in erosion potential look-up function
   double m_dTotPotentialPlatformErosionOnProfiles;
   double m_dTotPotentialPlatformErosionBetweenProfiles;
   double m_dProfileMaxSlope;
   double m_dMaxBeachElevAboveSWL;
   double m_dCliffErosionResistance;
   double m_dNotchDepthAtCollapse;
   double m_dNotchBaseBelowSWL;
   double m_dCliffDepositionA;
   double m_dCliffDepositionPlanviewWidth;
   double m_dCliffTalusMinDepositionLength;
   double m_dMinCliffTalusHeightFrac;
   double m_dThisIterCliffCollapseErosionFineUncons;
   double m_dThisIterCliffCollapseErosionSandUncons;
   double m_dThisIterCliffCollapseErosionCoarseUncons;
   double m_dThisIterCliffCollapseErosionFineCons;
   double m_dThisIterCliffCollapseErosionSandCons;
   double m_dThisIterCliffCollapseErosionCoarseCons;
   double m_dThisIterUnconsSandCliffDeposition;
   double m_dThisIterUnconsCoarseCliffDeposition;
   double m_dCoastNormalRandSpaceFact;
   double m_dDeanProfileStartAboveSWL;
   double m_dMissingValue;
   double m_dWaveDataWrapHours;
   double m_dThisIterTopElevMax;
   double m_dThisIterTopElevMin;
   double m_dThisiterUnconsFineInput;
   double m_dThisiterUnconsSandInput;
   double m_dThisiterUnconsCoarseInput;
   double m_dStartIterSuspFine;
   double m_dStartIterUnconsFine;
   double m_dStartIterUnconsSand;
   double m_dStartIterUnconsCoarse;
   double m_dStartIterConsFine;
   double m_dStartIterConsSand;
   double m_dStartIterConsCoarse;

   // These grand totals are all long doubles. The aim is to minimize rounding errors when many very small numbers are added to a single much larger number, see e.g. http://www.ddj.com/cpp/184403224
   long double m_ldGTotPotentialPlatformErosion;
   long double m_ldGTotFineActualPlatformErosion;
   long double m_ldGTotSandActualPlatformErosion;
   long double m_ldGTotCoarseActualPlatformErosion;
   long double m_ldGTotPotentialSedLostBeachErosion;
   long double m_ldGTotActualFineLostBeachErosion;
   long double m_ldGTotActualSandLostBeachErosion;
   long double m_ldGTotActualCoarseLostBeachErosion;
   long double m_ldGTotSandSedLostCliffCollapse;
   long double m_ldGTotCoarseSedLostCliffCollapse;
   long double m_ldGTotCliffCollapseFine;
   long double m_ldGTotCliffCollapseSand;
   long double m_ldGTotCliffCollapseCoarse;
   long double m_ldGTotCliffTalusFineToSuspension;
   long double m_ldGTotCliffTalusSandDeposition;
   long double m_ldGTotCliffTalusCoarseDeposition;
   long double m_ldGTotCliffCollapseFineErodedDuringDeposition;
   long double m_ldGTotCliffCollapseSandErodedDuringDeposition;
   long double m_ldGTotCliffCollapseCoarseErodedDuringDeposition;
   long double m_ldGTotPotentialBeachErosion;
   long double m_ldGTotActualFineBeachErosion;
   long double m_ldGTotActualSandBeachErosion;
   long double m_ldGTotActualCoarseBeachErosion;
   long double m_ldGTotSandBeachDeposition;
   long double m_ldGTotCoarseBeachDeposition;
   long double m_ldGTotSuspendedSediment;
   long double m_ldGTotSandDepositionDiff;
   long double m_ldGTotCoarseDepositionDiff;
   long double m_ldGTotFineSedimentInput;
   long double m_ldGTotSandSedimentInput;
   long double m_ldGTotCoarseSedimentInput;

   string m_strCMEDir;
   string m_strCMEIni;
   string m_strMailAddress;
   string m_strDataPathName;
   string m_strRasterGISOutFormat;
   string m_strVectorGISOutFormat;
   string m_strInitialBasementDEMFile;
   string m_strInitialLandformFile;
   string m_strInterventionClassFile;
   string m_strInterventionHeightFile;
   string m_strInitialSuspSedimentFile;
   string m_strShapeFunctionFile;
   string m_strTideDataFile;
   string m_strLogFile;
   string m_strOutPath;
   string m_strOutFile;
   string m_strPalFile;
   string m_strGDALBasementDEMDriverCode; // Basement DEM (raster)
   string m_strGDALBasementDEMDriverDesc;
   string m_strGDALBasementDEMProjection;
   string m_strGDALBasementDEMDataType;
   string m_strGDALLDriverCode; // Initial landform class (raster)
   string m_strGDALLDriverDesc;
   string m_strGDALLProjection;
   string m_strGDALLDataType;
   string m_strGDALICDriverCode; // Initial intervention class (raster)
   string m_strGDALICDriverDesc;
   string m_strGDALICProjection;
   string m_strGDALICDataType;
   string m_strGDALIHDriverCode; // Initial intervention class (raster)
   string m_strGDALIHDriverDesc;
   string m_strGDALIHProjection;
   string m_strGDALIHDataType;
   string m_strGDALIWDriverCode; // Initial water depth (raster)
   string m_strGDALIWDriverDesc;
   string m_strGDALIWProjection;
   string m_strGDALIWDataType;
   string m_strGDALISSDriverCode; // Initial suspended sediment (raster)
   string m_strGDALISSDriverDesc;
   string m_strGDALISSProjection;
   string m_strGDALISSDataType;
   string m_strOGRDWWVDriverCode; // Initial deep water wave stations (vector)
   string m_strOGRDWWVGeometry;
   string m_strOGRDWWVDataType;
   string m_strOGRSedInputDriverCode; // Sediment input event locations (vector)
   string m_strOGRSedInputGeometry;
   string m_strOGRSedInputDataType;
   string m_strOGRFloodDriverCode; // Flood input locations (point or vector)
   string m_strOGRFloodGeometry;
   string m_strOGRFloodDataType;
   string m_strGDALRasterOutputDriverLongname;
   string m_strGDALRasterOutputDriverExtension;
   string m_strOGRVectorOutputExtension;
   string m_strRunName;
   string m_strDurationUnits;
   string m_strDeepWaterWaveStationsShapefile;
   string m_strDeepWaterWavesTimeSeriesFile;
   string m_strSedimentInputEventShapefile;
   string m_strSedimentInputEventTimeSeriesFile;
      // string m_strLevel;
   string m_strFloodLocationShapefile;

   struct RandState
   {
      unsigned long s1, s2, s3;
   } m_ulRState[NRNG];

   time_t m_tSysStartTime;
   time_t m_tSysEndTime;

   ofstream OutStream;
   ofstream SeaAreaTSStream;
   ofstream StillWaterLevelTSStream;
   ofstream ErosionTSStream;
   ofstream CliffCollapseErosionTSStream;
   ofstream CliffCollapseDepositionTSStream;
   ofstream CliffCollapseNetTSStream;
   ofstream BeachErosionTSStream;
   ofstream BeachDepositionTSStream;
   ofstream BeachSedimentChangeNetTSStream;
   ofstream SedLoadTSStream;
   ofstream FloodSetupSurgeTSStream;
   ofstream FloodSetupSurgeRunupTSStream;

   vector<bool> m_bConsChangedThisIter;
   vector<bool> m_bUnconsChangedThisIter;

   vector<int> m_VnProfileToSave;

   //! ID for deep water wave station, this corresponds with the ID in the wave time series file
   vector<int> m_VnDeepWaterWaveStationID;

   //! ID for sediment input location, this corresponds with the ID in the sediment input time series file
   vector<int> m_VnSedimentInputLocationID;

   //! ID for flood location
   vector<int> m_VnFloodLocationID;

   //! Savitzky-Golay shift index for the coastline vector(s)
   vector<int> m_VnSavGolIndexCoast;

   vector<unsigned long>
      m_VulProfileTimestep,
      m_VlDeepWaterWaveValuesAtTimestep; // Calculate deep water wave values at these timesteps

   vector<double> m_VdSliceElev;

   //! For erosion potential lookup
   vector<double> m_VdErosionPotential;

   //! For erosion potential lookup
   vector<double> m_VdDepthOverDB;

   //! Savitzky-Golay filter coefficients for the coastline vector(s)
   vector<double> m_VdSavGolFCRWCoast;

   //! Savitzky-Golay filter coefficients for the profile vectors
   vector<double> m_VdSavGolFCGeomProfile;

   //! Tide data: one record per timestep, is the change (m) from still water level for that timestep
   vector<double> m_VdTideData;

   //! X co-ordinate (grid CRS) for deep water wave station
   vector<double> m_VdDeepWaterWaveStationX;

   //! Y co-ordinate (grid CRS) for deep water wave station
   vector<double> m_VdDeepWaterWaveStationY;

   //! This-iteration wave height at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationHeight;

   //! This-iteration wave orientation at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationAngle;

   //! This-iteration wave period at deep water wave station
   vector<double> m_VdThisIterDeepWaterWaveStationPeriod;

   //! Time series of wave heights at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationHeight;

   //! Time series of wave orientation at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationAngle;

   //! Time series of wave period at deep water wave station
   vector<double> m_VdTSDeepWaterWaveStationPeriod;

   //! X co-ordinate (grid CRS) for sediment input event
   vector<double> m_VdSedimentInputLocationX;

   //! X co-ordinate (grid CRS) for sediment input event
   vector<double> m_VdSedimentInputLocationY;

   //! X co-ordinate (grid CRS) for total water level flooding
   vector<double> m_VdFloodLocationX;

   //! X co-ordinate (grid CRS) for total water level flooding
   vector<double> m_VdFloodLocationY;

   vector<string> m_VstrInitialFineUnconsSedimentFile;
   vector<string> m_VstrInitialSandUnconsSedimentFile;
   vector<string> m_VstrInitialCoarseUnconsSedimentFile;
   vector<string> m_VstrInitialFineConsSedimentFile;
   vector<string> m_VstrInitialSandConsSedimentFile;
   vector<string> m_VstrInitialCoarseConsSedimentFile;
   vector<string> m_VstrGDALIUFDriverCode;
   vector<string> m_VstrGDALIUFDriverDesc;
   vector<string> m_VstrGDALIUFProjection;
   vector<string> m_VstrGDALIUFDataType;
   vector<string> m_VstrGDALIUSDriverCode;
   vector<string> m_VstrGDALIUSDriverDesc;
   vector<string> m_VstrGDALIUSProjection;
   vector<string> m_VstrGDALIUSDataType;
   vector<string> m_VstrGDALIUCDriverCode;
   vector<string> m_VstrGDALIUCDriverDesc;
   vector<string> m_VstrGDALIUCProjection;
   vector<string> m_VstrGDALIUCDataType;
   vector<string> m_VstrGDALICFDriverCode;
   vector<string> m_VstrGDALICFDriverDesc;
   vector<string> m_VstrGDALICFProjection;
   vector<string> m_VstrGDALICFDataType;
   vector<string> m_VstrGDALICSDriverCode;
   vector<string> m_VstrGDALICSDriverDesc;
   vector<string> m_VstrGDALICSProjection;
   vector<string> m_VstrGDALICSDataType;
   vector<string> m_VstrGDALICCDriverCode;
   vector<string> m_VstrGDALICCDriverDesc;
   vector<string> m_VstrGDALICCProjection;
   vector<string> m_VstrGDALICCDataType;

   //! Pointer to the raster grid object
   CGeomRasterGrid* m_pRasterGrid;

   //! The coastline objects
   vector<CRWCoast> m_VCoast;

   // vector<CRWCoast> m_VFloodWaveSetup;
   vector<CRWCoast> m_VFloodWaveSetupSurge;
   vector<CRWCoast> m_VFloodWaveSetupSurgeRunup;

   //! Pointers to coast polygon objects
   vector<CGeomCoastPolygon*> m_pVCoastPolygon;

   //! Edge cells
   vector<CGeom2DIPoint> m_VEdgeCell;

   //! The grid edge that each edge cell belongs to
   vector<int> m_VEdgeCellEdge;

   //! The location to compute the total water level for flooding
   vector<int> m_VCellFloodLocation;

   //! Sediment input events
   vector<CSedInputEvent*> m_pVSedInputEvent;

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
   int nSaveProfile(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*);
   bool bWriteProfileData(int const, int const, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<CGeom2DIPoint>* const, vector<double> const*) const;
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
   int nAssignAllCoastalLandforms(void);
   int nAssignNonCoastlineLandforms(void);
   int nDoAllPropagateWaves(void);
   int nDoAllShorePlatFormErosion(void);
   int nDoAllWaveEnergyToCoastLandforms(void);
   int nDoCliffCollapse(int const, CRWCliff*, double&, double&, double&, double&, double&);
   int nDoCliffCollapseDeposition(int const, CRWCliff const*, double const, double const, double const, double const);
   int nUpdateGrid(void);

   // Lower-level simulation routines
   void FindAllSeaCells(void);
   int FindAllInundatedCells(void);
   void FloodFillSea(int const, int const);
   void FloodFillLand(int const, int const);
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
   static bool bCheckForIntersection(CGeomProfile *const, CGeomProfile* const, int&, int&, double&, double&, double&, double&);
   void MergeProfilesAtFinalLineSegments(int const, int const, int const, int const, int const, double const, double const, double const, double const);
   void TruncateOneProfileRetainOtherProfile(int const, int const, int const, double const, double const, int const, int const, bool const);
   int nInsertPointIntoProfilesIfNeededThenUpdate(int const, int const, double const, double const, int const, int const, int const, bool const);
   void TruncateProfileAndAppendNew(int const, int const, int const, vector<CGeom2DPoint> const*, vector<vector<pair<int, int>>> const*);
   void RasterizeProfile(int const, int const, vector<CGeom2DIPoint>*, vector<bool>*, bool&, bool&, bool&, bool&, bool&);
   static void CalcDeanProfile(vector<double>*, double const, double const, double const, bool const, int const, double const);
   static double dSubtractProfiles(vector<double> const*, vector<double> const*, vector<bool> const*);
   void RasterizeCliffCollapseProfile(vector<CGeom2DPoint> const*, vector<CGeom2DIPoint>*) const;
   int nCalcPotentialPlatformErosionOnProfile(int const, int const);
   int nCalcPotentialPlatformErosionBetweenProfiles(int const, int const, int const);
   void ConstructParallelProfile(int const, int const, int const, int const, int const, vector<CGeom2DIPoint>* const, vector<CGeom2DIPoint>*, vector<CGeom2DPoint> *);
   double dCalcBeachProtectionFactor(int const, int const, double const);
   void FillInBeachProtectionHoles(void);
   void FillPotentialPlatformErosionHoles(void);
   void DoActualPlatformErosionOnCell(int const, int const);
   double dLookUpErosionPotential(double const) const;
   static CGeom2DPoint PtChooseEndPoint(int const, CGeom2DPoint const*, CGeom2DPoint const*, double const, double const, double const, double const);
   int nGetCoastNormalEndPoint(int const, int const, int const, CGeom2DPoint const*, double const, CGeom2DPoint*, CGeom2DIPoint*);
   int nLandformToGrid(int const, int const);
   int nCalcWavePropertiesOnProfile(int const, int const, int const, vector<double>*, vector<double>*, vector<double>*, vector<double>*, vector<bool>*);
   int nGetThisProfileElevationVectorsForCShore(int const, int const, int const, vector<double>*, vector<double>*, vector<double>*);
   int nCreateCShoreInfile(int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, int const, double const, double const, double const, double const, double const, double const, double const, double const, vector<double> const*, vector<double> const*, vector<double> const*);
   int nReadCShoreOutput(int const, string const*, int const, int const, vector<double> const*, vector<double>*);   
   static void InterpolateCShoreOutput(vector<double> const*, int const, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double> const*, vector<double>*, vector<double>*, vector<double>*, vector<double>*);
   static double dCalcWaveAngleToCoastNormal(double const, double const, int const);
   void CalcCoastTangents(int const);
   void InterpolateWavePropertiesBetweenProfiles(int const, int const, int const);
   void InterpolateWaveHeightToCoastPoints(int const);
   // void InterpolateWavePropertiesToCells(int const, int const, int const);
   void ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const, int const);
   static double dCalcCurvature(int const, CGeom2DPoint const*, CGeom2DPoint const*, CGeom2DPoint const*);
   void CalcD50AndFillWaveCalcHoles(void);
   int nDoAllShadowZones(void);
   static bool bOnOrOffShoreAndUpOrDownCoast(double const, double const, int const, bool&);
   static CGeom2DIPoint PtiFollowWaveAngle(CGeom2DIPoint const*, double const, double&);
   // int nFindAllShadowZones(void);
   int nFloodFillShadowZone(int const, CGeom2DIPoint const*, CGeom2DIPoint const*, CGeom2DIPoint const*);
   void DoShadowZoneAndDownDriftZone(int const, int const, int const, int const);
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
   // int nEstimateErosionOnPolygon(int const, int const, double const, double&, double&, double&);
   // int nEstimateUnconsErosionOnParallelProfile(/*int const, int const,*/ int const, int const, /* int const, */ int const, vector<CGeom2DIPoint> const*, vector<double> const*, double&, double&, double&, double&, double&);
   int nDoParallelProfileUnconsErosion( int const, int const,  int const, int const, int const,  int const,  int const, vector<CGeom2DIPoint> const*, vector<double> const*, double&, double&, double&);
   // void EstimateUnconsErosionOnCell(int const, int const, int const, double const, double&, double&, double&);
   void ErodeCellBeachSedimentSupplyLimited(int const, int const, int const, int const, double const, double&);
   // int nEstimateMovementUnconsToAdjacentPolygons(int const, int const);
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
   // double dExtCRSXToGridCentroidX(double const) const;
   // double dExtCRSYToGridCentroidY(double const) const;
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
   // vector<CGeom2DPoint> VGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DPoint PtGetPerpendicular(CGeom2DPoint const*, CGeom2DPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(CGeom2DIPoint const*, CGeom2DIPoint const*, double const, int const);
   static CGeom2DIPoint PtiGetPerpendicular(int const, int const, int const, int const, double const, int const);
   static CGeom2DPoint PtAverage(CGeom2DPoint const*, CGeom2DPoint const*);
   static CGeom2DPoint PtAverage(vector<CGeom2DPoint>*);
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
   CGeomLine LSmoothCoastSavitzkyGolay(CGeomLine*, int const, int const) const;
   CGeomLine LSmoothCoastRunningMean(CGeomLine*) const;
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
   static vector<string> *VstrSplit(string const*, char const, vector<string>*);
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
   void WritePolygonActualMovement(int const, vector<vector<int> > const&);

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
