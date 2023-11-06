/*!
 * \class CGeomCell
 * \brief Geometry class for the cell objects which comprise the raster grid
 * \details TODO This is a more detailed description of the CGeomCell class.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2023
 * \copyright GNU General Public License
 *
 * \file cell.h
 * \brief Contains CGeomCell definitions
 *
 */

#ifndef CELL_H
#define CELL_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "cell_landform.h"
#include "cell_layer.h"
#include "raster_grid.h"

class CGeomCell
{
   friend class CSimulation;

private:
   //! Switch to indicate if this is a sea cell, contiguous with other sea cells
   bool m_bInContiguousSea;

   //!
   bool m_bInContiguousFlood;

   //!
   bool m_bIsInActiveZone;

   //!
   bool m_bCoastline;

   //!
   bool m_bFloodLine;

   //!
   bool m_bWaveFlood;

   //!
   bool m_bCheckCell;

   //!
   bool m_bCheckFloodCell;

   //!
   bool m_bShadowBoundary;

   //!
   bool m_bPossibleCoastStartCell;

   //!
   bool m_bPossibleFloodStartCell;

   //!
   bool m_bFloodBySetupSurge;

   //!
   bool m_bFloodBySetupSurgeRunup;

   //!
   int m_nBoundingBoxEdge;

   //!
   int m_nPolygonID;

   //!
   int m_nCoastlineNormal;

   //!
   int m_nShadowZoneNumber;

   //!
   int m_nDownDriftZoneNumber;

   //! As used in erosion calcs (really just for display purposes)
   double m_dLocalConsSlope;

   //! Elevation of basement surface (m)
   double m_dBasementElevation;

   //! Depth of still water (m), is zero if not inundated
   double m_dSeaDepth;

   //! Total depth of still water (m) since beginning of simulation (used to calc average)
   double m_dTotSeaDepth;

   //! Wave height (m)
   double m_dWaveHeight;

   //! Total wave height (m) (used to calc average)
   double m_dTotWaveHeight;

   //! Wave orientation
   double m_dWaveAngle;

   //! Wave period (s)
   double m_dWavePeriod;

   //! Total wave orientation  (used to calc average)
   double m_dTotWaveAngle;

   //! Wave height if this is a deep water cell
   double m_dDeepWaterWaveHeight;

   //! Wave orientation if this is a deep water cell
   double m_dDeepWaterWaveAngle;

   //! Wave period if this is a deep water cell
   double m_dDeepWaterWavePeriod;

   //! Only meaningful if in zone of platform erosion. 0 is fully protected; 1 = no protection
   double m_dBeachProtectionFactor;

   //! Suspended sediment as depth equivalent (m)
   double m_dSuspendedSediment;

   //! Total depth of suspended sediment (m) since simulation start (used to calc average)
   double m_dTotSuspendedSediment;

   //! Depth of sediment on the shore platform that could be eroded this timestep, if no supply-limitation
   double m_dPotentialPlatformErosionThisIter;

   //! Total depth of sediment eroded from the shore platform, if no supply-limitation
   double m_dTotPotentialPlatformErosion;

   //! Depth of sediment actually eroded from the shore platform this timestep
   double m_dActualPlatformErosionThisIter;

   //! Total depth of sediment actually eroded from the shore platform
   double m_dTotActualPlatformErosion;

   //! Depth of fine sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
   double m_dCliffCollapseFineThisIter;

   //! Depth of sand sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
   double m_dCliffCollapseSandThisIter;

   //! Depth of coarse sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
   double m_dCliffCollapseCoarseThisIter;

   //! Total depth of fine sediment (consolidated and unconsolidated) removed via cliff collapse
   double m_dTotFineCliffCollapse;

   //! Total depth of sand sediment (consolidated and unconsolidated) removed via cliff collapse
   double m_dTotSandCliffCollapse;

   //! Total depth of coarse sediment (consolidated and unconsolidated) removed via cliff collapse
   double m_dTotCoarseCliffCollapse;

   //! Depth of unconsolidated sand sediment deposited as a result of cliff collapse this timestep
   double m_dTalusSandDepositionThisIter;

   //! Total depth of unconsolidated sand sediment deposited as a result of cliff collapse
   double m_dTotTalusSandDeposition;

   //! Depth of unconsolidated coarse sediment deposited as a result of cliff collapse this timestep
   double m_dTalusCoarseDepositionThisIter;

   //! Total depth of unconsolidated coarse sediment deposited as a result of cliff collapse
   double m_dTotTalusCoarseDeposition;

   //! Depth of unconsolidated beach sediment that could be eroded this timestep, if no supply-limitation
   double m_dPotentialBeachErosionThisIter;

   //! Total depth of unconsolidated beach sediment eroded; if no supply-limitation
   double m_dTotPotentialBeachErosion;

   //! Depth of unconsolidated beach sediment actually eroded this timestep
   double m_dActualBeachErosionThisIter;

   //! Total depth of unconsolidated beach sediment actually eroded
   double m_dTotActualBeachErosion;

   //! Depth of unconsolidated beach sediment deposited this timestep
   double m_dBeachDepositionThisIter;

   //! Total depth of unconsolidated beach sediment deposited
   double m_dTotBeachDeposition;

   //! d50 of unconsolidated sediment on top layer with unconsolidated sediment depth > 0
   double m_dUnconsD50;

   //! Height of intervention structure
   double m_dInterventionHeight;

   //! This cell's landform data
   CRWCellLandform m_Landform;

   // Initialize these as empty vectors
   //! Number of layers NOT including the basement. Layer 0 is the lowest
   vector<CRWCellLayer> m_VLayerAboveBasement;

   //! Number of layer-top elevations (inc. that of the basement, which is m_VdAllHorizonTopElev[0]); size 1 greater than size of m_VLayerAboveBasement
   vector<double> m_VdAllHorizonTopElev;

public:
    static CGeomRasterGrid *m_pGrid;

    CGeomCell();
    ~CGeomCell(void);

    void SetInContiguousSea(void);
    bool bIsInContiguousSea(void) const;

    void SetInContiguousFlood(void);
    void UnSetInContiguousFlood(void);
    void SetFloodBySetupSurge(void);
    bool bIsFloodBySetupSurge(void) const;
    void SetFloodBySetupSurgeRunup(void);
    bool bIsFloodBySetupSurgeRunup(void) const;
    bool bIsInContiguousFlood(void) const;

    // void SetActualBeachErosionEstimated(void);
    // bool bGetActualBeachErosionEstimated(void) const;

    void SetInActiveZone(bool const);
    bool bIsInActiveZone(void) const;
    bool bPotentialPlatformErosion(void) const;
    //    bool bActualPlatformErosion(void) const;
    void SetAsCoastline(bool const);
    bool bIsCoastline(void) const;
    void SetAsFloodLine(bool const);
    bool bIsFloodLine(void) const;

    void SetProfile(int const);
    int nGetProfile(void) const;
    bool bIsProfile(void) const;

    void SetShadowZoneBoundary(void);
    bool bIsShadowZoneBoundary(void) const;

    void SetBoundingBoxEdge(int const);
    int nGetBoundingBoxEdge(void) const;
    bool bIsBoundingBoxEdge(void) const;

    void SetPossibleCoastStartCell(void);
    bool bIsPossibleCoastStartCell(void) const;

    void SetPossibleFloodStartCell(void);
    bool bIsPossibleFloodStartCell(void) const;

    void SetPolygonID(int const);
    int nGetPolygonID(void) const;

    CRWCellLandform* pGetLandform(void);

    void SetWaveFlood(void);
    bool bIsElevLessThanWaterLevel(void) const;

    void SetCheckCell(void);
    bool bIsCellCheck(void) const;

    void SetCheckFloodCell(void);
    void UnSetCheckFloodCell(void);
    bool bIsCellFloodCheck(void) const;

    void SetLocalConsSlope(double const);
    double dGetLocalConsSlope(void) const;

    void SetBasementElev(double const);
    double dGetBasementElev(void) const;
    bool bBasementElevIsMissingValue(void) const;

    double dGetVolEquivSedTopElev(void) const;
    double dGetSedimentTopElev(void) const;
    double dGetSedimentPlusInterventionTopElev(void) const;
    double dGetOverallTopElev(void) const;

    bool bIsInundated(void) const;
    double dGetIterSWL(void) const;
    double dGetIterTotWaterLevel(void) const;
    bool bIsSeaIncBeach(void) const;
    void SetSeaDepth(void);
    double dGetSeaDepth(void) const;
    void InitCell(void);
    double dGetTotSeaDepth(void) const;

    void SetWaveHeight(double const);
    double dGetWaveHeight(void) const;
    double dGetTotWaveHeight(void) const;
    void SetWaveAngle(double const);
    double dGetWaveAngle(void) const;
    double dGetTotWaveAngle(void) const;

    void SetCellDeepWaterWaveHeight(double const);
    double dGetCellDeepWaterWaveHeight(void) const;
    void SetCellDeepWaterWaveAngle(double const);
    double dGetCellDeepWaterWaveAngle(void) const;
    void SetCellDeepWaterWavePeriod(double const);
    double dGetCellDeepWaterWavePeriod(void) const;

    void SetWaveValuesToDeepWaterWaveValues(void);

    void SetBeachProtectionFactor(double const);
    double dGetBeachProtectionFactor(void) const;

    void SetSuspendedSediment(double const);
    void AddSuspendedSediment(double const);
    double dGetSuspendedSediment(void) const;
    double dGetTotSuspendedSediment(void) const;

    int nGetTopNonZeroLayerAboveBasement(void) const;
    int nGetTopLayerAboveBasement(void) const;

    double dGetConsSedTopForLayerAboveBasement(int const) const;
    CRWCellLayer* pGetLayerAboveBasement(int const);
    void AppendLayers(int const);
    void CalcAllLayerElevsAndD50(void);
    int nGetLayerAtElev(double const) const;
    double dCalcLayerElev(const int);

    double dGetTotConsFineThickConsiderNotch(void) const;
    double dGetTotUnconsFineThickness(void) const;
    double dGetTotConsSandThickConsiderNotch(void) const;
    double dGetTotUnconsSandThickness(void) const;
    double dGetTotConsCoarseThickConsiderNotch(void) const;
    double dGetTotUnconsCoarseThickness(void) const;

    double dGetTotConsThickness(void) const;
    double dGetTotUnconsThickness(void) const;
    double dGetTotAllSedThickness(void) const;

    void SetPotentialPlatformErosion(double const);
    double dGetPotentialPlatformErosion(void) const;
    double dGetTotPotentialPlatformErosion(void) const;

    void SetActualPlatformErosion(double const);
    double dGetActualPlatformErosion(void) const;
    double dGetTotActualPlatformErosion(void) const;

    void IncrCliffCollapseErosion(double const, double const, double const);
    double dGetThisIterCliffCollapseErosionFine(void) const;
    double dGetThisIterCliffCollapseErosionSand(void) const;
    double dGetThisIterCliffCollapseErosionCoarse(void) const;
    double dGetTotCliffCollapseFine(void) const;
    double dGetTotCliffCollapseSand(void) const;
    double dGetTotCliffCollapseCoarse(void) const;

    void AddSandTalusDeposition(double const);
    double dGetThisIterCliffCollapseSandTalusDeposition(void) const;
    double dGetTotSandTalusDeposition(void) const;
    void AddCoarseTalusDeposition(double const);
    double dGetThisIterCliffCollapseCoarseTalusDeposition(void) const;
    double dGetTotCoarseTalusDeposition(void) const;

    void SetPotentialBeachErosion(double const);
    double dGetPotentialBeachErosion(void) const;
    double dGetTotPotentialBeachErosion(void) const;
    void SetActualBeachErosion(double const);
    double dGetActualBeachErosion(void) const;
    double dGetTotActualBeachErosion(void) const;
    //    bool bActualBeachErosionThisIter(void) const;

    void IncrBeachDeposition(double const);
    double dGetBeachDeposition(void) const;
    double dGetTotBeachDeposition(void) const;
    //    bool bBeachDepositionThisIter(void) const;

    bool bBeachErosionOrDepositionThisIter(void) const;

    double dGetUnconsD50(void) const;

    void SetInterventionClass(int const);
    int nGetInterventionClass(void) const;
    void SetInterventionHeight(double const);
    double dGetInterventionHeight(void) const;
    double dGetInterventionTopElev(void) const;

    void SetShadowZoneNumber(int const);
    int nGetShadowZoneNumber(void) const;
    bool bIsinThisShadowZone(int const) const;
    bool bIsinAnyShadowZone(void) const;
    void SetDownDriftZoneNumber(int const);
    int nGetDownDriftZoneNumber(void) const;
};
#endif // CELL_H
