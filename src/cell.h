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
    bool
        m_bInContiguousSea,                  // Is a sea cell, contiguous with other sea cells
        m_bInContiguousFlood,
        m_bIsInActiveZone,
        m_bCoastline,
        m_bFloodLine,
        m_bWaveFlood,
        m_bCheckCell,
        m_bCheckFloodCell,
        m_bShadowBoundary,
        m_bPossibleCoastStartCell,
        m_bPossibleFloodStartCell,
        m_bFloodBySetupSurge,
        m_bFloodBySetupSurgeRunup;

    int
        m_nBoundingBoxEdge,
        m_nPolygonID,
        m_nCoastlineNormal,
        m_nShadowZoneNumber,
        m_nDownDriftZoneNumber;

    double
        m_dLocalConsSlope,                   // As used in erosion calcs (really just for display purposes)
        m_dBasementElevation,                // Elevation of basement surface (m)
        m_dSeaDepth,                         // Depth of still water (m), is zero if not inundated
        m_dTotSeaDepth,                      // Total depth of still water (m) since beginning of simulation (used to calc average)
        m_dWaveHeight,                       // Wave height (m)
        m_dTotWaveHeight,                    // Total wave height (m) (used to calc average)
        m_dWaveAngle,                        // Wave orientation
        m_dWavePeriod,                       // Wave period (s)
        m_dTotWaveAngle,                     // Total wave orientation  (used to calc average)
        m_dDeepWaterWaveHeight,              // Wave height if this is a deep water cell
        m_dDeepWaterWaveAngle,               // Wave orientation if this is a deep water cell
        m_dDeepWaterWavePeriod,              // Wave period if this is a deep water cell
        m_dBeachProtectionFactor,            // Only meaningful if in zone of platform erosion. 0 is fully protected, 1 = no protection
        m_dSuspendedSediment,                // Suspended sediment as depth equivalent (m)
        m_dTotSuspendedSediment,             // Total depth of suspended sediment (m) since simulation start (used to calc average)
        m_dPotentialPlatformErosionThisIter, // Depth of sediment on the shore platform that could be eroded this timestep, if no supply-limitation
        m_dTotPotentialPlatformErosion,      // Total depth of sediment eroded from the shore platform, if no supply-limitation
        m_dActualPlatformErosionThisIter,    // Depth of sediment actually eroded from the shore platform this timestep
        m_dTotActualPlatformErosion,         // Total depth of sediment actually eroded from the shore platform
        m_dCliffCollapseFineThisIter,        // Depth of fine sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
        m_dCliffCollapseSandThisIter,        // Depth of sand sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
        m_dCliffCollapseCoarseThisIter,      // Depth of coarse sediment (consolidated and unconsolidated) removed via cliff collapse this timestep
        m_dTotFineCliffCollapse,             // Total depth of fine sediment (consolidated and unconsolidated) removed via cliff collapse
        m_dTotSandCliffCollapse,             // Total depth of sand sediment (consolidated and unconsolidated) removed via cliff collapse
        m_dTotCoarseCliffCollapse,           // Total depth of coarse sediment (consolidated and unconsolidated) removed via cliff collapse
        m_dTalusSandDepositionThisIter,      // Depth of unconsolidated sand sediment deposited as a result of cliff collapse this timestep
        m_dTotTalusSandDeposition,           // Total depth of unconsolidated sand sediment deposited as a result of cliff collapse
        m_dTalusCoarseDepositionThisIter,    // Depth of unconsolidated coarse sediment deposited as a result of cliff collapse this timestep
        m_dTotTalusCoarseDeposition,         // Total depth of unconsolidated coarse sediment deposited as a result of cliff collapse
        m_dPotentialBeachErosionThisIter,    // Depth of unconsolidated beach sediment that could be eroded this timestep, if no supply-limitation
        m_dTotPotentialBeachErosion,         // Total depth of unconsolidated beach sediment eroded, if no supply-limitation
        m_dActualBeachErosionThisIter,       // Depth of unconsolidated beach sediment actually eroded this timestep
        m_dTotActualBeachErosion,            // Total depth of unconsolidated beach sediment actually eroded
        m_dBeachDepositionThisIter,          // Depth of unconsolidated beach sediment deposited this timestep
        m_dTotBeachDeposition,               // Total depth of unconsolidated beach sediment deposited
        m_dUnconsD50,                        // d50 of unconsolidated sediment on top layer with unconsolidated sediment depth > 0
        m_dInterventionHeight;               // Height of intervention structure

    // This cell's landform data
    CRWCellLandform m_Landform;

    // Initialize these as empty vectors
    vector<CRWCellLayer> m_VLayerAboveBasement; // Number of layers NOT including the basement. Layer 0 is the lowest
    vector<double> m_VdAllHorizonTopElev;       // Number of layer-top elevations (inc. that of the basement, which is m_VdAllHorizonTopElev[0]); size 1 greater than size of m_VLayerAboveBasement

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

    CRWCellLandform *pGetLandform(void);

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
    CRWCellLayer *pGetLayerAboveBasement(int const);
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
