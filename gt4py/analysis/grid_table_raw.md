# Grids shared between WaSiM modules: raw output of coupling.py (unverified, see ../coupling.md)

## Cross-module writes

Grids a module writes although another module (or the input) owns them.

| Grid | Internal name | Created by | Written by | Element writes | Pointer(s) |
| --- | --- | --- | --- | --- | --- |
| `_Reale_Verdunstung_` | real_transpiration | (input grid) | abfluss.cc | 2 | EI, ETR, conc, qu |
| `_ponding_storage_` | ponding_storage | unsatzon.cc | abfluss.cc | 2 | POND |
| `_ponding_storage_top_` | ponding_storage_top | (input grid) | abfluss.cc | 2 | POND |
| `_river_links_` | river_links | (input grid) | defltuzr.cc | 7 | RIVLNK, p |
| `_Bodenspeicher_` | soil_storage | (input grid) | etpetr.cc | 4 | SB |
| `_Globalstrahlung_` | global_radiation | (input grid) | etpetr.cc | 1 | RG |
| `_Nacht_Temperatur_` | temperature_night | expokorr.cc | etpetr.cc | 1 | Tnacht |
| `_Nettostrahlung_` | net_radiation | (input grid) | etpetr.cc | 1 | RG |
| `_SnowEvaporation_onGlac_` | SnowEvaporation_onGlacier | snowmod.cc | etpetr.cc | 6 | EVRSGL |
| `_SolidSnowStorage_onGlacier_` | SolidSnowStorage_onGlacier | snowmod.cc | etpetr.cc | 5 | SSNOWGL |
| `_ponding_storage_` | ponding_storage | unsatzon.cc | etpetr.cc | 2 | POND |
| `_ponding_storage_top_` | ponding_storage_top | (input grid) | etpetr.cc | 2 | POND |
| `_Basisabfluss_` | baseflow | unsatzon.cc | gwflow.cc | 2 | DELTAGW |
| `_Landnutzung_` | landuse | defltuzr.cc, unsatzon.cc | interzep.cc | 13 | EIP, ETP, ETR, LAI, LUSE, RSC, RSI, VBG, Z0 |
| `_Niederschlag_` | precipitation | (input grid) | irrigate.cc | 1 | P |
| `_irrigation_amounts_` | irrigation_amounts | unsatzon.cc | irrigate.cc | 16 | IRRIG |
| `_Niederschlag_` | precipitation | (input grid) | meteoint.cc | 5 | PRECIP |
| `_Niederschlag_` | precipitation | (input grid) | siltingup.cc | 3 | P |
| `_SURFSTORSILTINGUP_` | SurfStorSiltingUp | (input grid) | siltingup.cc | 3 | SSSU |
| `_Schneeausfluss_` | snowcover_outflow | snowmod.cc | siltingup.cc | 3 | P |
| `_Albedo_` | albedo | (input grid) | snowmod.cc | 4 | Albedo |
| `_GlacierCodes_` | GlacierCodes | (input grid) | snowmod.cc | 14 | GLACID |
| `_GlacierizedCells_` | GlacierizedCells | (input grid) | snowmod.cc | 8 | ICE |
| `_Niederschlag_` | precipitation | (input grid) | snowmod.cc | 1 | P |
| `_SnowDensity_` | SnowDensity | unsatzon.cc | snowmod.cc | 2 | refDensSnow |
| `_SnowLayerThickness_` | SnowLayerThickness | unsatzon.cc | snowmod.cc | 3 | refThickSnow |
| `_SnowWaterEquivalentLiquid_` | SnowWaterEquivalentLiquid | unsatzon.cc | snowmod.cc | 3 | refSWEL |
| `_SnowWaterEquivalentSolid_` | SnowWaterEquivalentSolid | unsatzon.cc | snowmod.cc | 5 | refSWES |
| `_Soil_Temperature_Stack_Snow_` | Temperature_Stack_Snow | unsatzon.cc | snowmod.cc | 2 | refTemperSnow |
| `_ice_firn_` | ice_firn | (input grid) | snowmod.cc | 8 | ICE |
| `_upperIndexSnow_` | upperIndexSnow | unsatzon.cc | snowmod.cc | 1 | refUpperIndexSnow |
| `_Direktabfluss_` | directflow | unsatzon.cc | surfrout.cc | 7 | QDgrid |
| `_MPIRegions_` | MPIRegions | defltuzr.cc | surfrout.cc | 1 | MPIREG |
| `_SURFSTORSILTINGUP_` | SurfStorSiltingUp | (input grid) | surfrout.cc | 6 | SSSU |
| `_Albedo_` | albedo | (input grid) | unsatzon.cc | 1 | ALBEDO |
| `_Bodenarten_` | soil_types | (input grid) | unsatzon.cc | 1 | SOIL |
| `_DIRECTRUNOFFSILTINGUP_` | DirectRunoffSiltingUp | siltingup.cc | unsatzon.cc | 2 | QDSU |
| `_Drainagedistance_` | drainage_distance | (input grid) | unsatzon.cc | 1 | DRAINDIST |
| `_EnergyBalance_` | EnergyBalance | snowmod.cc | unsatzon.cc | 3 | EBAL |
| `_Gefaelle_` | slope_angle | (input grid) | unsatzon.cc | 2 | SLP |
| `_GlacierizedCellsOld_` | GlacierizedCellsOld | snowmod.cc | unsatzon.cc | 1 | ICEOLD |
| `_GlacierizedCells_` | GlacierizedCells | (input grid) | unsatzon.cc | 1 | ICE |
| `_HeatFlux_Lower_Boundary_Condition_` | HeatFlux_Lower_Boundary_Condition | (input grid) | unsatzon.cc | 1 | QLOWBND |
| `_IceWedgeLowerDepth_` | IceWedgeLowerDepth | (input grid) | unsatzon.cc | 2 | ICEWEDGELOWERDEPTH |
| `_IceWedgeUpperDepth_` | IceWedgeUpperDepth | (input grid) | unsatzon.cc | 2 | ICEWEDGEUPPERDEPTH |
| `_LAKE_CODES_` | lake_codes | (input grid) | unsatzon.cc | 1 | LAKE |
| `_MPIRegions_` | MPIRegions | defltuzr.cc | unsatzon.cc | 1 | MPIREG |
| `_Niederschlag_` | precipitation | (input grid) | unsatzon.cc | 4 | P |
| `_PermafrostThawDepth_` | PermafrostThawDepth | (input grid) | unsatzon.cc | 1 | THAWDEPTH |
| `_Potential_Evaporation_` | potential_evaporation | etpetr.cc | unsatzon.cc | 1 | EVAP |
| `_Reale_Evaporation_` | real_evaporation | etpetr.cc | unsatzon.cc | 4 | EVAR |
| `_Regenanteil_` | rain_fraction | snowmod.cc | unsatzon.cc | 1 | P_RAIN |
| `_SB_1_REL_grid_` | SB_1_REL_grid | (input grid) | unsatzon.cc | 4 | SB_1_REL |
| `_SURFSTORSILTINGUP_` | SurfStorSiltingUp | (input grid) | unsatzon.cc | 2 | QDSU |
| `_Schneeausfluss_` | snowcover_outflow | snowmod.cc | unsatzon.cc | 4 | P |
| `_Schneespeicher_Total_` | snow_storage_total | snowmod.cc | unsatzon.cc | 2 | SSTO |
| `_Snow_out_` | melt_from_snowcover | snowmod.cc | unsatzon.cc | 3 | QSNOW_OUT |
| `_SurfaceEvaporationResistance_` | SurfaceEvaporationResistance | (input grid) | unsatzon.cc | 1 | RSE |
| `_Temperatur_` | temperature | (input grid) | unsatzon.cc | 1 | T |
| `_aquifer_thickness_` | aquifer_thickness_ | gwflow.cc | unsatzon.cc | 1 | gw_thick |
| `_clay_depth_` | clay_depth | (input grid) | unsatzon.cc | 1 | CLAY |
| `_colmation_frozen_` | gw_1_colmation_frozen | (input grid) | unsatzon.cc | 1 | GW1COLMFROZEN |
| `_colmation_unfrozen_` | gw_1_colmation_unfrozen | (input grid) | unsatzon.cc | 1 | GW1COLMUNFROZEN |
| `_gw_kolmation_` | gw_kolmation_ | (input grid) | unsatzon.cc | 3 | GW1COLMATION |
| `_ice_firn_` | ice_firn | (input grid) | unsatzon.cc | 1 | ICE |
| `_irrigation_codes_` | irrigation_codes | (input grid) | unsatzon.cc | 1 | IRRIGCODE |
| `_kolmation_` | kolmation | (input grid) | unsatzon.cc | 1 | KOLM |
| `_maxFieldCapacity_mm_` | maxFieldCapacity_mm | (input grid) | unsatzon.cc | 4 | FCMAX |
| `_max_ponding_storage_` | max_ponding_storage | (input grid) | unsatzon.cc | 2 | FILLPOR, MAXPOND |
| `_minFieldCapacity_mm_` | minFieldCapacity_mm | (input grid) | unsatzon.cc | 2 | FCMIN |
| `_river_depth_` | river_depth | (input grid) | unsatzon.cc | 1 | RIVDEP |
| `_river_links_` | river_links | (input grid) | unsatzon.cc | 1 | RIVLNK |
| `_river_width_` | river_width | (input grid) | unsatzon.cc | 1 | RIVWIT |
| `_snowstorage_liquid_` | snowstorage_liquid | snowmod.cc | unsatzon.cc | 2 | SLIQ |
| `_snowstorage_solid_` | snowstorage_solid | snowmod.cc | unsatzon.cc | 2 | SSNO |
| `_surftemp_` | surftemp | snowmod.cc | unsatzon.cc | 1 | SURFTEMP |
| `_usableFieldCapacity_mm_` | usableFieldCapacity_mm | (input grid) | unsatzon.cc | 4 | FCACT |
| `_usableFieldCapacity_perc_` | usableFieldCapacity_perc | (input grid) | unsatzon.cc | 3 | FCACTPERC |

## Per module

### meteoint.cc (meteo interpolation + precipitation correction)

Creates: `_Niederschlagskorrektur_`

Uses from the grid table: `_CoordX_`, `_CoordY_`, `_MPIHelper_`, `_MPIRegions_`, `_Niederschlag_` (writes), `_RelCellArea_`, `_ScalingGrid_`, `_Temperatur_`, `_Windgeschwindigkeit_`

### expokorr.cc (radiation/temperature correction)

Creates: `_Expositionskorrektur_`, `_IEFF_MAX_`, `_Korrigierte_Temperatur_`, `_Nacht_Temperatur_`, `_Schatten_`

Uses from the grid table: `_Exposition_`, `_Gefaelle_`, `_Hoehenmodell_`, `_MPIHelper_`, `_MPIRegions_`

### etpetr.cc (evapotranspiration)

Creates: `_AspectCorrectedGlobalRadiation_`, `_Cloudiness_`, `_Potential_Evaporation_`, `_Potentielle_Verdunstung_`, `_Reale_Evaporation_`, `_Snow_Evaporation_`

Uses from the grid table: `_Albedo_`, `_Bodenspeicher_` (writes), `_GlacierizedCells_`, `_Globalstrahlung_` (writes), `_Hoehenmodell_`, `_LAKE_CODES_`, `_Landnutzung_`, `_Luftdruck_`, `_MPIHelper_`, `_MPIRegions_`, `_Nacht_Temperatur_` (writes), `_Nettostrahlung_` (writes), `_Nutzbare_Feldkapazitaet_`, `_Potentielle_Verdunstung_`, `_SnowEvaporation_onGlac_` (writes), `_SolidSnowStorage_onGlacier_` (writes), `_Sonnenscheindauer_`, `_SurfaceEvaporationResistance_`, `_Windgeschwindigkeit_`, `_ice_firn_`, `_ponding_storage_` (writes), `_ponding_storage_top_` (writes)

### snowmod.cc (snow and glaciers)

Creates: `_EnergyBalance_`, `_EnergyBalance_onGlac_`, `_GlacierCodesOld_`, `_GlacierMassBalance_`, `_GlacierizedCellsOld_`, `_LiquidSnowStorage_onGlacier_`, `_OldGlacierMassBalance_`, `_Regenanteil_`, `_Schneealter_`, `_Schneeanteil_`, `_Schneeausfluss_`, `_Schneespeicher_Total_`, `_Schneetage_`, `_SnowCoverDaysGrid_`, `_SnowEvaporation_onGlac_`, `_SnowFreeDaysGrid_`, `_SnowMelt_onGlacier_`, `_Snow_out_`, `_SolidSnowStorage_onGlacier_`, `_TotalSnowStorage_onGlacier_`, `_deposition_`, `_elev_order_col_`, `_elev_order_row_`, `_firn_WE_stack_`, `_input_mass_`, `_max_deposition_`, `_melt_from_firn_`, `_melt_from_ice_`, `_mobile_mass_`, `_snowstorage_liquid_`, `_snowstorage_solid_`, `_snowtemp_`, `_snowtemp_onGlac_`, `_surftemp_`, `_surftemp_onGlac_`

Uses from the grid table: `_Albedo_` (writes), `_AspectCorrectedGlobalRadiation_`, `_Cloudiness_`, `_Dampfdruck_`, `_EnergyBalance_` (writes), `_EnergyBalance_onGlac_` (writes), `_Gefaelle_`, `_GlacierCodes_` (writes), `_GlacierizedCells_` (writes), `_Korrigierte_Temperatur_`, `_LAKE_CODES_`, `_Landnutzung_`, `_Luftfeuchte_`, `_MPIHelper_`, `_MPIRegions_`, `_Nacht_Temperatur_`, `_Niederschlag_` (writes), `_PermafrostThawDepth_`, `_Schneealter_` (writes), `_Schneespeicher_Total_` (writes), `_Schneetage_` (writes), `_SnowCoverDaysGrid_` (writes), `_SnowDensity_` (writes), `_SnowFreeDaysGrid_` (writes), `_SnowLayerThickness_` (writes), `_SnowWaterEquivalentLiquid_` (writes), `_SnowWaterEquivalentSolid_` (writes), `_Snow_Canopy_IR_`, `_Snow_Canopy_RAD_`, `_Snow_Canopy_RH_`, `_Snow_Canopy_Temperature_`, `_Snow_Canopy_Unload_solid_`, `_Snow_Canopy_Wind_`, `_Snow_Evaporation_`, `_Soil_Temperature_Stack_Snow_` (writes), `_SurfaceTemperature_`, `_Temperatur_`, `_Temperatur_14_`, `_ThawCoeffPermaFrost_`, `_Windgeschwindigkeit_`, `_debris_on_glaciers_`, `_depoindex_`, `_ice_firn_` (writes), `_max_deposition_` (writes), `_ponding_storage_`, `_ponding_storage_top_`, `_slidefraction1_`, `_slidefraction2_`, `_slidefraction3_`, `_slidefraction4_`, `_snowstorage_liquid_` (writes), `_snowstorage_solid_` (writes), `_snowtemp_` (writes), `_snowtemp_onGlac_` (writes), `_surftemp_` (writes), `_surftemp_onGlac_` (writes), `_upperIndexSnow_` (writes)

### snow_canopy.cc (snow canopy interception)

Creates: -

Uses from the grid table: -

### interzep.cc (interception)

Creates: `_Canopy_Dummy_`, `_Canopy_Energy_Balance_`, `_Canopy_Temperature_`, `_Canopy_Temperature_Old_`, `_Interzeptionsspeicher_`, `_Snow_Canopy_Age_`, `_Snow_Canopy_Albedo_`, `_Snow_Canopy_IR_`, `_Snow_Canopy_Load_`, `_Snow_Canopy_RAD_`, `_Snow_Canopy_RH_`, `_Snow_Canopy_Sublim_`, `_Snow_Canopy_Temperature_`, `_Snow_Canopy_Unload_liquid_`, `_Snow_Canopy_Unload_solid_`, `_Snow_Canopy_Wind_`

Uses from the grid table: `_Albedo_`, `_Cloudiness_`, `_GlacierizedCells_`, `_Hoehenmodell_`, `_Interzeptionsspeicher_`, `_Landnutzung_` (writes), `_Luftfeuchte_`, `_MPIHelper_`, `_MPIRegions_`, `_Niederschlag_`, `_Schneeausfluss_`, `_Schneespeicher_Total_`, `_Windgeschwindigkeit_`, `_ice_firn_`, `_surftemp_`

### irrigate.cc (irrigation)

Creates: -

Uses from the grid table: `_ETR2ETP_`, `_MPIHelper_`, `_MPIRegions_`, `_Max_Pumping_`, `_Niederschlag_` (writes), `_SBHIGH_grid_`, `_SBLOW_grid_`, `_SB_1_grid_`, `_Unsatzonmodell_`, `_Wurzeltiefe_`, `_irrigation_amounts_` (writes), `_irrigation_codes_`

### siltingup.cc (silting up)

Creates: `_DIRECTRUNOFFSILTINGUP_`, `_KINETIC_ENERGY_`, `_RAINPAUSE_`, `_SILTING_UP_DISPOSITION_`, `_TIMESINCELASTSOILTILLAGE_`

Uses from the grid table: `_Bodenarten_`, `_Landnutzung_`, `_MPIHelper_`, `_MPIRegions_`, `_Niederschlag_` (writes), `_SURFSTORSILTINGUP_` (writes), `_Schneeausfluss_` (writes)

### surfrout.cc (surface routing)

Creates: `_CRITICAL_SHEAR_STRESS_`, `_PLANCURVATURE_`, `_SURFACEROUTING_FLUXES_`, `_SURFACEROUTING_HEAD_`, `_SURFACEROUTING_HEAD_OLD_`, `_SURFACEROUTING_MANNING_COEFF_`, `_SURFACEROUTING_QIN_`, `_SURFACEROUTING_QOUT_`, `_SURFACEROUTING_ROUGHNESS_X_`, `_SURFACEROUTING_ROUGHNESS_Y_`, `_SURFACEROUTING_VELOCITY_`, `_SURFACEROUTING_WEIGHTS_`

Uses from the grid table: `_Bodenarten_`, `_CRITICAL_SHEAR_STRESS_` (writes), `_CellSizeX_`, `_CellSizeY_`, `_Direktabfluss_` (writes), `_Fliesszeiten_`, `_Interflow_`, `_Landnutzung_`, `_MPIHelper_`, `_MPIRegions_` (writes), `_PLANCURVATURE_` (writes), `_SURFACEROUTING_QIN_` (writes), `_SURFACEROUTING_QOUT_` (writes), `_SURFSTORSILTINGUP_` (writes), `_river_depth_`, `_river_links_`, `_river_width_`

### unsatzon.cc (unsaturated zone (Richards))

Creates: `_AVGDD_`, `_AVINDAVG_`, `_AVINHLP_`, `_AVINMAX7DAVG_`, `_AVINSCAL_`, `_AVNMINTOT_`, `_AVNMIN_`, `_AVPGP_`, `_AVSNOWFREE_`, `_AVSON_`, `_AVTDAVG_`, `_AnnualTemperature_`, `_BALANCE_`, `_Basisabfluss_`, `_Capillar_Uprise_`, `_DailyTemperature_`, `_Direktabfluss_`, `_Drainageabfluss_`, `_Drainagelevel_`, `_ETR2ETP_`, `_ExcessIceFactor_`, `_FreezingDegreeDays_`, `_GWTableExtern_`, `_Grundwasserlevel_`, `_Grundwasserneubildung_`, `_Grundwasserstand_`, `_Infiltrationsexzess_`, `_Infiltrationssaettigung_`, `_Interflow_`, `_Landnutzung_`, `_Max_Pumping_`, `_Percolation_`, `_ROOTgrid_`, `_SBHIGH_grid_`, `_SBLOW_grid_`, `_SB_1_grid_`, `_SB_2_grid_`, `_SnowDensity_`, `_SnowLayerThickness_`, `_SnowWaterEquivalentLiquid_`, `_SnowWaterEquivalentSolid_`, `_Soil_SE_Stack_`, `_Soil_Temperature_Stack_`, `_Soil_Temperature_Stack_Snow_`, `_THROOT_grid_`, `_T_Lower_Boundary_Condition_`, `_Tangens_Beta_`, `_ThawDepthTMod_`, `_ThawingDegreeDays_`, `_Theta_`, `_Theta_Snow_`, `_Theta_in_GWLEVEL_`, `_Total_Infiltration_`, `_frozenSoilInGroundwater_`, `_geodetic_height_`, `_geodetic_height_Snow_`, `_gwflow_`, `_hydraulic_head_`, `_hydraulic_head_Snow_`, `_irrigation_amounts_`, `_macropore_infiltration_`, `_number_of_substeps_`, `_ponding_storage_`, `_river_exfiltration_`, `_river_infiltration_`, `_soil_bottom_fluxes_`, `_soil_bottom_fluxes_Snow_`, `_upperIndexSnow_`

Uses from the grid table: `_Albedo_` (writes), `_Bodenarten_` (writes), `_DIRECTRUNOFFSILTINGUP_` (writes), `_Drainagedistance_` (writes), `_Drainagehoehen_`, `_ETR2ETP_` (writes), `_EnergyBalance_` (writes), `_EnergyBalance_onGlac_`, `_FreezingDegreeDays_` (writes), `_Gefaelle_` (writes), `_GlacierizedCellsOld_` (writes), `_GlacierizedCells_` (writes), `_HeatFlux_Lower_Boundary_Condition_` (writes), `_IceWedgeLowerDepth_` (writes), `_IceWedgeUpperDepth_` (writes), `_LAKE_CODES_` (writes), `_Landnutzung_` (writes), `_MPIHelper_`, `_MPIRegions_` (writes), `_Niederschlag_` (writes), `_PermafrostThawDepth_` (writes), `_Potential_Evaporation_` (writes), `_Reale_Evaporation_` (writes), `_Regenanteil_` (writes), `_SB_1_REL_grid_` (writes), `_SURFSTORSILTINGUP_` (writes), `_Schneeausfluss_` (writes), `_Schneespeicher_Total_` (writes), `_Snow_out_` (writes), `_SurfaceEvaporationResistance_` (writes), `_T_Lower_Boundary_Condition_` (writes), `_Temperatur_` (writes), `_aquifer_thickness_` (writes), `_clay_depth_` (writes), `_colmation_frozen_` (writes), `_colmation_unfrozen_` (writes), `_gw_head_`, `_gw_kolmation_` (writes), `_gw_storage_coeff_`, `_gwflow_` (writes), `_ice_firn_` (writes), `_irrigation_codes_` (writes), `_kolmation_` (writes), `_maxFieldCapacity_mm_` (writes), `_max_ponding_storage_` (writes), `_minFieldCapacity_mm_` (writes), `_river_depth_` (writes), `_river_links_` (writes), `_river_width_` (writes), `_snowstorage_liquid_` (writes), `_snowstorage_solid_` (writes), `_surftemp_` (writes), `_surftemp_onGlac_`, `_usableFieldCapacity_mm_` (writes), `_usableFieldCapacity_perc_` (writes)

### gwflow.cc (groundwater)

Creates: `_aquifer_thickness_`

Uses from the grid table: `_Basisabfluss_` (writes), `_Bodenarten_`, `_CellSizeX_`, `_CellSizeY_`, `_GlacierizedCells_`, `_Grundwasserneubildung_`, `_Grundwasserstand_`, `_Hoehenmodell_`, `_LAKE_CODES_`, `_MPIHelper_`, `_MPIRegions_`, `_PermafrostThawDepth_`, `_SURFSTORSILTINGUP_`, `_ThawDepthTMod_`, `_aquifer_thickness_` (writes), `_frozenSoilInGroundwater_`, `_ice_firn_`, `_ponding_storage_`

### abfluss.cc (routing, lakes, reservoirs)

Creates: `_PrecOnLakes_`

Uses from the grid table: `_Basisabfluss_`, `_Hoehenmodell_`, `_LAKE_CODES_`, `_MPIHelper_`, `_MPIRegions_`, `_PrecOnLakes_`, `_Reale_Verdunstung_` (writes), `_Zonengrid_Bodenmodell_`, `_max_ponding_storage_`, `_ponding_storage_` (writes), `_ponding_storage_top_` (writes)

### dataassim.cc (data assimilation)

Creates: -

Uses from the grid table: `_MPIHelper_`, `_MPIRegions_`

### defltuzr.cc (initialization (zones, tables, MPI))

Creates: `_CellSizeX_`, `_CellSizeY_`, `_CoordX_`, `_CoordY_`, `_HydrologicResponseUnits_`, `_Landnutzung_`, `_MPIHalo_`, `_MPIHelper_`, `_MPIRegions_`, `_MasterHydrologicResponseUnits_`, `_RelCellArea_`, `_regression_regions_1_`, `_regression_regions_2_`, `_regression_regions_3_`, `_regression_weights_1_`, `_regression_weights_2_`, `_regression_weights_3_`

Uses from the grid table: `_CellSizeX_` (writes), `_CellSizeY_` (writes), `_HydrologicResponseUnits_` (writes), `_Landnutzung_` (writes), `_MPIHelper_` (writes), `_MPIRegions_` (writes), `_RelCellArea_` (writes), `_THROOT_grid_`, `_river_links_` (writes)

### modapuzr.cc (application (time loop, wiring))

Creates: -

Uses from the grid table: -

### teilmod.cc (submodel base, statistics)

Creates: -

Uses from the grid table: `_MPIHelper_`, `_MPIRegions_`, `_RelCellArea_`, `_Zonengrid_Bodenmodell_`
