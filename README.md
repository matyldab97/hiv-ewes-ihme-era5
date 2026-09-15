# People living with HIV exposure to extreme weather events in sub-Saharan Africa

Analysis code for estimating the number of people living with HIV (PLHIV) in
sub-Saharan Africa (SSA) exposed to heatwaves, extreme rainfall and drought in 2018,
combining ERA5 reanalysis with gridded HIV estimates from the Institute for
Health Metrics and Evaluation (IHME).

The pipeline produces:

- extreme weather indicators at 0.25° resolution for 2018, with
  thresholds and reference distributions defined over 1991–2020;
- 1000 Monte Carlo draws of gridded PLHIV counts propagating the published
  uncertainty intervals;
- the number and proportion of PLHIV in each exposure and severity category,
  continent-wide and by country;
- Getis-Ord Gi\* hotspots of HIV prevalence and of each indicator, and the
  overlap between them.

---

## Repository layout

```
00_inputs/        what the ERA5 input must be, and a check that it is
01_setup/         paths, configuration and the two shared libraries
02_indicators/    the extreme weather indicators
03_population/    Monte Carlo simulation of gridded PLHIV
04_exposure/      indicator screening, severity classification, threshold exposure
05_hotspots/      Gi* inputs, transformation choice, Gi*, overlap analysis
06_supplementary/ optional analyses that feed the supplements only
```

Every script is run directly, in the order given below, except the four in
`01_setup/`, which are sourced by the others rather than run.

---

## Requirements

R ≥ 4.1, with GDAL, GEOS and PROJ available to `terra` and `sf`.

```r
install.packages(c(
  "terra", "sf", "spdep", "exactextractr", "ncdf4",
  "SPEI", "zoo", "moments", "Matrix",
  "dplyr", "tidyr",
  "ggplot2", "scales", "gridExtra", "png", "openxlsx"
))
```

| Package | Used for | Minimum |
| --- | --- | --- |
| `terra` | all raster work | 1.7 |
| `sf` | boundary polygons and the nearest-country assignment | 1.0 |
| `spdep` | the neighbour graph and the Gi\* statistic | 1.2 |
| `exactextractr` | area-weighted aggregation from the fine population grid | 0.9 |
| `ncdf4` | required by `terra::writeCDF()` in step 02, not loaded by name | -- |
| `SPEI` | the gamma and log-logistic fits behind SPI and SPEI | 1.8 |
| `zoo` | rolling accumulation windows | -- |
| `moments` | skewness and kurtosis in step 17 | -- |
| `Matrix` | the sparse weight matrix in step 12 | -- |
| `dplyr`, `tidyr` | table assembly | -- |
| `ggplot2`, `scales` | figures | -- |
| `gridExtra`, `png` | the multi-panel screening pages in step 13 | -- |
| `openxlsx` | workbook output | -- |

`grid`, `tools` and `utils` ship with R and need no installation.

Obtaining the ERA5 input additionally needs a Copernicus Climate Data Store
account. The download and the reduction from hourly to daily fields are done
once, outside this repository.

### Project root

All paths resolve from a single root directory, `EWE_ROOT`, so that no script
contains an absolute path. Set it before running anything, either with
`Sys.setenv(EWE_ROOT = "/path/to/project")` or by creating an untracked file
`01_setup/00_paths_local.R` that contains `EWE_ROOT <- "/path/to/project"`. Every
script is run from the repository root, for example
`Rscript 02_indicators/07_index_spi.R`.

### Expected layout under `EWE_ROOT`

```
data/era5/<var>/            daily ERA5 NetCDF, one folder per variable
data/plhiv/                 IHME gridded PLHIV counts, fine resolution
data/prevalence/            IHME gridded HIV prevalence, fine resolution
data/boundaries/gadm.gpkg   GADM polygons, subset to the 43 SSA countries
results/indices/            extreme weather indicator rasters
results/montecarlo/         PLHIV draws and uncertainty intervals
results/exposure/figures/   exposure figures, JPEG and PDF
results/exposure/tables/    exposure workbooks and CSV
results/screening/          indicator screening workbooks
results/screening/EWE_verification_figures/   one page per candidate indicator
results/gistar/inputs/      Gi* input stacks (augmented, ready, reports, transformation)
results/gistar/outputs/     Gi* rasters and CSVs
results/gistar/plots/       Gi* maps, by reference frame and by topic
tmp/                        scratch space for terra
```

`results/` and `tmp/` are created automatically.

---

## Run order

Steps 01–11 build the indicators. Step 12 is independent of them and may be run
in parallel. Steps 13 onward require both.

| # | Script | Requires | Produces |
| --- | --- | --- | --- |
| 01 | `00_inputs/01_era5_inputs.R` | daily ERA5 NetCDF under `data/era5/` | verification that the ERA5 input is complete and in raw units |
| -- | `01_setup/00_paths.R` | -- | path constants (sourced, not run) |
| -- | `01_setup/01_config.R` | 00 | variable table, `load_year()` (sourced) |
| -- | `01_setup/02_lib_indices.R` | 01 | drought index engine (sourced) |
| -- | `01_setup/03_lib_shared.R` | 01 | helpers shared by steps 12 to 19 (sourced) |
| 02 | `02_indicators/02_build_soil_moisture.R` | swvl1–3 | depth-weighted soil moisture |
| 03 | `02_indicators/03_temperature_thresholds.R` | tmax | TX90p and TX95p day-of-year thresholds |
| 04 | `02_indicators/04_heatwaves.R` | 03 | HWN, HWF, HWM, HWMF, HWNM |
| 05 | `02_indicators/05_extreme_rainfall.R` | pr | R95p and R99p thresholds, day counts and totals |
| 06 | `02_indicators/06_severity_excess_percentile.R` | 03, 05 | HWES, HWPD, rainfall excess and percentile deviation |
| 07 | `02_indicators/07_index_spi.R` | pr | SPI at 1 and 3 month accumulation; unfittable cell-months set to zero |
| 08 | `02_indicators/08_index_spei.R` | pr, pet | SPEI at the same scales |
| 09 | `02_indicators/09_index_sri.R` | ro | SRI at the same scales; unfittable cells set to zero |
| 10 | `02_indicators/10_index_sma.R` | 02 | SMA at the same scales |
| 11 | `02_indicators/11_drought_month_count.R` | 07–10 | annual drought-month counts per index and threshold |
| 12 | `03_population/12_monte_carlo_plhiv.R` | IHME, GADM | 1000 PLHIV draws, uncertainty intervals, country layer; written for females, males and the two combined |
| 13 | `04_exposure/13_indicator_screening.R` | 02–12 | one page set per candidate indicator, with distribution summaries |
| 14 | `04_exposure/14_exposure_classification.R` | 02–12 | PLHIV by severity category, per indicator |
| 15 | `04_exposure/15_sevenday_exposure.R` | 02–12 | PLHIV exposed to ≥7 days or ≥1 month per hazard |
| 16 | `05_hotspots/16_gistar_inputs.R` | 02–12 | transformed Gi\* input stack |
| 17 | `05_hotspots/17_transformation_selection.R` | 16, stage 1 | transformation selection grid |
| 18 | `05_hotspots/18_gistar.R` | 16 | Gi\* z-scores, BH-adjusted p-values, categories |
| 19 | `05_hotspots/19_gistar_postprocessing.R` | 18, 12 | hotspot overlap maps, PLHIV in overlaps |
| 20 | `06_supplementary/20_day_threshold_sweep.R` | 04, 05, 12 | supplementary: exposure at every day threshold from 1 to 14 |

The analysis covers the combined population aged 15 to 59. Step 12 models the
eighteen age and sex groups internally and aggregates them before writing, so a
sex-stratified analysis would need step 12 to write the groups separately.

### Figures

Every figure the pipeline writes is reported in the manuscript or the
supplements. Steps 14 and 15 write six between them:

| File | Shows |
| --- | --- |
| `exposure_severity_by_country_main` | share of each country's PLHIV in each severity class, three main indicators, all 43 countries |
| `exposure_severity_by_country_sensitivity` | the same for the six sensitivity indicators |
| `exposure_very_high_map_main` | share of each country's PLHIV in the Very High class, as a map, three main indicators |
| `exposure_very_high_map_sensitivity` | the same for the six sensitivity indicators |
| `sevenday_exposure_map_main` | share of each country's PLHIV exposed for at least seven days, three main indicators |
| `sevenday_exposure_map_sensitivity` | the same for the six sensitivity indicators |

The first four come from step 14 and the last two from step 15, which writes
each map as JPEG and as PDF. Step 19 additionally writes the hotspot overlap
panels described under Gi* plot folders below.

### Gi* plot folders

| Folder | Holds |
| --- | --- |
| `plots/SSA/` | Gi* z scores and hotspot classes, thresholds taken across the whole of SSA |
| `plots/COUNTRY/` | the same maps with thresholds taken within each country |
| `plots/prevalence_estimates/` | the Gi* result repeated for the three HIV prevalence surfaces: the mean and the lower and upper bounds of its interval |
| `plots/overlap_vulnerability/main/` | where an EWE hotspot meets an HIV prevalence hotspot, for the three main indicators |
| `plots/overlap_vulnerability/robustness/` | the same for the six sensitivity indicators, and at the 90, 95 and 99 per cent thresholds |

---

## Conventions

**Indicator naming.** Directories and file names use the index name alone:
`SPI`, `SPEI`, `SRI`, `SMA`, `TX90p`, `TX95p`, `R95p`, `R99p`, `Heatwaves`.
The soil moisture index is named `SMA` throughout, matching the manuscript.

**Masking.** Indicator rasters are written unmasked. The analysis domain is
applied downstream, in steps 14–16, using the modelled PLHIV footprint taken
directly from the Monte Carlo output. Applying a coarser mask earlier, by
nearest-neighbour resampling, drops coastal cells.

**Sign convention.** Drought severity sums are stored as signed (negative)
values and multiplied by −1 wherever they enter a classification or a Gi\*
input, so that higher values denote more severe conditions for every indicator.

**Point estimates.** All point estimates are the mean across the 1000 draws;
uncertainty intervals are the 2.5th and 97.5th percentiles.
