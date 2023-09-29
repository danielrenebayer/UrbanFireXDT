# Output description

The simulation outputs multiple files, that are described below.
Folders are marked in **bold**.

| File name and level                                 | Description |
| ---                                                 | -----       |
| `expansion-matrix-abs-values.csv`                   | Expansion matrix filled with absolute numbers      |
| `expansion-per-cu.csv`                              | List the expansion that happens per control unit   |
|                                                     |                                                    |
| **`no param vari`**                                 | Output directory if no param variation is selected |
| &nbsp;&nbsp;\| `XXXX-substation-time-series.csv`    | Time series of the load, wind/o.s. PV feedin, etc. aggregated on substation level |
| &nbsp;&nbsp;\| `metrics-per-cu.csv`                 | Computed metrics per control unit (like SSR, SCR, NPV, total grid demand, ...) including the concrete parameters for the simulated PV, battery and the charging station component |
| &nbsp;&nbsp;\| `parameter-settings.csv`             | Parameter settings for the concrete simulation run (very useful for parameter variations) |
| &nbsp;&nbsp;\| `sim-added-roof-sections-per-cu.csv` | A list of all sim. added roof sections per control unit |
| &nbsp;&nbsp;\| `substation-detailed-time-series.csv`| Addition information, aggregated on substation level, that is not contained in the first file |
| &nbsp;&nbsp;\| `build_and_run_info.txt`             | Information on the simulation run and the used program  |
|                                                     |                                                         |
| **`param vari XXXX`**                               | Output directory if a parameter variation is selected   |
| &nbsp;&nbsp;\| `build_and_run_info.txt`             | Information on the simulation run and the used program  |
| &nbsp;&nbsp;\| **`variation index XXXX`**           | There is one subfolder for every variation combination, the concrete parameters can be found in `parameter-settings.csv` |
| &nbsp;&nbsp;\|&nbsp;&nbsp;\| *rest of the content*  | see folder `no param vari` | |

