% Output description

# General description of the output files

The simulation outputs multiple files, that are described below.
Folders are marked in **bold**.

| File name and level                                 | Description |
| ---                                                 | -----       |
| `expansion-matrix-abs-values.csv`                   | Expansion matrix filled with absolute numbers      |
| `expansion-per-cu.csv`                              | List the expansion that happens per control unit   |
|                                                     |                                                    |
| **`no param vari`**                                 | Output directory if no param variation is selected |
| &nbsp;&nbsp;\| `metrics-per-cu.csv`                 | Computed metrics per control unit (like SSR, SCR, NPV, total grid demand, ...) including the concrete parameters for the simulated PV, battery and the charging station component |
| &nbsp;&nbsp;\| `weekly-metrics-per-cu.csv`          | Computed metrics per control unit on a weekly level. See `metrics-per-cu.csv` for details. |
| &nbsp;&nbsp;\| `XXXX-substation-time-series.csv`    | Time series of the load, wind/o.s. PV feedin, etc. aggregated on substation level |
| &nbsp;&nbsp;\| `substation-detailed-time-series.csv`| Time series of the additional details on substation level, like residential load (includes feed-in) and demand (only residential grid demand, no surplus feed-in accumulated) |
| &nbsp;&nbsp;\| `parameter-settings.csv`             | Parameter settings for the concrete simulation run (very useful for parameter variations) |
| &nbsp;&nbsp;\| `sim-added-roof-sections-per-cu.csv` | A list of all sim. added roof sections per control unit |
| &nbsp;&nbsp;\| `substation-detailed-time-series.csv`| Addition information, aggregated on substation level, that is not contained in the first file |
| &nbsp;&nbsp;\| `build_and_run_info.txt`             | Information on the simulation run and the used program  |
| &nbsp;&nbsp;\| `metrics-per-ev.csv`                 | Computed metrics per EV                                 |
| **`param vari XXXX`**                               | Output directory if a parameter variation is selected   |
| &nbsp;&nbsp;\| `build_and_run_info.txt`             | Information on the simulation run and the used program  |
| &nbsp;&nbsp;\| **`variation index XXXX`**           | There is one subfolder for every variation combination, the concrete parameters can be found in `parameter-settings.csv` |
| &nbsp;&nbsp;\|&nbsp;&nbsp;\| *rest of the content*  | see folder `no param vari` |
|                                                     |                            |
| _Optional output per cu_:                           | _Files only generated if selected by specified command line arguments. See command line parameter `--cu-output` for details_ |
| &nbsp;&nbsp;\| `XXXX-CU-time-series.csv`            | Output of all power flows inside the control units **in one file** |
| &nbsp;&nbsp;\| **`ts-per-cu`**                      | Output of all power flows inside the control units **in separated files (grouped by substation)** |
| &nbsp;&nbsp;\|&nbsp;&nbsp;\| `YYYY-AllCUs-ts.csv`   |                            |
|                                                     |                            |
| _Further optional output_:                          | _Files only generated if selected by specified command line arguments. See command line parameter `--cu-output` for details_ |
| ev-details.csv                                      | Detailed time series for every simulated EV |



# Detailed description per file

## metrics-per-cu.csv and weekly-metrics-per-cu.csv

| Column name | Description |
| ---         | ------      |
| UnitID      | The control unit ID of the current record |
| SCR         | The self-consumption rate over the complete simulated time span |
| SSR         | The self-sufficiency rate over the complete simulated time span |
| NPV         | The net-present value (including installation costs) over the complete simulated time span [not in weekly output] |
| ALR         | Array-to-load ratio (see Nyholm et al., 2016 in Applied Energy, https://doi.org/10.1016/j.apenergy.2016.08.172) [not in weekly output] |
| BDR         | Battery-to-Demand ratio (see Nyholm et al., 2016) [not in weekly output]  |
| RBC         | Relative battery capacity (see Nyholm et al., 2016) [not in weekly output] |
| Sum of demand [kWh]              | The demand sum of electricity demand of the real smart meters and all simulated consumers (heat pump and charging station) in kWh |
| Sum of MU demand [kWh]           | The demand sum of the measurement units, i.e., the real smart meters in kWh |
| Sum of self-consumed e. [kWh]    | The sum of locally produced and self-consumed electricity in kWh |
| Sum of PV-generated e. [kWh]     | The sum of locally produced PV energy in kWh |
| Sum of grid feed-in [kWh]        | The sum of electricity that was fed into the grid in kWh |
| Sum of grid demand [kWh]         | The sum of electricity that was demanded from the grid in kWh |
| BS EFC                           | The battery equivalent full cycles (if not battery is simulated, this value defaults to 0.0) |
| BS n\_ts\_empty                  | The number of time steps where the battery was empty (if not battery is simulated, this value defaults to 0.0) |
| BS n\_ts\_full                   | The number of time steps where the battery was full (if not battery is simulated, this value defaults to 0.0) |
| BS total E withdrawn [kWh]       | The sum of electricity that was withdrawn from the battery in kWh (if not battery is simulated, this value defaults to 0.0) - Attention: The battery might not be empty at the end of the simulation procedure |
| Sum of HP demand [kWh]           | The sum of demanded electricity of the heat pump in kWh (if no heat pump is simulated, this value defaults to 0.0) |
| Sum of CS demand [kWh]           | The sum of demanded electricity of the EV charging station in kWh (if no EV charging station is simulated, this value defaults to 0.0) |
| Peak grid demand [kW]            | The peak grid demand in kW over the complete simulated time span |
| Emissions cbgd [kg CO2eq]        | The CO2e emissions caused by grid demand including upstream emissions (but no installation emissions for PV panels, battery storage or heat pumps) |
| Avoided emissions [kg CO2eq]     | The avioded CO2e emissions caused by local production (and possibly storing) of renewable energy |
| Electricity cons. costs [CU]     | The total electricity costs for consumption from the grid over the complete simulation horizon in the general currency unit (CU) |
| Avoided electricity cons. costs [CU] | The sum of electricity costs that were saved by local self-consumption over the complete simulation horizon |
| Feed-in revenue [CU]             | The summed revenue of feed-in to the grid in the general currency unit (CU) |
| Sim. PV max P [kWp]              | The installed power in kW (resp. kWp) of the simulated PV installation, or 0.0 if none is present |
| Sim. BS P [kW]                   | The installed power in kW of the simulated battery storage system, or 0.0 if none is present |
| Sim. BS E [kWh]                  | The installed capacity in kWh of the simulated battery storage system, or 0.0 if none is present |
| n EVs                            | The number of EVs with their home at this given unit |
| Sim. CS max P [kW]               | The installed power in kW of the simulated EV charging station, or 0.0 if none is present |
| Simulated PV                     | Indicates whether this control unit contains a simulated PV installation     |
| Simulated BS                     | Indicates whether this control unit contains a simulated battery storage     |
| Simulated HP                     | Indicates whether this control unit contains a simulated heat pump           |
| Simulated CS                     | Indicates whether this control unit contains a simulated EV charging station |
| n errors in cntrl                | Number of error happend in controller optimization (e.g., infeasible optimization problem) |
| n errors cntrl cmd appl          | Number of error happend during the application of the control commands in at least one component (e.g., bound violations in heat pump or EV charging station / EVFSM) |

In the case of weekly metrics file, an additional column giving the week number (named `Week number`) is added.

## XXXX-substation-time-series.csv

| Column name           | Description |
| ---                   | ------      |
| Timestep              | The time step ID of the current record |
| _Per substation YYYY_ |             |
| &nbsp;&nbsp; YYYY     | Active power in kW at substation YYYY |
| pv_feedin_kW          | Sum of feed-in from control units caused by PV installations in kW |
| bs_feedin_kW          | Sum of feed-in from control units caused by battery storage systems in kW |
| chp_feedin_kW         | Sum of feed-in from control units caused by CHPs in kW |
| wind_feedin_kW        | Sum of feed-in from control units caused by wind turbines in kW |
| unknown_feedin_kW     | Sum of feed-in from control units in kW where the source is unknown or not attributable to any generation technology |
| OverallBatterySOC     | Mean SOC over all simulated battery storage systems on control unit level |
| total_load            | Total active power in kW in the grid |


## substation-detailed-time-series.csv

The columns of this output follows the following nomenclature:

- **load:** Sum of all current virtual smart meter measurements per time step 
- **demand:** Sum of the positive current virtual smart meter measurements per time step

| Column name             | Description |
| ---                     | ------      |
| Timestep                | The time step ID of the current record |
| _Per substation YYYY_   |             |
| &nbsp;&nbsp; YYYY_resident_load_kW   | Sum of residential load in kW of all residential buildings connected to substation YYYY |
| &nbsp;&nbsp; YYYY_resident_demand_kW | Sum of residential demand in kW of all residential buildings connected to substation YYYY |
| total_residential_load   | Total residential load in kW   |
| total_residential_demand | Total residential demand in kW |


## XXXX-CU-time-series.csv or YYYY-AllCUs-ts.csv

| Column name           | Description |
| ---                   | ------      |
| Timestep              | The time step ID of the current record |
| ControlUnitID         | The control unit ID of the current record |
| Load_vSmartMeter_kW   | The power at the virtual smart meter in kW (positive values denote a demand, negative values denote a feed-in) |
| Load_rSmartMeters_kW  | The power summed over all real smart meters in kW (positive values denote a demand, negative values denote a feed-in) |
| Load_self_produced_kW | The self-produced power that is consumed in the current time step in kW at the given unit - Please note: If excess power was fed into the battery in the previous steps, this will only be taken into account at the moment of discharge |
| PVFeedin_simulated_kW | The current PV production in kW (if present, else 0.0) |
| BS_SOC                | The state of charge of the connected battery (if present, else 0.0) |
| BS_load_kW            | The current power of the battery storage in kW (if present, else 0.0; positive values denote battery charging, negative values denote discharging) |
| HP_load_kW            | The current power of the heat pump in kW (if present, else 0.0) |
| CS_load_kW            | The current power of the EV charging station in kW (if present, else 0.0) |
| CS_n_EVs_conn         | The number of home-parking vehicles that are currently connected to the charging station |
| CS_n_EVs_not_conn     | The number of home-parking vehicles that are currently not connected to the charging station |


## metrics-per-ev.csv


| Column name                | Description |
| ---                        | ------      |
| CarID                      | The ID of the car as defined in the input data |
| Driving distance [km]      | The total distance the EV is driven in the simulation in km |
| E used for driving [kWh]   | The sum of electric energy required by the car in the simulation in kWh |
| Home-charged E [kWh]       | The sum of electric energy that has been charged at home in kWh. Can be higher than `E used for driving [kWh]`, especially in the case of bidirectional charging |
| Home-discharged E [kWh]    | The sum of electric energy that has been discharged from the EV battery at home in kWh (only important for the case of bidirectional charging) |

