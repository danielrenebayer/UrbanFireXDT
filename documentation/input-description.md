% Input description



# General overview of the input data

- The simulation requires a sqlite3-database with the structure of the system that should be simulated (default filename is SystemStructure.db).
- The existing smart meter data of the individual measurement units must be stored in a subfolder called 'SeparatedSmartMeterData' as individual comma-seperated values (csv) files.
- The file name of the CSV files must equal the measurement unit ID (key: `list_of_measurement_units.MeUID`) as they are defined in the system structure database.


## The general structure of the system structure database

The system structure database requires the following tables:

| Table name                | Description               |
| ---                       | --------                  |
| time_indices              | List of time indices      |
| list_of_substations       | List of substations       |
| list_of_control_units     | List of control units     |
| list_of_measurement_units | List of measurement units |
| global_profiles_pv        | Global PV feedin (normalized) |
| global_profiles_pv_info   | Information about how many global PV feedin time series exist per orientation |
| global_profile_wind       | Global wind farm feedin (normalized) |
| global_profiles_heatpumps | Global heap pump demand profiles (normalized) |
| address_data              | List of all known addresses with yearly energy demand in case of heat pump heating |
| heat_demand_per_location  | Information about the heat demand that is converted using gas consumption data (if available) |
| address_roof_data         | List of all roof sections of all buildings (with orientation and size) |
| residual_grid_load        | Time series of the residual grid load, that is not measured by any smart meter, but still present |
| electricity_emissions     | Time series of the CO2 emissions (optional table, is used, if present, otherwise ignored) |
| electricity_prices        | Time series of the hourly electricity prices (optional table, is used, if present) |


## Structure of the smart meter data / csv files

Each smart meter csv file must have the following columns in this order:

| Column name   | Description |
| ---           | -------     |
| TimestepID    | The id of the time step (key: `time_indices.TimestepID`) |
| Value_Demand  | Reading of the demand in kWh at the given time step |
| Status_Demand | *ignored*   |
| Value_Feedin  | Reading of the feed-in in kWh at the given time step |
| Status_Feedin | *ignored*   |





# Details on the database tables

Documentation of the configuration database in detail


## Table time_indices

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | ID of the timestep, must start with 1 and must be numbered consecutively |
| UTC_time      | TIMESTAMP | UTC time stamp of the beginning of the time step in the format YYYY-MM-DD HH:MM:SS |
| local_time    | TIMESTAMP | Local time stamp of the beginning of the time step in the format YYYY-MM-DD HH:MM:SS |
| local_time_zone | TEXT/VARCHAR(4) | Time zone name of the local time as text (might change e.g. through daylight saving time) |


## Table list_of_substations

| Column name   | Type     | Description       |
| ---           | ---      | -------           |
| substation_id | INTEGER  | ID of the substation, may start with 1 and may be numbered consecutively (not mandatory anymore) |
| substation_name | TEXT   | Name of the substation |


## Table list_of_control_units

| Column name   | Type     | Description       |
| ---           | ---      | -------           |
| UnitID        | INTEGER  | ID of the control unit, may start with 1 and may be numbered consecutively (not mandatory anymore) |
| substation_id | INTEGER  | ID of the substation to which this control unit is connected to [Foreign key] |
| LocID         | INTEGER  | ID of the location where the control unit is located at [Foreign key]  |
| has_cs        | INTEGER  | Holds the value 1, iff there is a charging station already connected to this control unit
| n_flats       | INTEGER  | The number of flats in a building / for a location (estimation or real value) |


## Table list_of_measurement_units

| Column name   | Type     | Description       |
| ---           | ---      | -------           |
| MeUID         | INTEGER  | ID of the measurement unit, may start with 1 and may be numbered consecutively (not mandatory anymore) |
| UnitID        | INTEGER  | ID of the control unit to which this measruement unit is assigned [Foreign key] |
| UnitID        | INTEGER  | ID of the control unit to which this measruement unit is assigned [Foreign key] |
| MeterPointID  | TEXT     | Name of the meter point, i.e. MPRN (Meter Point Reference Number in UK) or MELO (Meter Location Number in Germany) |
| has_demand    | INTEGER  | Holds the value 1, iff this measurement unit shows a demand at least at one point in the simulation time |
| has_feedin    | INTEGER  | Holds the value 1, iff this measurement unit shows a feedin at least at one point in the simulation time |
| has_pv_residential | INTEGER | Holds the value 1, iff a residential PV installation is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_pv_open_space  | INTEGER | Holds the value 1, iff an open-space PV installation is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_bess      | INTEGER | Holds the value 1, iff a battery is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_hp        | INTEGER | Holds the value 1, iff a heat pump is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_chp       | INTEGER | Holds the value 1, iff a CHP is connected (exclusivley and not exclusivley) to the measuremt unit |
| LocID         | INTEGER | ID of the location where the measurement unit is located at [Foreign key] |
| has_wind      | INTEGER | Holds the value 1, iff a wind farm is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_biomass   | INTEGER | Holds the value 1, iff a bio mass power plant is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_evcs      | INTEGER | Holds the value 1, iff a residential EV charging station is connected (exclusivley and not exclusivley) to the measuremt unit |
| has_public_evcs | INTEGER | Holds the value 1, iff a public EV charging station is connected (exclusivley and not exclusivley) to the measuremt unit |


## Table global_profiles_pv

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | The time step for which the dataset is valid [Foreign key] |
| Value_Feedin  | REAL      | The normalized feed-in value at the given timestep |
| Orientation   | TEXT / VARCHAR(2) | The orientation of the given feedin value |
| SameOrientationTimeSeriesIndex | INTEGER | If there are more time series for one orientation, this number gives the id (starting with 0) of the time series with the same direction |


## Table global_profiles_pv_info

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| orientation   | TEXT / VARCHAR(2) | The orientation |
| number_of_ts  | INTEGER   |Number of time series for the given orientation |


## Table global_profile_wind
| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | The time step for which the dataset is valid [Foreign key] |
| wind_profile_value |REAL  | The normalized feed-in value at the given timestep |


## Table global_profiles_heatpumps

This table holds the normalized demand values of every heat pump time series in kW.
The normalization takes place on an annual level, so that every time series results in a annual demand sum of 1000 kWh / a.

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | The time step for which the dataset is valid [Foreign key] |
| ShiftableDemand_kW   | REAL | The normalized demand value of the heat pump demand at this time step that can be shifted in theory |
| UnshiftableDemand_kW | REAL | The normalized demand value of the heat pump demand that cannot be shifted |
| TimeSeriesIndex | INTEGER | The index of the heat pump time series (important, if there is more than one time series) |


## Table address_data

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| LocID         | INTEGER   | ID of the location (starting with 0) |
| n_buildings   | INTEGER   | Number of buildings (for which geodata is available) |
| has_residential_buildings | INTEGER | Holds the value 1, iff the location holds a residential building, and can thus be regarded as residential location |
| max_volume    | REAL      | The volume of the (biggest) building on the given location - used for aproximating the heat demand if no heat demand is given in table heat_demand_per_location |


## Table heat_demand_per_location

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| LocID         | INTEGER   | ID of the location |
| MeanHeatEnergy_kWh | REAL | Mean heat demand of the location in kWh per year / mean over all regarded years |


## Table address_roof_data

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| LocID         | INTEGER   | ID of the referenced location [Foreign key] |
| Area_in_m2    | REAL      | Area of the roof section in square meter |
| Orientation   | TEXT / VARCHAR(2) | The orientation of the given roof section; for every orientation there has to exist a correspoding time series in the tabel `global_profiles_pv` |


## Table residual_grid_load

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | ID of the timestep, must start with 1 and must be numbered consecutively |
| P_residual_gridload | REAL | Residual grid load in kW during time step TimestepID |


## Table electricity_emissions

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | ID of the timestep, must have the same alignment as time_indices.TimestepID |
| emissions_g_kWh | REAL    | CO2eq. emissions of one kWh of electiricty that is taken from the transmission grid |


## Table electricity_prices

| Column name   | Type      | Description       |
| ---           | ---       | -------           |
| TimestepID    | INTEGER   | ID of the timestep, must have the same alignment as time_indices.TimestepID |
| local_price   | REAL      | The local (dynamic) price for residential households for a dynamic tariff in ct/kWh |
| spotmarket_price | REAL   | The (dynamic) price at the spot market (without vat or taxes) in ct/kWh |
