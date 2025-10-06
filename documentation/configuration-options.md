% Description of the JSON configuration file

---
fontsize: 9pt
---

# General structure

The simulation requires a configuration file formatted as JSON.
A config file can have several scenario definitions.
The scenario to be simulated is selected via a command line argument when starting the simulation.

The JSON file holds a dictionary with three keys / sections on the top level:

| Number | Section / key name  | Description |
| ---    | ---                 | ----------- |
|  1     | Default Scenario Values | A general section defining default values for the configuration variables that apply to all scenarios (if not overwritten) |
|  2     | Scenarios           | A list of scenarios, where every scenario is a dictionary with key-value entries for defining the configuration variables for this scenario, only |
|  3     | Parameter Variation | A list of possible parameter variation that can be combined with every scenario |

```{.json .numberLines}
{
    "Default Scenario Values": [],
    "Scenarios": [
        {
            "id": 1,
            ...
        }
    ],
    "Parameter Variation": []
}
```

# Details on the available configuration variables


## Group 1: General scenario definitions

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| id                    | integer | yes      |                      | 0 - inf         | ID of the defined scenario - must be unique in the present file | no |
| inherits from         | integer | no       |                      | 0 - inf         | Define another scenario ID from which the current scenario should inherit the defined values | no |
| expansion id          | integer | yes      |                      |                 | ID of the scenario which defines the expansion | no |
| start                 | string  | yes      |                  | YYYY-mm-dd HH:MM:SS | Start date in the format YYYY-mm-dd HH:MM:SS - if the data start after start date, simulation will wait for the data to start | no |
| end                   | string  | yes      |                  | YYYY-mm-dd HH:MM:SS | End date in the format YYYY-mm-dd HH:MM:SS - if the data end earlier, simulation will stop there | no |
| tariff feed-in per kWh| float   | yes      |                      |                 |             | no |
| tariff demand per kWh | float   | yes      |                      |                 | The demand per kWh energy taken from the grid - This parameter is ignored if a time series on the demand prices is given | no |
| emissions per kWh     | float   | no       | 100                  | 0.0 - inf       | The emissions in CO2eq per kWh of energy taken from the grid - This parameter is ignored if a time series on emissions is given | no |
| net present value discount rate | float | yes |                   |                 | The discount rate for the computation of the metric `NPV` in the metrics-per-cu.csv. If `expansion PV sizing mode` or `expansion BS capacity computation mode` is set to `Optimized`, this value is also used there. | no |
| net present value time horizon in years | float | yes |           |                 | The time horizon for the NPV calculation. See also the remarks for parameter `net present value discount rate`. | no |
| installation cost PV per kWp  | float | yes |                     |                 | Installation costs for a simulated PV installation per peak-power. Required for the NPV and EAC calculation. | no |
| installation cost BS per kWh  | float | yes |                     |                 | Installation costs for a simulated BS installation per kWh of capacity. Required for the NPV and EAC calculation. | no |
| installation cost HP per kW   | float | no  |                     |                 | Installation costs for a simulated HP installation per power. Required for the NPV and EAC calculation. | no |
| installation cost CS per unit | float | no  |                     |                 | Installation costs for a simulated charging station installation per charging unit. Required for the NPV and EAC calculation. Mind: A control unit might have multiple chargers per charging station. | no |

## Group 2: Data sources and output definition

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| time steps per hour   | integer | yes      |                      | 1 to infinity   | Number of time steps per hour - **must align with the input data!** | no |
| data input path       | string  | yes      |                      |                 | Path (relative or absolute) to the input data | no |
| data output path      | string  | yes      |                      |                 | Path (relative or absolute) to the directory where the output should be saved | no |
| database name         | string  | no       | "SystemStructure.db" |                 | Name of the database that contains the system structure | no |
| ev data path          | string  | no       |                      |                 | Path (relative or absolute) to the EV driving profiles and other EV related data | no |
| use emission time series ia | bool | no    | true                 |                 | Use the emission time series for electricity demand from the grid (if it is available in the data) if set to true (default) | no |
| use prices time series ia   | bool | no    | true                 |                 | Use the electricity prices time series (if it is available in the data) if set to true (default) | no |


## Group 3: Settings defining the simulatively added components

### Group 3, Part A: Rooftop PV options

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| expansion PV sizing mode              | string | no  | `MaxAvailableRoofArea` | - `MaxAvailableRoofArea`<br>- `StaticPVSize`<br>- `Optimized` | Defines the actual sizing mode for new PV installations.<br>- `MaxAvailableRoofArea` (default) uses the maximum available roof area given the restictions of the next variables <br>- `StaticPVSize` adds a static PV installation to all expanded control units with equal size and the same orientation <br>- `Optimized` uses a LP to determine the optimal PV size respecting the individual profile and the orientation of the roofs. This option only works if the `controller mode` is set to `opti with perfect forecast` and the `control horizon in ts` is set to the complete simulation time span and `control update freq in ts` is also set to this value. Finally, the option `controller optimization target` must be set to `electricity costs`. | no |
| expansion PV min kWp for section usage | float | no  |  0      |  | Minimal size of a PV installation on a roof section, to include this section | yes |
| expansion PV max inst kWp per section  | float | no  | -1      |  | Maximal size of a PV installation in kW per roof section | yes |
| expansion PV max inst kWp per unit     | float | no  | -1      |  | Maximal size of the total PV installation in kW per control unit | yes |
| expansion PV kWp per roof area in m2   | float | yes |         |  | kWp per roof area (how it is in the three dimensional space) | yes |
| expansion PV kWp static                | float | only if `expansion PV sizing mode` is set to `StaticPVSize` |  |  | if `expansion PV sizing mode` is set to `StaticPVSize`, the simulation ignores all roof information and just add a static pv installation per roof of the given size | yes |
| expansion PV static mode profile orientation | string | no |   |  | Profile orientation if static mode is selected | no |
| expansion PV static mode profile index | integer | no |        |  | Profile index if static mode is selected. Especially usefull if an fixed orientation is given | no |

### Group 3, Part B: Battery storage configuration

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| expansion BS P in kW  | float   | yes      |     |    |  Power in kW of the simulatively added batteries. Ignored, if 'expansion BS power computation mode' is set to 'Use E:P-ratio'. | yes |
| expansion BS E in kWh | float   | yes      |  |  | Capacity in kWh of the simulatively added batteries | yes |
| expansion BS initial SOC | float | yes     |  |  | Initial SOC of the simulatively added batteries | yes |
| expansion BS E:P ratio   | float | only if ‘expansion BS power computation mode’ == ‘Use E:P ratio’ |  |  | E:P-ratio for setting the power of simulatively added batteries dependent of the capacity (interesting for parameter variations). This is the inverse of the so-called C-Rate. | yes |
| expansion BS power computation mode | string | yes |  | - `Power as given`<br>- `Use E:P-ratio` | If set to `Use E:P ratio` the variable ‘expansion BS P in kW’ will be ignored, instead $$ P = \frac{E}{expansion BS E:P ratio} $$ | no |
| expansion BS capacity computation mode | string | no | `const` | `const`<br>- `use PV power`<br>- `use mean annual consumption`<br>- `use mean annual consumption with heat pump`<br>- `optimized` | Defines the rule for the calculation of the battery capacity.<br>If set to `const`, the value of `expansion BS E in kWh` will be used.<br>If set to `use PV power`, the capacity of the added storage is the kWp of the PV installation multiplied with the parameter `use PV power` if a simulated PV installation is present. Else it uses the same parameter as `const`.<br>If set to `use mean annual consumption`, it uses the mean annual consumption of the measurement units as defined in the data given.<br>If set to `optimized`, the optimal size is computed base on a optimization. Important: Check option `expansion PV sizing mode` for details. | no |
| expansion BS max capacity | float | no | -1 | 0.0 (excl.) to inf. or -1 | Maxium capacity of an added battery storage (very important for sizing with other modes than constant, but always applied). Defaults to -1.0, meaning there is no limit set. | no |
| expansion BS capacity sizing factor for PV | float | no | 1 | 0.0 (exclusive) to inf. | The sizing parameter of the battery storage if the computation mode is ‘use PV power’ | no |
| expansion BS efficiency in | float | no | 1 | 0.0 to 1.0 | Efficiency of the battery for charging in percent | 
| expansion BS efficiency out | float | no | 1 | 0.0 to 1.0 | Efficiency of the battery for discharging in percent | 
| expansion BS self-discharge per ts | float | no | 0 | 0.0 to 1.0 | Self-discharging rate in percent per time step | 
| expansion BS power for SOC 0 | float | no | 0 | 0.0 to inf | Power consumption of the battery controller if SOC is 0 in kW |  
| expansion BS power for SOC 1 | float | no | 0 | 0.0 to inf | Power consumption of the battery controller if SOC is 1 in kW | 

### Group 3, Part C: Heat pump configuration

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| HP flexibility in ts  | unsigned int | no  | 1                    | 0 to inf        | The number of time steps an existing heat pump profile can be shifted (in both temporal directions) | no |
| th. E to HP el. E conversion factor | float | yes  |              |                 | Factor for converting a thermal energy per building into a heat pump electricity consumption | no |
| Heat consumption apo building V: slope     | float   | no | 0.0   | 0.0-inf         | Parameter for estimating the building heat consumption (in kWh thermal) if no data is given for a building in the data (param: coef. of the linear regression) | no |
| Heat consumption apo building V: intercept | float   | no | 0.0   | 0.0-inf         | Parameter for estimating the building heat consumption (in kWh thermal) if no data is given for a building in the data (param: intercept of the linear regression) | no |

### Group 3, Part D: Individual EV and EV charging station configuration

Notice the difference between individual EVs and the EV charging station (CS): There is only one CS component per control unit, but one single CS can be the home for arbitrary EVs.

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| EV plugin probability | float   | no       | 0.25                 | 0.0 to 1.0      | The probability of plugin in an EV when arraving at home (special cases for low SOC are hard coded) | no |
| EV battery size kWh   | float   | no       | 30.0                 | 0.0 to inf      | The capacity of a simulated EV | no |
| EV consumption kWh per km | float | no     |  0.2                 | 0.0 to inf      | The electricity consumption of an EV for driving 1 km | no |
| EV max charging power | float   | no     |  11.0                  | 0.0 to inf      | The maximum charging power of a simulated EV | no |
| EV charging efficiency | float  | no     |   1.0                  | 0.0 to 1.0      | The efficiency for charging the EV battery   | no |
| CS max power kW       | float   | no       | -1.0                 | 0.0 (ex) to inf | The maximum charging power per charging station component, i.e., per building (regardless of the number of connected EVs). Take care on multi-residential buildings as it holds for them as well. This value will be ignored if the `controller mode` is not set to the rule-based strategy. | no |

### Group 3, Part E: Options applying to all components

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| expansion profile selection | string | yes |             |  as in data / random     | Controls how the profiles (if there are more than one) should be assigned to a new created component | no |

### Group 3, Part F: Other component options

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| open space PV kWp     | float   | yes      |                      | 0.0 to infinity | Additional open space PV power in kWp         | no |
| open space wind kWp   | float   | yes      |                      | 0.0 to infinity | Additional open space wind power in kWp       | no |


## Group 4: Selection of control units for the addition of simulated components

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| break SAC loop if limit reached | bool | no | true                |                 | Should the SAC loop be stopped for an individual combination (like PV + HP) if one of the limits is reached (either PV or HP) (even though HP components should still be added)? | no |
| CU selection mode for comp. add. | string | yes |   | - `as in data`<br>- `random`<br>- `best SSR`<br>- `best NPV`<br>- `use list` | Controls how to select the control units for expansion.<br>`Random`: Shuffle the list of units randomly<br>`As in data`: Use the same sorting present in data (i.e., sorted by unit IDs)<br>In case of `best SSR` and `best NPV` a simulation is executed first, where all SSR and NPV values for every control unit are computed, then they are sorted and finally this sorted list is taken<br>In the case of `use list`, only the units in the given list (see option `cu selection: list of unit IDs to expand`) are used for expansion | no |
| cu selection: list of unit IDs to expand | list of unsigned longs | no | - | - | If the option `CU selection mode for comp. add.` is set to `use list`, this list defines the control unit IDs that can be selected for expansion. Otherwise, this parameter is ignored. | no |
| select buildings with given heat demand only | bool | no | no |                 | Select only buildings where the heat consumption is given in the input database, table ‘heat\_demand\_per\_location’. If set to true, parameters `Heat consumption apo building V: slope` and `Heat consumption apo building V: intercept` will be ignored. | no |
| annual heat demand limit for selection | float | no | -1   | -1 (for no limit) or a value between 1 and inf. | Select only buildings where the annual heat consumption is lower or equal than the given limit (in thermal kWh),  set to -1 (default) if no limit should be choosen | no |
| select only residential buildings | bool | no | no                |                 | If set to yes, only control units representing residential properties will be selected for SAC addition | no |
| expansion PV max total kWp addition    | float | no  | no max. set |  | Maximal size of all added PV installations over all control units | no |
| expansion BS max total E addition | float | no | -1 | -1 (for no limit) or a value >= 0.0 | Upper limit for the accumulated capacity of all added resid. battery storage systems. If this limit has been reached, no further components will be added to other control units, even if the expansion matrix tells so. | no
| expansion BS max total P addition | float | no | -1 | -1 (for no limit) or a value >= 0.0 | Upper limit for the accumulated power of all added resid. battery storage systems. If this limit has been reached, no further components will be added to other control units, even if the expansion matrix tells so. | no
| expansion HP max total addition | unsigned long | no | unset | value >= 0 | Upper limit for the total number of added heat pumps. Once the parameter has been set for a scenario from which it is inherited, it can no longer be set to -1! | no |
| expansion EV max total addition | unsigned long | no | unset | value >= 0 | Upper limit for the total number of added EVs. Once the parameter has been set for a scenario from which it is inherited, it can no longer be set to -1! | no |
| expansion CS max EV per CU      | unsigned long | no | unset | value > 0  | Maximum number of EVs per control unit. Units with more EVs will be excluded. Usefull for optimizations. Once the parameter has been set for a scenario from which it is inherited, it can no longer be set to -1! | no |


## Group 4a: Selection / Deselection of the substations for addition

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| substations without HP addition | list of string | no |           |                 | Substations, where no HP addition can take place (e.g., because of historic buildings) | no |
| substations without PV addition | list of string | no |           |                 | Substations, where no PV addition can take place (e.g., because of historic buildings) | no |
| substations open-space PV potential | dict (substation name: kWp) | no |     |      | Open-space PV potential at a given substation, if an emtpy string (“”) is given as key, this value will be applied for all substations (except special rules are defined AFTER the emtpy string | no |
| substations wind potential | dict (substation name: kW) | no |    |  | TODO DESCRIPTION ! | no |
| substations BS potential   | dict (substation name: [E in kWh, P in kW]) | no |  |  | TODO DESCRIPTION ! | no |


## Group 5: Control strategy settings

| Config Parameter Name | Type    | Required | Default (if not set) | Possible Values | Description | Can be used for parameter variations |
| ------                | ---     | ---      | ---                  | ---             | ------      | ---                                  |
| controller mode       | string  | no       | rule-based | - `rule-based`<br>- `opti with perfect forecast` | Defines the control strategy used inside the control units: rule-based control to optimize (using a linear program) PV self-consumption or an optimization with a perfect forecast | no |
| control horizon in ts | unsigned int | no  | 24                   | 1 to infinity   | If a optimization is used a control strategy, if defines the number of time steps that the optimization looks into the future (i.e., the horizon) | yes |
| control update freq in ts | unsigned int | no | 1                 | 1 to infinity   | If a optimization is used a control strategy, if defines after how many time steps the optimization is executed again inside the control units (the default value of 1 means, that the optimization is executed every time step) | yes |
| controller allow bs charging from grid | string | no | `no grid charging`        | - `no grid charging`<br>- `only grid charging`<br>- `grid charging and discharging` | Controls whether the battery can be charged from the grid and if the battery can also be discharged into the grid. This option is considered only if a `controller mode` is selected other than 'rule-based'. | no |
| controller optimization target | string | no | `electricity costs` | - `electricity costs`<br>- `peak load`<br>- `emissions` | The target of the optimization, if `controller mode` is set to a value other than 'rule-based'.<br>Use `electricity costs` to minimize the costs for electricity consumption (minus revenue for feedin) per control unit.<br>Use `peak load` to minimize the peak load per control unit.<br>Use `emissions` to minimize the CO2 emissions caused by grid demand.<br>If `expansion PV sizing mode` or `expansion BS capacity computation mode` is set to `Optimized`, this parameter must be set to `electricity costs`, and in this case, it also considers (discounted) installation costs for PV and/or BESS respecting the global parameters `net present value time horizon in years` and `net present value discount rate` for discounting. | no |
