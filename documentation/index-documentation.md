# Documentation of UrbanFireXDT

Welcome to the documentation of the **Urban Future Integrated Renewable Energy systems eXploration Digital Twin (UrbanFireXDT)**.

This documentation is divided into two major parts:

---

## Operational Documentation

For *users* of the simulation framework:

- [configuration-options.md](configuration-options.md):
  Explains all available configuration options that can be defined inside the json config file.
- [input-description.md](input-description.md):
  Describes the structure of the most important input file: the structure of the modeled energy system.
- [output-description.md](output-description.md):
  Describes all generated output files and the meaning of all columns in these files.

To finally start the framework, a user requires the following parts:
- A database containing the system structure of the simulated energy system
- The separated (electricity) smart meter files
- A configuration file describing the future development of the energy system (i.e., the simulation target)
- A scenario file referenced in the configuration file describing the targets of the simulateively added components to control units.

A minimal example can be found in the [test scenarios](https://github.com/danielrenebayer/UrbanFireXDT/tree/main/test) in the corresponding repository.

---

## Developer Documentation

For *developers* extending or modifying the simulation framework:

- **Architecture Overview:**
  The code follows a modular component design (\ref Substation "Substations", \ref ControlUnit "Control Units", \ref MeasurementUnit "Measurement Units", etc.).
- **Code Reference:**
  All namespaces, classes, and methods are documented under the *Classes* and *Files* tabs of this documentation.

An overview of all major classes (but not all classes) of the code is given in the image below ([click here](./ClassDiagram.svg) for a SVG-version).

![Class diagram](./ClassDiagram.png)

---

## Repository Structure

For informations on the repository structure, see the [README.md](https://github.com/danielrenebayer/UrbanFireXDT) file on the top level.


