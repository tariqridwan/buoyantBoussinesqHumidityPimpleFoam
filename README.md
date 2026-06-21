# buoyantBoussinesqHumidityPimpleFoam ✎
A transient `OpenFOAM` solver with specific humidity transport for buoyant 🍃, turbulent flow 🌪 of incompressible fluids based on the Boussinesq approximation.<br />
Accompanied by a moisture flux boundary 🚧 condition: `atmTurbulentMoistureFluxHumidity`.

## Documentation 📜
Documentation [➚](https://github.com/tariqridwan/buoyantBoussinesqHumidityPimpleFoam/blob/main/doc/Solver_BC_equations.pdf) on the governing equations, boundary condition, and their implementation in OpenFOAM.<br />
Follow the provided test cases to see what to add in `0/sh`, `constant/transportProperties` and in `system` files.

## Compatibility 🧩
The solver and boundary condition are tested 🛠️ and verified on the ESI (.com) versions: `v2312`, `v2406`, `v2506` and `v2512`. Not tested in the Foundation (.org) `vXX` versions.<br />
Follow `Installation_instruction.txt` for installation.

## Funding 💶
This work was funded by the **Ministerio de Ciencia, Innovación y Universidades (MICIU)** of Spain 🇪🇸 through the **Formación de Personal Investigador (FPI)** fellowship 🎓 number [`PRE2020-095284`](https://www.aei.gob.es/en/announcements/announcements-finder/ayudas-contratos-predoctorales-formacion-doctores-2020) for the national project: *"Turbulence and large coherent structures in the atmospheric boundary layer: Fundamental aspects for parametrizations of cloud formation and for wind-energy applications"* [(TABL4CW)](https://futur.upc.edu/28980349).

## Acknowledgements 🤝🏼
The developers gratefully acknowledge the computational resources and technical support provided by the **Red Española de Supercomputación (RES)** and **Barcelona Supercomputing Center (BSC-CNS)**. CFD simulations associated with the development and testing of this solver were performed on the **MareNostrum 5** 🖥️ supercomputer under allocation [`IM-2025-3-0017`](https://www.bsc.es/res-intranet/files/resolution/resolucioncomitedeacceso3erperiodo2025.pdf).
