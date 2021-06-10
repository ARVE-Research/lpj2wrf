# lpj2wrf

Simple Fortran program and helper scripts and namelists for converting [LPJ-LMfire](https://github.com/ARVE-Research/LPJ-LMfire) dynamic global vegetation model output to land surface boundary conditions for the Weather Research and Forecasting [WRF](https://github.com/wrf-model/WRF) regional climate model.

The program generates the following variables as output (for input to WRF):

|no.|variable|
|---|---|
|1.|land cover fraction (category)|
|2.|dominant land cover type (category)|
|3.|FPAR (monthly)|
|4.|LAI (monthly)|
|5.|albedo (monthly)|
|6.|soil temperature (annual mean)|
|7.|snow albedo (annual max)|
