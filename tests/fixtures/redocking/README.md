# Redocking regression data

CCD definitions (`92H`, `BUA`, `S6Q`) are RCSB Chemical Component Dictionary
records. BUA exercises hydroxyl vs carbonyl typing; S6Q exercises substituted
pyrrole nitrogen typing before native Chimera AddH.

`1XLU_ISP_1001A.hydrogenated.mol2` and
`3VF3_0GS_401A.hydrogenated.mol2` are small, pre-charge Chimera exports from the
configured Alzheimer validation pairs (2026-10-27), with crystallographic heavy
coordinates and native hydrogen placement. They reproduce cases where Chimera's
charge heuristic gave 0 instead of +1 (observed ISP fragment) and +2 instead of 0
(0GS), respectively. They are topology/charge regression fixtures, not docking
acceptance references or claims of a preferred physiological protonation state.

The older 5NN0 fixtures exercise export validation and rejection of lost hydrogen
inventory. Tests run offline using these files and never modify the saved project.
