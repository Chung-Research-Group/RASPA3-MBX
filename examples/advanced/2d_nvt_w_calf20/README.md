# H₂O–CO₂ NVT-W window in CALF-20

Run `raspa3` from this directory. The conventional force field describes four H₂O and four CO₂ molecules in a 3 × 3 × 3 CALF-20 supercell at 295 K and a reservoir pressure of 97,000 Pa. The input runs 2 initialization cycles, 2 equilibration cycles and 20 production cycles as a short functionality test.

Each component sets its minimum and maximum macrostate to 4. Translation, rotation and reinsertion sample configurations within the fixed population. The sampler evaluates insertion and deletion proposals for both species and retains their unbiased virtual transition weights at the window boundaries.

Read `tmmc/tmmc_statistics.txt`. The single production row gives the state index, H₂O and CO₂ populations, the stay weight, H₂O deletion and insertion weights, CO₂ deletion and insertion weights, bias, histogram estimate of ln π, retained visit count and mean potential energy in native units. The five collection weights sum to the number of production attempts, up to floating-point rounding.

The CIF, force field and component definitions come from `tmmc/nvt_w/base` in Chung-Research-Group/RASPA3-MBX at commit `0140ee08d18f90e42eb12a409b3aeb6410c68147`. The simulation input selects its (4, 4) window and shortens the run for this example.
