# %%
"""
Creates a number of scenario folders incl. Uranos.cfg and GeometryConfig.dat,
based on a loop over certain variables, such as a range of soil moisture.
Also creates a file containing a list of these paths.
"""
import os
import numpy as np

# %%
# 1. Set your paths here
work_folder   = '/work/username/URANOS/my_project'
output_folder = '/work/username/URANOS/my_project'
spectrum_file = '/software/easybuild-broadwell/software/URANOS/1.07-foss-2020b/lib/IncomingSpectrum.root'
endf_folder   = '/software/easybuild-broadwell/software/URANOS/1.07-foss-2020b/lib/ENDF/'

# %%
# 2. Set your Uranos configuration here.
uranoscfg = """%s/%s/
%s/%s/
%s
default
%s
40000000	 Number of Neutrons
1000000.000	 Dimension [mm]
500000.000	 Beam Radius [mm]
531	 #Neutrons Refresh Rate
0.000005000	 Lower THL
1.000000000	 Higher THL
90000.000	 Detector Radius [mm]
0.000	 Detector Pos X [mm]
0.000	 Detector Pos Y [mm]
0	 Detector Absorbing [bool]
0	 Detector File Output [bool]
%.3f	 Soil Volumetric Water Fraction
5.0	 Air Humidity [g/m3]
0.50	 Soil Porosity
1013.000	 Atmospheric Depth [g/cm2]
5.000	 Cutoff Rigidity [GeV]
0	 #Neutrons Precalulated Spectrum [Power]
0	 Radial Source [bool]
0	 Volume Source [bool]
1	 HE Cascade Model [bool]
270.000	 Downward Scotoma Angle [0-360]
360.000	 Downward Acceptance Angle [0-360]
0	 Only Record in Material [bool]
0	 Record Material No [bool]
0	 Exclude Multiple Scattering Recording [bool]
0	 Track all Layers [bool]
1	 Layer Physics Model [bool]
1	 Detector Physics Model [bool]
0	 Cylindrical Detector [bool]
1	 Spherical Detector [bool]
1	 ROOT Output [bool]
0	 Separate Folder Each Export [bool]
0	 Export Epithermal Map [bool]
0	 Export Epithermal Data [bool]
0	 Export Intermediat Map [bool]
0	 Export Intermediat Data [bool]
0	 Export Fast Map [bool]
0	 Export Fast Data [bool]
0	 Export Detector Energy Map [bool]
1	 Export Detector Energy Data [bool]
0	 Export Detector Origins Map [bool]
0	 Export Detector Origins Data [bool]
0	 Export Detector Distance Data [bool]
0	 Export Detector Layer Distance Data [bool]
0	 Export x Map [bool]
0	 Export x Data [bool]
0	 Export Travel Distance Graph [bool]
0	 Export Detector Layer File Output [bool]
1	 Export Detector Physics Ext. Mod. [bool]
0	 Use Layer Maps [bool]
0	 Sideways Tracking [bool]
0.000	 Sideways Tracking y +/-Cutout [mm]
0	 Clear Every x Neutrons [bool]
1	 Clear Every x Neutrons Number
1	 Refresh Rate Auto Update [bool]
1.000	  Refresh Rate Auto Update Time [s]
260909.000
-154688.000
0
0
0
1.000
0.000
0	 Set Clearing to Refresh Rate [bool]
0	 Export Neutron Track Data [bool]
0	 Export High Resolution Neutron Track Data [bool]
0	 Detector Sheet along y-Axis [bool]
0	 Detector Sheet along y-Axis [bool]
0.000	 Detector Sheet Length [mm]
0	 use Domain Cutoff Factor [bool]
2.000	 Domain Cutoff Factor [float]
0	 use Domain Cutoff Distance [bool]
0	 Domain Cutoff Distance [m]
0	 Detector Neutron Track File Output [bool]
0	 All Neutron Track File Output [bool]
13	 Energy Display Range for Birds-Eye View [int]
3.000	0.000	1.000	 Plant Gas Density Multiplicator [kg/m3]	 Plant dry density [g/cm3] 	 Plant water density [g/cm3] [float]
0	 Reflective Boundary Conditions [bool]
1	 Periodic Boundary Conditions [bool]
0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000 Element Density in Soil #19 [1e-6 g/cm3]

""" 

geometrycfg = """2
4
6
-1000	920	11
-80	30	11
-50	47.5	11
-2.5	0.5	11
-2	2	11
0	3.0	20

"""

# %%
# 3. Generate Scenatios

# Create a number of scenarios with different soil moisture values
soil_moistures = np.array([1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80])

scenario_names = []
for sm in soil_moistures:
    
    scenario_name = 'sm%02d' % sm
    scenario_names.append(scenario_name)
    scenario_path = '%s/%s/' % (work_folder, scenario_name)

    # Create folders
    os.makedirs(scenario_path, exist_ok=True)
    # Write Uranos.cfg
    with open(scenario_path + 'Uranos.cfg', 'w') as f:
        f.write(uranoscfg % (work_folder, scenario_name, work_folder, scenario_name, spectrum_file, endf_folder, sm/100))
    # Write UranosGeometryConfig.dat
    with open(scenario_path + 'UranosGeometryConfig.dat', 'w') as f:
        f.write(geometrycfg)

# %%
# 4. Write file with list of paths to each scenario

with open('%s/scenario.files' % work_folder, 'w', newline='\n') as f:
    for scenario_name in scenario_names:
        f.write('%s/%s/Uranos.cfg\n' % (work_folder, scenario_name))
