# %%
"""
Collects results from the finished URANOS scenarios
and summarieses them in a CSV file.
"""
# %%
# 1. Specify work folder which includes scenario folders
work_folder = '/work/username/URANOS/my_project/'

# %%
import numpy as np
import pandas
from glob import glob
from scipy.stats import sem
from lib.uranos import URANOS

# %%
# 2. Specify list of scenarios here (same as in prepare_scenarios.py)

print('Crawling scenarios...', flush=True)
soil_moistures = np.array([1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80])

list_of_rows = []
for sm in soil_moistures:
    scenario_name = 'sm%02d' % sm
    mypath = work_folder + '%s/' % scenario_name
    print(mypath, flush=True)

    # If scenario folder contains a CSV file (i.e., simulation is finished)
    if glob(mypath + '*.csv'):
        # Read URANOS data
        U = URANOS(folder=mypath, scaling=2, hum=5)
        U = U.read_density('densityMapSelected*', pad=True)
        #U = U.read_origins('detectorOrigins*', pad=True)
        U = U.read_root('uranosRawHistos*', show_vars=False)
        U = U.read_root_var('detectorDistanceDepth2;1')
        D86 = U.D86_from_root()
        #U.plot(image='Density')

        # Add more analysis here, e.g. footprint

        newdf = pandas.DataFrame([
            [sm, np.mean(U.Density), np.std(U.Density), sem(U.Density, axis=None), D86]
        ], columns=['sm','N','N_std','N_sem','D86'])

        list_of_rows.append(newdf)

# %%
# 3. Combine tables and export to CSV
if len(list_of_rows)>0:
    df = pandas.concat(list_of_rows, ignore_index=True)
    print("\n")

    df.to_csv(work_folder + 'results.csv')
