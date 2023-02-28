# URANOS Setup

How to setup and run URANOS jobs with examples and helpful scripts.

![grafik](https://user-images.githubusercontent.com/7942719/221894176-f12bdb56-5929-498f-a0ae-088d0ffd7a5f.png)

## 1. Graphical User Interface (GUI)

1. Double-click on `UranosGUI.exe`,
2. Set the folders and parameters,
3. Save config
4. Start the simulation
5. Repeat from 1 with another scenario with different parameters.

This can be cumbersome for large sets of scenarios.

## 2. Command line (noGUI)

Uranos can be run without GUI:

    UranosGUI.exe noGUI path/to/Uranos.cfg

There are several ways to create scenarios and jobs, here are a few examples:

## 3. Prepare CFG

To prepare several scenario folders with its Uranos.cfg files, adapt `prepare_scenarios.py` to your needs and execute it. E.g., it can create a bunch of scenarios with different soil moisture values.

If you have an existing set of folders, and just want to change path names or Uranos version paths, you can adapt and run `change_paths.py`.

## 4a. Run URANOS noGUI on Windows

Assuming you have folder with subfolders for each scenario. You can now:
- create batch files with a list of scenarios using `generate_batch_files.py`, or
- create a list of files including a list of scenarios and run everything with python.
  1. Adapt and run `generate_scenario_list.py -s 100` to create 100 files including equal number of scenarios, e.g.,
     ```
     uranos_scenario_list-01.txt
     uranos_scenario_list-02.txt
     uranos_scenario_list-03.txt
     ...
     ```
     The number of files should be chosen based on the number of available COU cores on your system.
  2. Adapt and run `python run_scenarios.py` in a separate command line window for every scenario list, e.g.
     ```
     python run_scenarios.py -f uranos_scenario_list-01.txt
     ```
  3. Each python script will forward the last line of the Uranos output to a file in the folder `monitor/`. Double-click `monitor.py` to show a command line window that updates reguarly to show the latest state of the files. E.g.:
     ```
     12% uslist-project-01.txt:    3/16 project/sm05   6.1% (-14.4h)
     98% uslist-project-02.txt:   15/16 project/sm15   1.1% (-12.4h)
     ...
     23% uslist-project-80.txt:    3/15 project/sm30  50.9% (- 6.2h)
     |>>>>>>>>>>                                                     |  10%
     Monitoring 80 processes:  127/1232 scenarios finished, average duration: 15.4h, total ETA:  8.9 days.
     ```
  The script can also be configured to automatically upload the current status to an FTP server.
  4. To monitor status of a remote job whlie staying on a local computer, edit `check_uranos.py` and set FTP settings, then run `python check_uranos.py`. Repeat on demand.
  
## 4b. Run URANOS noGUI on a Linux Cluster

Assuming you have folder with subfolders for each scenario. You can now:

1. Adapt and run `submit_scenarios.sub` to submit the list of 99 scenarios given by `scenarios.files` to a SLURM job system:
  ```
  sbatch -a 1-99 submit_scenarios.sub
  ```
2. The running jobs will create `*.out` files in a log folder with current URANOS output messages. Execute `./check_uranos.sh` to see the current progress of all jobs.
  ```
  95.988 % completed  (0:34:49) /work/URANOS/log/123-01.out
     100 % : deleting histos... /work/URANOS/log/123-02.out
  ...
  72.894 % completed  (5:09:57) /work/URANOS/log/123-99.out
  ```

## 5. Quick-collect results into a CSV

Adapt and run `collect-results.py` to quickly walk through all finished scenarios and extract infos like average neutron counts, soil moisture, and penetration depth. They will be summarized in a CSV file:
```
id, sm,     N, N_std, N_sem,  D86
 0,  1, 26.45,  5.57, 0.005, 1.09
 1,  5, 18.94,  4.52, 0.009, 1.03
...
99, 30, 10.62,  3.74, 0.010, 0.85
```
