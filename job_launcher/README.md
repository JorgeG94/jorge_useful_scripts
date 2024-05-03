Given an xyz file, generates inputs and job scripts for qchem, nwchem, orca. 
Also launches jobs and collects results.

Example usage: 48 cores on the gadi normal queue
```
python3 experiment_launcher.py --name w40 -i w40.xyz -j pbs -q normal -c 48
```
