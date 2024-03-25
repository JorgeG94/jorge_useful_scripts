#!/bin/bash


# Loop through all .pbs files and submit them using qsub
for script in *.pbs; do
    qsub "$script"
done

