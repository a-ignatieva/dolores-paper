#!/bin/bash

for invl in 50000 100000 200000
do
	for samp in 50
	do
		for i in {0..99}
		do
			cd /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations
			mkdir slim_${invl}_${samp}_${i}
			cd slim_${invl}_${samp}_${i}
			pwd
			cp /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/inversion_testing.py .
			cp /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/run-dolores.py .
			cp /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/run-invclust.R .
			cp /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/slim-sim-with-inversion.slim .
			cp /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations/slim-sim-neutral.slim .
			python -m inversion_testing ${invl} ${samp} ${i}
			python -m run-dolores -C chr0 -n slim_${invl}_${samp}_${i}_relate -t . -s HomSap -g dummy_map.txt -c 0.01 -m 0.05 -u 0
			Rscript run-invclust.R ${invl} slim_${invl}_${samp}_${i}
			cd /Users/ignatiev/Dropbox/projects/branch-durations/slim-simulations
			pwd
		done
	done
done