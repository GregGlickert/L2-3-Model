# To do list for batch module
## Goal of this module is to simplify the process for running several jobs with different parameters easily without having to activity monitor the simulation. Goal use would be to at night queue up 10s or 100s of simulations with different parameters and then in the morning have them all organized for analysis.


* ## How to currently use, First set up blocks in main.py. Then submit start_blocks.sh on SLURM. You can also just do python main.py. However when actually using this module blocks may take a while so it can be good to run main.py on the background on SLURM. 

* ## Need to add feature where two json files (think CS and CP) can be changed at the same time

* ## Need to be able to change a new json file. This way you can queue up loads of runs changing different jsons indivually without having to run main.py more than once. 


* ## While this is super not needed could add simAgent support https://github.com/tjbanks/SimAgentMPI/tree/master and make a super cool application for running BMTK jobs on HPC. This could be a LOT of worth so don't know if worth it but SimAgent is cool.