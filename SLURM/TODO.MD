# To do list for batch module
## Goal of this module is to simplify the process for running several jobs with different parameters easily without having to activity monitor the simulation. Goal use would be to at night queue up 10s or 100s of simulations with different parameters and then in the morning have them all organized for analysis.

* ## Could change the monitor to check if the job is running instead of finished. This would allow a lot more jobs to be running at once and could be very useful on Expanse or other large HPCs.

* ## add feature to change synaptic parameters. This will need to be done right before a block is ran and then modify again before the second block is RAN. Not when the blocks are added to the block queue. Only when the simulation is started for one block should the next blocks synapses be changes. Look at seedSweep.py for some ideas

* ## Network-Parameters/run_network.py has two new functions that save all synaptic properties of all synapses into the output dir of the simulation. depending on how the seed sweep is added may want to move those functions into the sweep class or something

* ## make sure output path for simulation results and reports is set to correct location

* ## While this is super not needed could add simAgent support https://github.com/tjbanks/SimAgentMPI/tree/master and make a super cool application for running BMTK jobs on HPC. 