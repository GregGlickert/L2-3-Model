from neuron import h
import time
import os

if os.environ.get("OUTPUT_DIR"):
    print(os.environ.get("OUTPUT_DIR"))
time.sleep(15)  # wait 15 seconds
h.nrnmpi_init()       # initialize MPI
pc = h.ParallelContext()
print('I am {} of {}'.format(pc.id(), pc.nhost()))
h.quit()              # necessary to avoid a warning message on parallel exit on some systems