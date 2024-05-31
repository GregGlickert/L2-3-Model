from neuron import h
import time
time.sleep(30)  # wait 30 seconds
h.nrnmpi_init()       # initialize MPI
pc = h.ParallelContext()
print('I am {} of {}'.format(pc.id(), pc.nhost()))
h.quit()              # necessary to avoid a warning message on parallel exit on some systems