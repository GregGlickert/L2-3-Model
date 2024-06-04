from block import SimulationBlock, SequentialBlockRunner
from seedSweep import seedSweep

# Define simulation cases
# simulation_cases (dict): Dictionary of simulation cases and their corresponding commands.
simulation_cases = {
    "baseline": "mpiexec nrniv -mpi -python batch_test1.py",
    "short": "mpiexec nrniv -mpi -python batch_test2.py",
    "long": "mpiexec nrniv -mpi -python batch_test3.py"
}

# Define block parameters
# could add other params like account to allow working on Expanse
block_params = {
    'time': '01:00:00',
    'partition': 'batch',
    'nodes': 1,
    'ntasks': 2,
    'mem': 1,
    'output_base_dir': 'output_test/'
}

# Define the number of blocks to create
num_blocks = 3

# Create a list to hold the blocks
blocks = []

# Create blocks with the defined simulation cases
for i in range(1, num_blocks + 1):
    block_name = f'block{i}'
    block = SimulationBlock(block_name, **block_params, simulation_cases=simulation_cases)
    blocks.append(block)

# Define the parameter to be changed and its values
param_name = 'initW'
param_values = [100, 200, 300]  # Example values for the parameter could also be a loop

# Define JSON file path and create seedSweep instance
json_file_path = 'test.json'
json_editor = seedSweep(json_file_path, param_name)


# Create a SequentialBlockRunner with the blocks and parameter changes
#runner = SequentialBlockRunner(blocks, json_editor, param_values,check_interval=10)
runner = SequentialBlockRunner(blocks,check_interval=10)

# Submit the blocks sequentially
runner.submit_blocks_sequentially()