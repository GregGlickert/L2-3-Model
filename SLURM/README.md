# How to use SLURM module

## Blocks
### A block is a set of simulations that will be ran with the same parameters, but a different simulation config. Currently the parameters supported to change are synaptic parameters defined in json files. First define the simulation cases which is the name of the case and then then command needed to run the simulation.
```python 
simulation_cases = {
    "baseline": "mpirun ./components/mechanisms/x86_64/special -mpi -python run_network.py simulation_config_baseline.json True",
    "short": "mpirun ./components/mechanisms/x86_64/special -mpi -python run_network.py simulation_config_short.json True",
    "long": "mpirun ./components/mechanisms/x86_64/special -mpi -python run_network.py simulation_config_long.json True"
}
```
### Then you need to define the block parameters. These are the parameters that SLURM will use to run the simulation
```python 
block_params = {
    'time': '02:00:00',
    'partition': 'batch',
    'nodes': 1,
    'ntasks': 24,
    'mem': 48,
    'output_base_dir': '../Run-Storage/Test',
    'account':'umc113'
}
```
### Then you can set up the blocks. There are more parameters in the SimulationBlock class so check that out.
```python 
from block import SimulationBlock
blocks = []
for i in range(1, num_blocks + 1):
    block_name = f'block{i}'
    block = SimulationBlock(block_name, **block_params, simulation_cases=simulation_cases)
    blocks.append(block)
```
### Next you will want to set up the Json editor. This will allow for synaptic parameters to be adjusted before each block is ran.
```python 
from seedSweep import seedSweep
param_name = 'initW'
param_values = [0.85,0.95,1.05,1.15,1.25] 
json_file_path = 'Json_path'
json_editor = seedSweep(json_file_path, param_name)
```
### Once the editor is set up then you can use the block runner. This will edit the synaptic parameters, and submit all cases in the block. It will then wait for all cases to finish(or start running depending on settings) before submitting the next block.
```python 
runner = SequentialBlockRunner(blocks, json_editor, param_values,check_interval=200)
runner.submit_blocks_sequentially()
```

