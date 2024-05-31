import sys
import os
import warnings
import synapses
from bmtk.simulator import bionet
from bmtk.simulator.bionet.pyfunction_cache import add_weight_function
from neuron import h

CONFIG = 'config.json'
USE_CORENEURON = False

import os
import json

def get_synaptic_params(path_to_syn_folder):
    """Gets values of all json files and puts them into one file"""
    combined_data = []
    # List all files in the input folder
    for filename in os.listdir(path_to_syn_folder):
        if filename.endswith('.json'):
            file_path = os.path.join(path_to_syn_folder, filename)
            
            # Open and read the JSON file
            with open(file_path, 'r') as file:
                data = json.load(file)
                # Append the filename and its data
                combined_data.append({
                    'filename': filename,
                    'data': data
                })
    return combined_data

def save_synaptic_params(data,path_to_output_dir):
    """Saves combined json data into one file"""
    with open(path_to_output_dir, 'w') as output_file:
        json.dump(data, output_file, indent=4)

def run(config_file=CONFIG, use_coreneuron=USE_CORENEURON):

    warnings.simplefilter(action='ignore', category=FutureWarning)

    # want to put report in output file so hope this wont get overwritten
    with open(config_file, 'r') as json_file:
        conf_dict = json.load(json_file)
        output_dir = conf_dict['manifest']['$OUTPUT_DIR']
        output_dir = output_dir + "/synaptic_report.json"
        syn_data = get_synaptic_params('components/synaptic_models/synapses_STP')

    # register synaptic weight function
    synapses.load(randseed=1111)
    add_weight_function(synapses.lognormal_weight, name='lognormal_weight')

    if use_coreneuron:
        import corebmtk
        conf = corebmtk.Config.from_json(config_file, validate=True)
    else:
        conf = bionet.Config.from_json(config_file, validate=True)

    conf.build_env()
    graph = bionet.BioNetwork.from_config(conf)

    if use_coreneuron:
        sim = corebmtk.CoreBioSimulator.from_config(
            conf, network=graph, gpu=False)
    else:
        sim = bionet.BioSimulator.from_config(conf, network=graph)

    '''
    # This calls insert_mechs() on each cell to use its gid as a seed
    # to the random number generator, so that each cell gets a different
    # random seed for the point-conductance noise
    cells = graph.get_local_cells()
    for cell in cells:
        cells[cell].hobj.insert_mechs(cells[cell].gid)
    '''

    # clear ecp temporary directory to avoid errors
    pc = h.ParallelContext()
    if pc.id() == 0:
        try:
            ecp_tmp = conf['reports']['ecp']['tmp_dir']
        except:
            pass
        else:
            if os.path.isdir(ecp_tmp):
                for f in os.listdir(ecp_tmp):
                    if f.endswith(".h5"):
                        try:
                            os.remove(os.path.join(ecp_tmp, f))
                        except Exception as e:
                            print(f'Failed to delete {f}. {e}')
    pc.barrier()


    sim.run()
    # must be ran after sim.run since that creates dir
    if pc.id() == 0:
        print(output_dir)
        print(syn_data)
        save_synaptic_params(syn_data,output_dir)

    bionet.nrn.quit_execution()


if __name__ == '__main__':
    for i, s in enumerate(sys.argv):
        if s in __file__:
            break

    if i < len(sys.argv) - 1:
        argv = sys.argv[i + 1:]
        for i in range(1, len(argv)):
            argv[i] = eval(argv[i])
        run(*argv)
    else:
        run()
