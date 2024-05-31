import os

class SlurmJob:
    """
    Class to represent a SLURM job.

    Attributes:
        job_name (str): The name of the job.
        output_file (str): The name of the output file for job logs.
        time (str): The maximum runtime for the job.
        partition (str): The partition to run the job on.
        nodes (int): The number of nodes required.
        ntasks (int): The number of tasks (or cores) required.
        commands (list): List of commands to be run in the job.
    """

    def __init__(self, job_name, output_file, time, partition, nodes, ntasks, commands):
        self.job_name = job_name
        self.output_file = output_file
        self.time = time
        self.partition = partition
        self.nodes = nodes
        self.ntasks = ntasks
        self.commands = commands

    def generate_script(self):
        """
        Generates the SLURM job script.

        Returns:
            str: The SLURM job script as a string.
        """
        command_str = "\n".join(self.commands)
        script = f"""#!/bin/bash
#SBATCH --job-name={self.job_name}
#SBATCH --output={self.output_file}
#SBATCH --time={self.time}
#SBATCH --partition={self.partition}
#SBATCH --nodes={self.nodes}
#SBATCH --ntasks={self.ntasks}

{command_str}
"""
        return script

    def save_script(self, directory='scripts'):
        """
        Saves the SLURM job script to a file.

        Args:
            directory (str): The directory to save the script in. Defaults to 'scripts'.

        Returns:
            str: The path to the saved script.
        """
        if not os.path.exists(directory):
            os.makedirs(directory)
        script_path = os.path.join(directory, f"{self.job_name}.sh")
        with open(script_path, 'w') as file:
            file.write(self.generate_script())
        return script_path
