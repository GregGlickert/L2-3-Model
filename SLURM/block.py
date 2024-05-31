from job import SlurmJob
from submit import submit_job
from monitor import check_job_status
import time

class SimulationBlock:
    """
    Class to represent a block of SLURM jobs.

    Attributes:
        block_name (str): The name of the block.
        base_output_file (str): The base name for the output files.
        time (str): The maximum runtime for the jobs.
        partition (str): The partition to run the jobs on.
        nodes (int): The number of nodes required.
        ntasks (int): The number of tasks (or cores) required.
        simulation_cases (dict): Dictionary of simulation cases and their corresponding commands.
    """

    def __init__(self, block_name, base_output_file, time, partition, nodes, ntasks, simulation_cases):
        self.block_name = block_name
        self.base_output_file = base_output_file
        self.time = time
        self.partition = partition
        self.nodes = nodes
        self.ntasks = ntasks
        self.simulation_cases = simulation_cases
        self.job_ids = []

    def create_jobs(self):
        """
        Creates SLURM jobs for each simulation case.

        Returns:
            list: A list of SlurmJob instances.
        """
        jobs = []
        for case_name, command in self.simulation_cases.items():
            job_name = f"{self.block_name}_{case_name}"
            output_file = f"{self.base_output_file}_{case_name}.txt"
            job = SlurmJob(job_name, output_file, self.time, self.partition, self.nodes, self.ntasks, [command])
            jobs.append(job)
        return jobs

    def submit_block(self):
        """
        Submits all SLURM jobs in the block.

        Returns:
            list: A list of job IDs for the submitted jobs.
        """
        jobs = self.create_jobs()
        self.job_ids = []
        for job in jobs:
            script_path = job.save_script()
            job_id = submit_job(script_path)
            self.job_ids.append(job_id)
        return self.job_ids

    def check_block_status(self):
        """
        Checks the status of all jobs in the block.

        Returns:
            bool: True if all jobs in the block are completed, False otherwise.
        """

        for job_id in self.job_ids:
            status = check_job_status(job_id)
            if status not in ['COMPLETED', 'FAILED', 'CANCELLED']: # could throw an error if failed or cancelled
                return False
        return True


class SequentialBlockRunner:
    """
    Class to handle submitting multiple blocks sequentially.

    Attributes:
        blocks (list): List of SimulationBlock instances to be run.
        checkDone (int): interval in seconds to check if jobs are complete or not
    """

    def __init__(self, blocks,checkDone=60):
        self.blocks = blocks
        self.checkDone = checkDone
    def submit_blocks_sequentially(self):
        """
        Submits all blocks sequentially, ensuring each block starts only after the previous block has completed.
        """
        for block in self.blocks:
            print(f"Submitting block: {block.block_name}")
            block.submit_block()
            while not block.check_block_status():
                print(f"Waiting for block {block.block_name} to complete...")
                time.sleep(self.checkDone)  # don't want this to constantly run we we will wait a bit 
            print(f"Block {block.block_name} completed.")
