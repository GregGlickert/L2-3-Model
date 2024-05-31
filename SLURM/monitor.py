import subprocess

def check_job_status(job_id):
    """
    Checks the status of a SLURM job.

    Args:
        job_id (str): The SLURM job ID.

    Returns:
        str: The job status.
    """
    result = subprocess.run(['scontrol', 'show', 'job', job_id], capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f"Error checking job status: {result.stderr}")

    job_info = result.stdout
    for line in job_info.split('\n'):
        if 'JobState=' in line:
            job_state = line.split('JobState=')[1].split()[0]
            return job_state
    raise Exception(f"Failed to retrieve job status for job ID: {job_id}")
