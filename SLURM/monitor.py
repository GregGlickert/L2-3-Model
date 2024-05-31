import subprocess

def check_job_status(job_id):
    """
    Checks the status of a SLURM job using scontrol command.

    Args:
        job_id (str): The job ID of the SLURM job.

    Returns:
        str: The status of the job.

    Raises:
        Exception: If there is an error in checking the job status.
    """
    result = subprocess.run(['scontrol', 'show', 'job', job_id], capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f"Error checking job status: {result.stderr}")
    
    # Parse the output to find the job state
    job_info = result.stdout.strip().split('\n')
    job_state = None
    for line in job_info:
        line = line.strip()  # Remove leading and trailing spaces
        if line.startswith('JobState'):
            job_state = line.split('=')[1].strip().split()[0]  # Extract only the state part
            break
    
    if job_state is None:
        raise Exception(f"Failed to retrieve job status for job ID: {job_id}")
    
    return job_state
