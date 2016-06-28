'''
misopy.cluster.slurm

@author: Aaron Kitzmiller
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
@contact: aaron_kitzmiller@harvard.edu
'''
import os, subprocess, traceback

from settings import load_settings



class SlurmClusterEngine():
    '''
    Run jobs on a Slurm cluster
    '''
    
    def __init__(self,settings_filename):
        '''
        Get the slurm job template file and load it
        '''
        self.settings = load_settings(settings_filename)
        if not 'slurm_template' in self.settings:
            raise Exception('slurm_template must be defined in settings to use Slurm')
        
        self.squeue_max_attempts = 10
        if 'squeue_max_attempts' in self.settings:
            self.squeue_max_attempts = int(self.settings['squeue_max_attempts'])
  
        template_filename = self.settings['slurm_template']
        if not os.path.exists(template_filename):
            raise Exception('Cannot find slurm template file %s' % template_filename)
        
        with open(template_filename,'r') as f:
            self.template = f.read()
            
        if self.template.strip() == '':
            raise Exception('slurm template file %s is empty' % template_filename)
        
            
    
    def run_on_cluster(self, cmd, job_name, cluster_output_dir,
                       cluster_scripts_dir=None,
                       queue_type=None,
                       settings_fname=None):
        pass
    
    def wait_on_job(self, job_id, cluster_cmd, delay=60):
        '''
        Wait until job is done.  Uses squeue first, then sacct.
        Runs squeue /sacct until either the job is done or until squeue_max_attempts is reached.
        Max attempts is needed to ensure that squeue information is available.
        '''
        squeue_cmd = 'squeue --noheader --format %%T -j %d' % job_id
        sacct_cmd  = 'sacct --noheader --format State -j %d.batch' % job_id

        done = False
        squeue_attempts = 0
        state = None

        while not done:

            # If we've tried squeue_max_attempts and gotten no information, then quit
            squeue_attempts += 1
            if squeue_attempts == self.squeue_max_attempts and state is None:
                raise Exception('Attempted to query squeue /sacct %d times and retrieved no result' % squeue_attempts)
            
            proc = subprocess.Popen(squeue_cmd, shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    stdin=subprocess.PIPE)
            output,err = proc.communicate()
            
            if proc.returncode != 0:
                # The whole command failed.  Weird.
                raise Exception('squeue command %s failed: %s' % (squeue_cmd,err))
            
            if output.strip() != '':
                state = output.strip()
            else:
                # Try sacct.  The job may be done and so disappeared from squeue
                proc = subprocess.Popen(sacct_cmd, shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        stdin=subprocess.PIPE)
                output,err = proc.communicate()
                
                if proc.returncode != 0:
                    # The whole command failed.  Weird.
                    raise Exception('sacct command %s failed: %s' % (sacct_cmd,err))
                
                if output.strip() != '':
                    state = output.strip()
                    
            if state is not None:
                if state in ["COMPLETED","COMPLETING","CANCELLED","FAILED","TIMEOUT","PREEMPTED","NODE_FAIL"]:
                    done = True
                
                
                
        
        
    
    def launch_job(self, cluster_cmd):
        '''
        Runs cluster command and returns the job id
        '''
        proc = subprocess.Popen(cluster_cmd, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE)
        # Read the job ID if it's a known cluster
        # submission system
        output,err = proc.communicate()
        if proc.returncode != 0 or err.strip() != '':
            raise Exception('Error launching job with %s: %s' % (cluster_cmd,err))
        
        job_id = output.strip().replace('Submitted batch job ','')
        try:
            job_id = int(job_id)
        except Exception as e:
            raise Exception('Returned job id %s is not a number ?!?!' % job_id)
        
        return job_id