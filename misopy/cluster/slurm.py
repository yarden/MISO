'''
misopy.cluster.slurm

@author: Aaron Kitzmiller
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
@contact: aaron_kitzmiller@harvard.edu
'''
import os

from settings import load_settings



class SlurmClusterEngine():
    '''
    Run jobs on a Slurm cluster
    '''
    
    def __init__(self,settings_filename):
        '''
        Get the slurm job template file and load it
        '''
        settings = load_settings(settings_filename)
        if not 'slurm_template' in settings:
            raise Exception('slurm_template must be defined in settings to use Slurm')
        
        template_filename = settings['slurm_template']
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
        pass
    
    def wait_on_jobs(self, job_ids, cluster_cmd, delay=120):
        """
        Wait on a set of job IDs.
        """
        pass    
    
    def launch_job(self, cluster_cmd, cmd_name):
        pass