'''
misopy.cluster

Cluster factory for different cluster types

@author: Aaron Kitzmiller
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
@contact: aaron_kitzmiller@harvard.edu
'''

import os, subprocess, traceback, time

from misopy.settings import load_settings, Settings
from misopy import misc_utils


def getClusterEngine(cluster_type,settings_fname):
    '''
    Returns the correct cluster engine
    '''
    ce = None
    if cluster_type == 'slurm':
        ce = SlurmClusterEngine(settings_fname)
    elif cluster_type == 'lsf':
        ce = LsfClusterEngine(settings_fname)
    elif cluster_type == 'sge':
        ce = SgeClusterEngine(settings_fname)
    else:
        raise Exception('Unknown cluster type %s' % cluster_type)
        
    return ce
    
class AbstractClusterEngine(object):
    '''
    Base class for cluster engines
    '''
    def __init__(self,settings_filename):
        '''
        Load settings
        '''
        self.settings = load_settings(settings_filename)
        
    def wait_on_jobs(self,job_ids, cluster_cmd,
                     delay=120.0):
        """
        Wait on a set of job IDs.
        """
        if len(job_ids) == 0:
            return
        num_jobs = len(job_ids)
        print "Waiting on a set of %d jobs..." %(num_jobs)
        curr_time = time.strftime("%x, %X")
        t_start = time.time()
        print "  - Starting to wait at %s" %(curr_time)
        completed_jobs = {}
        for job_id in job_ids:
            if job_id in completed_jobs:
                continue
            self.wait_on_job(job_id, cluster_cmd)
            print "  - Job ", job_id, " completed."
            completed_jobs[job_id] = True
        curr_time = time.strftime("%x, %X")
        t_end = time.time()
        print "Jobs completed at %s" %(curr_time)
        duration = ((t_end - t_start) / 60.) / 60.
        print "  - Took %.2f hours." %(duration)
        
        
    def wait_on_job(self, job_id, delay):
        
        raise Exception('Must implement wait on job')
    
    
class LsfClusterEngine(AbstractClusterEngine):
    '''
    Run jobs on an LSF cluster
    '''
    
    def make_bash_script(self,filename, cmd, crate_dir=None):
        """
        Make an executable bash script out of the given command.
        """
    #    os.system('ls %s' %(filename))
        if crate_dir == None:
            crate_dir = \
                os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
        f = open(filename, 'w')
        f.write("#!/bin/bash\n")
        f.write("export PATH=$PATH:%s\n" %(crate_dir))
        f.write("source ~/.bash_profile\n")
        f.write("cd %s\n" %(crate_dir))
        #write_cluster_preface(f)
        f.write(cmd + "\n")
        f.close()
        os.system('chmod +x \"%s\"' %(filename))
        
    
    def run_on_cluster(self, cmd, job_name, cluster_output_dir,
                       cluster_scripts_dir=None,
                       queue_type=None):
        '''
        Composes job script and launches job
        '''
        print "Submitting job: %s" %(job_name)
        queue_name = None
    
        # Load command name from settings file
        cmd_name = self.settings.get_cluster_command()
    
        if queue_type == "long":
            queue_name = self.settings.get_long_queue_name()
        elif queue_type == "short":
            queue_name = self.settings.get_short_queue_name()
        else:
            print "Warning: Unknown queue type: %s" %(queue_type)
            queue_name = queue_type
        
        if queue_type is None:
            print "  - queue type: unspecified"
        else:
            print "  - queue type: %s" %(queue_type)
        if queue_name is None:
            print " - queue name unspecified"
        else:
            print " - queue name: %s" %(queue_name)
            
        misc_utils.make_dir(cluster_output_dir)
        if cluster_scripts_dir == None:
            cluster_scripts_dir = os.path.join(cluster_output_dir,
                                               'cluster_scripts')
            misc_utils.make_dir(cluster_scripts_dir)
        scripts_output_dir = os.path.join(cluster_output_dir,
                                          'scripts_output')
        misc_utils.make_dir(scripts_output_dir)
        scripts_output_dir = os.path.abspath(scripts_output_dir)
        cluster_call = 'bsub -o \"%s\" -e \"%s\"' %(scripts_output_dir,
                                                  scripts_output_dir)
        # Add queue type if given one
        if queue_name != None:
            cluster_call += ' -q \"%s\"' %(queue_name)
            
        script_name = os.path.join(cluster_scripts_dir,
                                         '%s_time_%s.sh' \
                                         %(job_name,
                                           time.strftime("%m-%d-%y_%H_%M_%S")))
        self.make_bash_script(script_name, cmd)
        cluster_cmd = cluster_call + ' \"%s\"' %(script_name)
        job_id = self.launch_job(cluster_cmd)
        return job_id
    
    
    
    def wait_on_job(self, job_id, delay=60):
        # Handle bsub
        while True:
            output = subprocess.Popen("bjobs %i" %(job_id),
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
            if len(output[0]) > 0:
                status = output[0].split()[10]
                if status == "DONE":
                    break
            else:
                # No jobs available
                break
            time.sleep(delay)
        time.sleep(delay)
        
        
    
    def launch_job(self, cluster_cmd):
        """
        Execute cluster_cmd and return its job ID if
        it can be fetched.
        """
        print "Executing: %s" %(cluster_cmd)
        proc = subprocess.Popen(cluster_cmd, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE)
        # Read the job ID if it's a known cluster
        # submission system
        output = proc.communicate()
        job_id = None
        if "is submitted to" in output[0]:
            job_id = int(output[0].strip().split()[1][1:-1])                
        return job_id
    


class SgeClusterEngine(AbstractClusterEngine):
    '''
    Run jobs on an SGE cluster
    '''
            
    def make_bash_script(self,filename, cmd, crate_dir=None):
        """
        Make an executable bash script out of the given command.
        """
    #    os.system('ls %s' %(filename))
        if crate_dir == None:
            crate_dir = \
                os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
        f = open(filename, 'w')
        f.write("#!/bin/bash\n")
        f.write("export PATH=$PATH:%s\n" %(crate_dir))
        f.write("source ~/.bash_profile\n")
        f.write("cd %s\n" %(crate_dir))
        #write_cluster_preface(f)
        f.write(cmd + "\n")
        f.close()
        os.system('chmod +x \"%s\"' %(filename))
        
    
    def run_on_cluster(self, cmd, job_name, cluster_output_dir,
                       cluster_scripts_dir=None,
                       queue_type=None):
        '''
        Composes job script and launches job
        '''
        print "Submitting job: %s" %(job_name)
        queue_name = None
    
        # Load command name from settings file
        cmd_name = self.settings.get_cluster_command()
    
        if queue_type == "long":
            queue_name = self.settings.get_long_queue_name()
        elif queue_type == "short":
            queue_name = self.settings.get_short_queue_name()
        else:
            print "Warning: Unknown queue type: %s" %(queue_type)
            queue_name = queue_type
        
        if queue_type is None:
            print "  - queue type: unspecified"
        else:
            print "  - queue type: %s" %(queue_type)
        if queue_name is None:
            print " - queue name unspecified"
        else:
            print " - queue name: %s" %(queue_name)
            
        misc_utils.make_dir(cluster_output_dir)
        if cluster_scripts_dir == None:
            cluster_scripts_dir = os.path.join(cluster_output_dir,
                                               'cluster_scripts')
            misc_utils.make_dir(cluster_scripts_dir)
        scripts_output_dir = os.path.join(cluster_output_dir,
                                          'scripts_output')
        misc_utils.make_dir(scripts_output_dir)
        scripts_output_dir = os.path.abspath(scripts_output_dir)
        cluster_call = 'qsub -o \"%s\" -e \"%s\"' %(scripts_output_dir,
                                                  scripts_output_dir)
        # Add queue type if given one
        if queue_name != None:
            cluster_call += ' -q \"%s\"' %(queue_name)
            
        script_name = os.path.join(cluster_scripts_dir,
                                         '%s_time_%s.sh' \
                                         %(job_name,
                                           time.strftime("%m-%d-%y_%H_%M_%S")))
        self.make_bash_script(script_name, cmd)
        cluster_cmd = cluster_call + ' \"%s\"' %(script_name)
        job_id = self.launch_job(cluster_cmd)
        return job_id
    
    
    
    def wait_on_job(self, job_id, delay=60):
        # Handle qsub
        while True:
            output = \
                subprocess.Popen("qstat %i" %(job_id),
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()
            if "Unknown Job" in output[1]:
                break
            time.sleep(delay)
        time.sleep(delay)
        
        
    
    def launch_job(self, cluster_cmd):
        """
        Execute cluster_cmd and return its job ID if
        it can be fetched.
        """
        print "Executing: %s" %(cluster_cmd)
        proc = subprocess.Popen(cluster_cmd, shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                stdin=subprocess.PIPE)
        # Read the job ID if it's a known cluster
        # submission system
        output = proc.communicate()
        job_id = None
        if "." in output[0][:-1] and ">" not in output[0]:
            job_id = int(output[0].split(".")[0])
        return job_id



class SlurmClusterEngine(AbstractClusterEngine):
    '''
    Run jobs on a Slurm cluster
    '''
    
    def __init__(self,settings_filename):
        '''
        Get the slurm job template file and load it
        '''
        super(SlurmClusterEngine,self).__init__(settings_filename)
        
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
        
        
        
    def make_bash_script(self,script_name,cmd):
        '''
        Use the template to write out a sbatch submission script
        '''
        scripttxt = self.template.format(cmd=cmd)
        with open(script_name,'w') as script:
            script.write(scripttxt + '\n')
        
        
    
    def run_on_cluster(self, cmd, job_name, cluster_output_dir,
                       cluster_scripts_dir=None,
                       queue_type=None,
                       settings_fname=None):
        '''
        Composes job script and launches job
        '''
        
        misc_utils.make_dir(cluster_output_dir)
        if cluster_scripts_dir == None:
            cluster_scripts_dir = os.path.join(cluster_output_dir,
                                               'cluster_scripts')
            misc_utils.make_dir(cluster_scripts_dir)
            
        scripts_output_dir = os.path.join(cluster_output_dir,
                                          'scripts_output')
        misc_utils.make_dir(scripts_output_dir)
        scripts_output_dir = os.path.abspath(scripts_output_dir)
        cluster_call = 'sbatch -D \"%s\"' %(scripts_output_dir)
        
        script_name = os.path.join(cluster_scripts_dir,
                                         '%s_time_%s.sh' \
                                         %(job_name,
                                           time.strftime("%m-%d-%y_%H_%M_%S")))
        self.make_bash_script(script_name, cmd)
        cluster_cmd = cluster_call + ' \"%s\"' %(script_name)
        job_id = self.launch_job(cluster_cmd)
        return job_id
    
    
    
    def wait_on_job(self, job_id, delay=10):
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
                    print state
                    
            time.sleep(10)         
        
        
    
    def launch_job(self, cluster_cmd):
        '''
        Runs cluster command and returns the job id
        '''
        print 'Running command %s' % cluster_cmd
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
    