##
## Settings with relevant directories
##

import misopy
from misopy.parse_csv import *
import ConfigParser
import os
import sys

miso_path = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))

class Settings(object):
    @classmethod
    def load(cls, path):
        """
        Reads in settings from a ConfigParser formatted file
        ignores section headers, so make sure each option is unique in the file
        returns a dictionary with all the options mapped to their values.
        """
        config = ConfigParser.ConfigParser()

        if path != None:
            cls.settings_path = path
        else:
            # Use default settings file if none was given
            cls.settings_path = os.path.join(miso_path,
                                             "settings",
                                             "miso_settings.txt")

        print "Loading settings from: %s" %(cls.settings_path)
        cls.parsed_settings = config.read(cls.settings_path)

        cls.global_settings = {}

        print "Settings: "
        for section in config.sections():
            for option in config.options(section):
                # Load cluster options as strings, without attempting to evaluate them
                # Avoids misinterpretation of words like "long" as a data type
                if section == "cluster":
                    cls.global_settings[option] = str(config.get(section, option))
                else:
                    cls.global_settings[option] = tryEval(config.get(section, option))
                print "  ", option, cls.global_settings[option]

        # Set directory paths specific to pipeline
        if 'pipeline_results_dir' in cls.global_settings:
            cls.global_settings['analysis_dir'] = os.path.join(cls.global_settings['pipeline_results_dir'],
                                                               'analysis')
            cls.global_settings['rna_events_dir'] = os.path.join(cls.global_settings['analysis_dir'],
                                                                 'rna_events')

    @classmethod
    def get_sampler_params(cls):
        """
        Return sampler parameters.
        """
        param_names = ['burn_in', 'lag', 'num_iters']
        opt_param_names = ['num_chains']

        # Default number of chains is 6
        sampler_params = {'num_chains': 6}

        for name in param_names:
            if name not in cls.global_settings:
                print "settings: ", cls.global_settings
                raise Exception, "Error: need %s parameter to be set in settings file." \
                      %(name)
            sampler_params[name] = cls.global_settings[name]
        # Record optional parameters
        for name in opt_param_names:
            if name in cls.global_settings:
                sampler_params[name] = cls.global_settings[name]
        return sampler_params
    

    @classmethod
    def get_cluster_command(cls):
        """
        Return the name of the command to use for cluster submission
        (e.g. 'qsub')
        """
        if 'cluster_command' in cls.global_settings:
            return cls.global_settings['cluster_command']
        else:
            return None

    @classmethod
    def get_long_queue_name(cls):
        """
        Return the name of the long queue (for long jobs.)
        """
        print "cls.global_settings: ", cls.global_settings
        if 'long_queue_name' in cls.global_settings:
            return cls.global_settings['long_queue_name']
        else:
            return None
        
    @classmethod
    def get_short_queue_name(cls):
        """
        Return the name of the short queue (for short jobs.)
        """
        if 'long_queue_name' in cls.global_settings:
            return cls.global_settings['short_queue_name']
        else:
            return None

        

    @classmethod
    def get_min_event_reads(cls):
        """
        Return minimum number of reads an event should have.
        """
        min_event_reads = cls.global_settings["min_event_reads"]
        return min_event_reads
        
        
    @classmethod
    def get_counts_dir(cls, event_type):
        """
        Return counts directory for given event type.
        """
        if 'rna_events_dir' in cls.global_settings:
            return os.path.join(cls.global_settings['rna_events_dir'],
                                event_type)
        return None

    @classmethod
    def get_counts_filename(cls, sample_label, event_type):
        """
        Return counts filename for a given sample and its type.
        """
        return os.path.join(cls.get_counts_dir(event_type),
                            '%s.counts' %(sample_label))

    @classmethod
    def get_filters(cls, event_type):
        pass
        
    @classmethod
    def get(cls):
        return cls.global_settings

    @classmethod
    def get_miso_exec(cls):
        if 'MISO_SHELL_EXEC' in os.environ:
            return os.environ['MISO_SHELL_EXEC']
        return sys.executable
    

def load_settings(settings_filename):
    Settings.load(settings_filename)
    return Settings.get()


def main():
    settings_filename = 'settings/miso_settings.txt'
    Settings.load(settings_filename)
    print Settings.get()
    print Settings.get_counts_dir('TandemUTR')

if __name__ == '__main__':
    main()
