##
## Parse plotting configuration files for sashimi_plot
##
import ConfigParser

def get_default_settings():
    """
    Return default settings for sashimi_plot.
    """
    settings = {"intron_scale": 30,
                "exon_scale": 1,
                "logged": False,
                "ymax": None,
                "show_posteriors": True,
                "number_junctions": True,
                "posterior_bins": 40,
                "gene_posterior_ratio": 5,
                "resolution": .5,
                "fig_width": 8.5,
                "fig_height": 11,
                "bar_posteriors": False,
                "junction_log_base": 10.,
                "reverse_minus": False,
                "font_size": 6}
    return settings

def parse_plot_settings(settings,
                        FLOAT_PARAMS=["intron_scale", "exon_scale", "ymax",
                                      "resolution", "fig_width", "fig_height",
                                      "font_size", "junction_log_base"],
                        INT_PARAMS=["posterior_bins", "gene_posterior_ratio"],
                        BOOL_PARAMS=["logged", "show_posteriors", "number_junctions",
                                     "reverse_minus", "bar_posteriors"]):
    """
    Populate a settings dictionary with the plotting parameters, parsed
    as the right datatype.
    """
    for section in config.sections():
        for option in config.options(section):
            if option in FLOAT_PARAMS:
                settings[option] = config.getfloat(section, option)
            elif option in INT_PARAMS:
                settings[option] = config.getint(section, option)
            elif option in BOOL_PARAMS:
                settings[option] = config.getboolean(section, option)
            else:
                settings[option] = config.get(section, option)
    settings["bam_files"] = json.loads(settings["bam_files"])
    settings["miso_files"] = json.loads(settings["miso_files"])
    
    if "colors" in settings:
        colors = json.loads(settings["colors"])
    else:
        colors = [None for x in settings["bam_files"]]
    if "bam_prefix" in settings:
        bam_files = [os.path.join(settings["bam_prefix"], x) \
                     for x in settings["bam_files"]]
    else:
        bam_files = settings["bam_files"]
    if "miso_prefix" in settings:
        miso_files = get_miso_output_files(event, chrom, settings)
    else:
        miso_files = settings["miso_files"]
    if "coverages" in settings:
        coverages = json.loads(settings["coverages"])
        coverages = map(float, coverages)
        # Normalize coverages per M
        coverages = [x / 1e6  for x in coverages]
    else:
        coverages = [1 for x in settings["bam_files"]]
    return settings
