import os
import pstats
from run_events_analysis import *

def start_miso():
    gff_dir = os.path.abspath(os.path.join("gff-events", "mm9", "genes",
                                           "indexed"))

    if not os.path.isdir(gff_dir):
        print "%s is not a valid indexed GFF directory."
        return
    
    read_len = 36
    bam_filename = os.path.abspath(os.path.join("test-output",
                                                "sam-output",
                                                "c2c12.Atp2b1.sorted.bam"))
    compute_all_genes_psi(gff_dir, bam_filename,
                          read_len, "profiler-test")

    
def main():
    output_file = "profile"
    print "Profiling MISO... outputting to: %s" %(output_file)
    import cProfile as profile
    profile.run('start_miso()', output_file)
    p = pstats.Stats(output_file)
    print "name: "
    print p.sort_stats('name')
    print "all stats: "
    p.print_stats()
    print "cumulative (top 10): "
    p.sort_stats('cumulative').print_stats(10)
    

if __name__ == '__main__':
    main()
    
