import os, sys, operator, subprocess
try:
    import simplejson as json
except:
    import json
import ConfigParser
import sam_utils
import GFF as gff_utils
from pylab import *
from matplotlib.patches import PathPatch 
from matplotlib.path import Path


# Plot MISO events using BAM files and posterior distribution files.
# If comparison files are available, plot bayes factors too.
def plot_density_single(tx_start, tx_end, gene_obj, mRNAs, strand, graphcoords,\
    graphToGene, bam_filename, axvar, paired_end=False, intron_scale=30,\
    exon_scale=4, color='r', ymax=None, logged=False, coverage=1,\
    number_junctions=True, resolution=.5, showXaxis=True, showYaxis=True,\
    showYlabel=True, font_size=6, junction_log_base=10):

    bamfile = sam_utils.load_bam_reads(bam_filename) 

    gene_reads = sam_utils.fetch_bam_reads_in_gene(bamfile, gene_obj.chrom,\
        tx_start, tx_end, gene_obj)
    reads, num_raw_reads = sam_utils.sam_parse_reads(gene_reads,\
        paired_end=paired_end)
    wiggle, jxns = readsToWiggle(reads, tx_start, tx_end)

    wiggle = 1e3 * wiggle / coverage
 
    if logged:
        wiggle = log10(wiggle + 1)
    
    maxheight = max(wiggle)
    if ymax is None:
        ymax = 1.1 * maxheight
    else:
        ymax = ymax
    ymin = -.5 * ymax 

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []
    prevx = 0
    tmpval = []
    for i in range(len(graphcoords)):
        tmpval.append(wiggle[i])
        if abs(graphcoords[i] - prevx) > resolution:
            compressed_wiggle.append(mean(tmpval))
            compressed_x.append(prevx)
            prevx = graphcoords[i]
            tmpval = []

    fill_between(compressed_x, compressed_wiggle,\
        y2=0, color=color, lw=0)
   
    sslists = []
    for mRNA in mRNAs:
        tmp = []
        for s, e in mRNA:
            tmp.extend([s, e])
        sslists.append(tmp)

    for jxn in jxns:
        leftss, rightss = map(int, jxn.split(":"))

        ss1, ss2 = [graphcoords[leftss - tx_start - 1],\
            graphcoords[rightss - tx_start]]

        mid = (ss1 + ss2) / 2
        h = -3 * ymin / 4
   
        numisoforms = 0
        for i in range(len(mRNAs)):
            if leftss in sslists[i] and \
                rightss in sslists[i]:
                numisoforms += 1
        if numisoforms > 0:
            if numisoforms % 2 == 0: # put on bottom 
                pts = [(ss1, 0), (ss1, -h), (ss2, -h), (ss2, 0)]
                midpt = cubic_bezier(pts, .5)
            else:                         # put on top 
                leftdens = wiggle[leftss - tx_start - 1]
                rightdens = wiggle[rightss - tx_start]

                pts = [(ss1, leftdens),\
                    (ss1, leftdens + h),\
                    (ss2, rightdens + h),\
                    (ss2, rightdens)]
                midpt = cubic_bezier(pts, .5)

            if number_junctions:
                text(midpt[0], midpt[1], '%s'%(jxns[jxn]),\
                    fontsize=6, ha='center', va='center', backgroundcolor='w')

            a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
            p = PathPatch(a, ec=color, lw=log(jxns[jxn] + 1) /\
                log(junction_log_base), fc='none')
            axvar.add_patch(p) 

    # Format plot.
    nyticks = 4
    nxticks = 5
    ylim(ymin, ymax)

    axvar.spines['left'].set_bounds(0, ymax)
    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')

    if showXaxis:
        axvar.xaxis.set_ticks_position('bottom')
        xlabel('Genomic coordinate (%s), "%s" strand'%(gene_obj.chrom,\
            strand), fontsize=font_size)
        xticks(linspace(0, max(graphcoords), nxticks),\
            [graphToGene[int(x)] for x in \
            linspace(0, max(graphcoords), nxticks)], fontsize=font_size)
    else:
        axvar.spines['bottom'].set_color('none')
        xticks([])

    if showYlabel:
        if logged:
            ylabel('RPKM $(log_{10})$', fontsize=font_size)
        else:
            ylabel('RPKM', fontsize=font_size)
    if showYaxis:
        axvar.yaxis.set_ticks_position('left')
        yticks(linspace(0, ymax, nyticks), ['%d'%(x) for x in \
            linspace(0, ymax, nyticks)], fontsize=font_size)
    else:
        axvar.spines['left'].set_color('none')
        yticks([])
    
    text(max(graphcoords), ymax, os.path.basename(bam_filename).split(".")[0],\
        fontsize=font_size, va='bottom', ha='right', color=color)
    xlim(0, max(graphcoords))


# Plot density for a series of bam files.
def plot_density(pickle_filename, event, bam_files, miso_files, out_f,\
    intron_scale=30, exon_scale=1, gene_posterior_ratio=5, posterior_bins=40,\
    colors=None, ymax=None, logged=False, show_posteriors=True, coverages=None,\
    number_junctions=True, resolution=.5, fig_width=8.5, fig_height=11,\
    font_size=6, junction_log_base=10, reverse_minus=False):

    # Parse gene pickle and get scaling information.
    tx_start, tx_end, exon_starts, exon_ends, gene_obj, mRNAs, strand = \
        parseGene(pickle_filename, event)
    graphcoords, graphToGene = getScaling(tx_start, tx_end, strand,\
        exon_starts, exon_ends, intron_scale, exon_scale, reverse_minus)

    nfiles = len(bam_files)
    figure(figsize=(fig_width, fig_height))
    suptitle(event, fontsize=10)
    
    for i in range(nfiles):
        if colors is not None:
            color = colors[i]
        else:
            color = None
        if coverages is not None:
            coverage = coverages[i]
        else:
            coverage = 1
        if i < nfiles - 1:
            showXaxis = False 
        else:
            showXaxis = True 

        bam_file = os.path.expanduser(bam_files[i])
        miso_file = os.path.expanduser(miso_files[i])
        ax1 = subplot2grid((nfiles + 2, gene_posterior_ratio), (i, 0),\
            colspan=gene_posterior_ratio - 1)
        plot_density_single(tx_start, tx_end, gene_obj, mRNAs, strand,\
            graphcoords, graphToGene, bam_file, ax1, paired_end=False,\
            intron_scale=intron_scale, exon_scale=exon_scale, color=color,\
            ymax=ymax, logged=logged, coverage=coverage,\
            number_junctions=number_junctions, resolution=resolution,\
            showXaxis=showXaxis, showYlabel=True, font_size=font_size,\
            junction_log_base=junction_log_base)

        if show_posteriors:
            try:
                ax2 = subplot2grid((nfiles + 2, gene_posterior_ratio),\
                    (i, gene_posterior_ratio - 1))

                if not os.path.isfile(miso_file):
                    print "Warning: MISO file %s not found" %(miso_file)
                    raise Exception

                print "Loading MISO file: %s" %(miso_file)
                plot_posterior_single(miso_file, ax2, posterior_bins,\
                    showXaxis=showXaxis, showYlabel=False, font_size=font_size)
            except:
                box(on=False)
                xticks([])
                yticks([])
                print "Posterior plot failed."

    # Draw gene structure
    ax = subplot2grid((nfiles + 2, gene_posterior_ratio), (nfiles, 0),\
        colspan=gene_posterior_ratio - 1, rowspan=2)
    plot_mRNAs(tx_start, mRNAs, strand, graphcoords, ax)

    subplots_adjust(hspace=.1, wspace=.7)
    savefig(out_f)


# Parse the pickled gene.
def parseGene(pickle_filename, event):
    if not os.path.isfile(pickle_filename):
        raise Exception, "Error: no filename %s" %(pickle_filename)
    gff_genes = gff_utils.load_indexed_gff_file(pickle_filename)

    if gff_genes == None:
        raise Exception, "Error: could not load genes from %s" \
              %(pickle_filename)

    exon_starts = []
    exon_ends = []
    mRNAs = []
    for gene_id, gene_info in gff_genes.iteritems():
        if event == gene_id:
            gene_obj = gene_info['gene_object']
            gene_hierarchy = gene_info['hierarchy']
            tx_start, tx_end = gff_utils.get_inclusive_txn_bounds(\
                gene_hierarchy[gene_id])

            for mRNA_id, mRNA_info in gene_hierarchy[gene_id]['mRNAs'].iteritems():
                mRNA = []
                for exon_id, exon_info in gene_hierarchy[gene_id]['mRNAs']\
                    [mRNA_id]['exons'].\
                    iteritems():

                    exon_rec = gene_hierarchy[gene_id]['mRNAs']\
                        [mRNA_id]['exons'][exon_id]['record']
                    strand = exon_rec.strand
                    exon_starts.append(exon_rec.start)
                    exon_ends.append(exon_rec.end)
                    mRNA.append(sorted([exon_rec.start, exon_rec.end]))

                mRNAs.append(mRNA)
            break

    mRNAs.sort(key=len)
    return tx_start, tx_end, exon_starts, exon_ends, gene_obj, mRNAs, strand

# Compute the scaling factor across various genic regions.
def getScaling(tx_start, tx_end, strand, exon_starts, exon_ends,\
    intron_scale, exon_scale, reverse_minus):
   
    exoncoords = zeros((tx_end - tx_start + 1))
    for i in range(len(exon_starts)):
        exoncoords[exon_starts[i] - tx_start : exon_ends[i] - tx_start] = 1

    graphToGene = {}
    graphcoords = zeros((tx_end - tx_start + 1), dtype='f')
    x = 0
    if strand == '+' or not reverse_minus:
        for i in range(tx_end - tx_start + 1):
            graphcoords[i] = x
            graphToGene[int(x)] = i + tx_start
            if exoncoords[i] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    else:
        for i in range(tx_end - tx_start + 1):
            graphcoords[-(i + 1)] = x
            graphToGene[int(x)] = tx_end - i + 1
            if exoncoords[-(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale
    
    return graphcoords, graphToGene


# Get wiggle and junction densities from reads.
def readsToWiggle(reads, tx_start, tx_end):

    read_positions, read_cigars = reads
    wiggle = zeros((tx_end - tx_start + 1), dtype='f')
    jxns = {}
    for i in range(len(read_positions)):
        pos, cigar = [read_positions[i], read_cigars[i]]
        if "N" not in cigar:
            rlen = int(cigar[:-1])
            s = max([pos - tx_start, 0])
            e = min([pos - tx_start + rlen, len(wiggle) - 1])
            wiggle[s : e] += 1. / rlen
        else:
            left, right = cigar.split("N")
            left, middle = map(int, left.split("M"))
            right = int(right[:-1])
            rlen = left + right
            s1 = pos - tx_start
            e1 = pos - tx_start + left
            s2 = pos + left + middle - tx_start
            e2 = pos + left + middle + right - tx_start

            # Include read coverage from adjacent junctions.
            if (e1 >= 0 and e1 < len(wiggle)) or (s1 >= 0 and s1 < len(wiggle)):
                wiggle[max([s1, 0]) : min([e1, len(wiggle)])] += 1. / rlen
            if (e2 >= 0 and e2 < len(wiggle)) or (s2 >= 0 and s2 < len(wiggle)):
                wiggle[max([s2, 0]) : min([e2, len(wiggle)])] += 1. / rlen

            # Plot a junction if both splice sites are within locus.
            leftss = pos + left
            rightss = pos + left + middle + 1
            if leftss - tx_start >= 0 and leftss - tx_start < len(wiggle) \
                and rightss - tx_start >= 0 and rightss - tx_start < \
                len(wiggle): 

                jxn = ":".join(map(str, [leftss, rightss]))
                try:
                    jxns[jxn] += 1 
                except:
                    jxns[jxn] = 1 

    return wiggle, jxns


# Draw the gene structure.
def plot_mRNAs(tx_start, mRNAs, strand, graphcoords, axvar):
  
    yloc = 0 
    exonwidth = .3
    narrows = 50

    for mRNA in mRNAs:
        for s, e in mRNA:
            s = s - tx_start
            e = e - tx_start
            x = [graphcoords[s], graphcoords[e], graphcoords[e], graphcoords[s]]
            y = [yloc - exonwidth / 2, yloc - exonwidth / 2,\
                yloc + exonwidth / 2, yloc + exonwidth / 2]
            fill(x, y, 'k', lw=.5, zorder=20)

        # Draw intron.
        axhline(yloc, color='k', lw=.5)

        # Draw intron arrows.
        spread = .2 * max(graphcoords) / narrows
        for i in range(narrows):
            loc = float(i) * max(graphcoords) / narrows
            x = [loc - spread, loc, loc - spread]
            y = [yloc - exonwidth / 5, yloc, yloc + exonwidth / 5]
            plot(x, y, lw=.5, color='k')

        yloc += 1 

    xlim(0, max(graphcoords)) 
    ylim(-.5, len(mRNAs) + .5)
    box(on=False)
    xticks([])
    yticks([]) 



# Plot a posterior probability distribution for a MISO event
def plot_posterior_single(miso_f, axvar, posterior_bins,\
    showXaxis=True, showYaxis=True, showYlabel=True, font_size=6):
  
    posterior_bins = int(posterior_bins) 
    psis = [] 
    for line in open(miso_f):
        if not line.startswith("#") and not line.startswith("sampled"):
            psi, logodds = line.strip().split("\t")
            psis.append(float(psi.split(",")[0]))
  
    ci = .95 
    alpha = 1 - ci
    lidx = int(round((alpha / 2) * len(psis)) - 1)
    # the upper bound is the (1-alpha/2)*n nth smallest sample, where n is
    # the number of samples
    hidx = int(round((1 - alpha / 2) * len(psis)) - 1)
    psis.sort()
    clow, chigh = [psis[lidx], psis[hidx]]
   
    y, x, p = hist(psis, linspace(0, 1, posterior_bins),\
        normed=True, facecolor='k', edgecolor='w', lw=.2) 
    axvline(clow, ymin=.33, linestyle='--', dashes=(1, 1), color='#CCCCCC', lw=.5)
    axvline(chigh, ymin=.33, linestyle='--', dashes=(1, 1), color='#CCCCCC', lw=.5)
    axvline(mean(psis), ymin=.33, color='r')

    ymax = max(y) * 1.5
    ymin = -.5 * ymax

    #text(1, ymax,\
    #    "$\Psi$ = %.2f"%(mean(psis)), va='top', ha='left',\
    #    fontsize=font_size)
    
    text(1, ymax,\
        "$\Psi$ = %.2f\n$\Psi_{0.05}$ = %.2f\n$\Psi_{0.95}$ = %.2f"%\
        (mean(psis), clow, chigh), fontsize=font_size, va='top', ha='left')

    ylim(ymin, ymax)
    axvar.spines['left'].set_bounds(0, ymax)
    axvar.spines['right'].set_color('none')
    axvar.spines['top'].set_color('none')
    axvar.spines['bottom'].set_position(('data', 0)) 
    axvar.xaxis.set_ticks_position('bottom')
    axvar.yaxis.set_ticks_position('left')
    
    if showXaxis:
        xticks(fontsize=font_size)
        xlabel("MISO $\Psi$", fontsize=font_size)
    else:
        [label.set_visible(False) for label in axvar.get_xticklabels()]

    nyticks = 4
    if showYaxis:
        yticks(linspace(0, ymax, nyticks),\
            ["%d"%(y) for y in linspace(0, ymax, nyticks)],\
            fontsize=font_size)
    else:
        yticks([])
    if showYlabel:
        ylabel("Frequency", fontsize=font_size, ha='right', va='center')


# Get points in a cubic bezier.
def cubic_bezier(pts, t):
    p0, p1, p2, p3 = pts
    p0 = array(p0)
    p1 = array(p1)
    p2 = array(p2)
    p3 = array(p3)
    return p0 * (1 - t)**3 + 3 * t * p1 * (1 - t) ** 2 + \
        3 * t**2 * (1 - t) * p2 + t**3 * p3

    
# A wrapper to allow reading from a file.
def plot_density_from_file(pickle_filename, event, settings_f, out_f):

    config = ConfigParser.ConfigParser()
    config.read(settings_f) 
   
    settings = {"intron_scale": 30,\
        "exon_scale": 1,\
        "logged": False,\
        "ymax": None,\
        "show_posteriors": True,\
        "number_junctions": True,\
        "posterior_bins": 40,\
        "gene_posterior_ratio": 5,\
        "resolution": .5,\
        "fig_width": 8.5,\
        "fig_height": 11,\
        "junction_log_base": 10.,\
        "reverse_minus": False,\
        "font_size": 6} 
    for section in config.sections():
        for option in config.options(section):
            if option in ["intron_scale", "exon_scale", "ymax", "resolution",\
                "fig_width", "fig_height", "font_size", "junction_log_base"]:
                settings[option] = config.getfloat(section, option)
            elif option in ["posterior_bins", "gene_posterior_ratio"]:
                settings[option] = config.getint(section, option)
            elif option in ["logged", "show_posteriors", "number_junctions",\
                "reverse_minus"]:
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
    chrom = event.split(":")[0]
    strand = event.split(":")[-1]
    if "miso_prefix" in settings:
        miso_files = [os.path.join(settings["miso_prefix"], x, chrom,\
            event) + ".miso" for x in settings["miso_files"]]
    else:
        miso_files = settings["miso_files"]
    if "coverages" in settings:
        coverages = json.loads(settings["coverages"])
        coverages = map(float, coverages)
        coverages = [x / 1e6  for x in coverages]
    else:
        coverages = [1 for x in settings["bam_files"]]

    plot_density(pickle_filename, event, bam_files, miso_files, out_f,\
        intron_scale=settings["intron_scale"],\
        exon_scale=settings["exon_scale"],\
        gene_posterior_ratio=settings["gene_posterior_ratio"],\
        posterior_bins=settings["posterior_bins"],\
        show_posteriors=settings["show_posteriors"],\
        logged=settings["logged"],\
        colors=colors, ymax=settings["ymax"], coverages=coverages,\
        number_junctions=settings["number_junctions"],\
        resolution=settings["resolution"],\
        fig_width=settings["fig_width"], fig_height=settings["fig_height"],\
        font_size=settings["font_size"],\
        junction_log_base=settings["junction_log_base"],\
        reverse_minus=settings["reverse_minus"])


