#!/bin/bash

# convert file
#python sam_to_bam.py --convert test-data/sam-data/c2c12.Atp2b1.sam indexed-sam

# Run MISO
#python run_events_analysis.py --compute-genes-psi indexed-gff/ indexed-sam/ --output-dir test-output/ \
#--read-len 36

### C-version
rm ~/Projects/fastmiso/MISO/test-output/chr10/ENSMUSG00000019943.miso

# Single-end
#echo "Testing single-end"
#python /Users/yarden/Projects/fastmiso/MISO/run_miso.py --compute-gene-psi "ENSMUSG00000019943" "/Users/yarden/Projects/fastmiso/MISO/indexed-gff/chr10/ENSMUSG00000019943.pickle" indexed-sam/c2c12.Atp2b1.sorted.bam /Users/yarden/Projects/fastmiso/MISO/test-output --read-len 36 --settings-filename /Users/yarden/Projects/fastmiso/MISO/settings/miso_settings.txt

rm ~/Projects/fastmiso/MISO/test-output/chr10/ENSMUSG00000019943.miso

# Paired-end
echo "Testing paired-end"
python /Users/yarden/Projects/fastmiso/MISO/run_miso.py --compute-gene-psi "ENSMUSG00000019943" "/Users/yarden/Projects/fastmiso/MISO/indexed-gff/chr10/ENSMUSG00000019943.pickle" indexed-sam/c2c12.Atp2b1.sorted.bam /Users/yarden/Projects/fastmiso/MISO/test-output --paired-end 250 30 --read-len 36 --settings-filename /Users/yarden/Projects/fastmiso/MISO/settings/miso_settings.txt


##
## Py version for comparison
##
#python /home/yarden/MISO/run_miso.py --compute-gene-psi "ENSMUSG00000019943" "/home/yarden/MISO/gff-events/mm9/genes/Atp2b1/indexed/chr10/ENSMUSG00000019943.pickle" /home/yarden/MISO/test-output/sam-output/c2c12.Atp2b1.sorted.bam /home/yarden/MISO/test-output/gene-psi-output --read-len 36 --settings-filename /home/yarden/MISO/settings/miso_settings.txt