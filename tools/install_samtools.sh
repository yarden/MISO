#! /bin/bash

VER=1.2
wget https://github.com/samtools/samtools/releases/download/${VER}/samtools-${VER}.tar.bz2
tar xjf samtools-${VER}.tar.bz2
cd samtools-${VER}
make
make install
