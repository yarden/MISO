#!/bin/bash

sudo rm -rf pysplicing/build
sudo python pysplicing/setup.py clean
make
cd pysplicing
sudo python setup.py install
