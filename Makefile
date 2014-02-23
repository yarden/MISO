bam_example: 
	g++ -I ./lib/bamtools-master/include/ -L ./lib/bamtools-master/src/ -o bam_example src/bam_example.cpp -lz -lbamtools

clean:
	rm bam_example