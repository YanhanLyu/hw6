sample.fasta - An example Fasta file that can be used to simulate reads

reads.txt - Output file for reads

simulate.sh will take in the following parameters (in this order):
• FASTA sequence file (string) - File containing a DNA sequence.
• Coverage (integer) - The coverage of the simulated sequencing experiment.
• Read length (integer) - The length of reads to simulate.
• Error rate (float) - The sequencing error rate (between 0 and 1).

Example: 
chmod a+x simulate.sh
 ./simulate.sh sample.fasta 30 50 0.01

 assemble.sh will take in the following parameters (in this order):
• Reads file (string) - A file of reads, as output by simulate.sh.
• k (int) - the size of k-mer to use when building a De Bruijn graph.

Example: 
chmod a+x assemble.sh
 ./assemble.sh reads.txt 49

 *** In order to run assemble, you must run simulator first. Since if not, then 
 the reads.txt file is not formed.

