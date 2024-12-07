# VRMOD
Vertebrate Regulatory MOdule Detector (VRMOD), an algorithm to accurately predict cis-regulatory modules in the vertebrate genomes.

To set VRMOD up:
1. Download the following code:
  VRMOD_V0.1.pl, 
  alphabet, 
  patser-v3e

2. Download the position weight matrix files: Mouse_phylonet_L30_Short_combined_Strategy2_final.tgz and unzip it.
   
3. modify the script VRMOD_V0.1.pl at the following lines to use your system-specific path to the files and director.
   
software path
my $path_patser = "/ref/gzlab/software/VRMOD_V0.1/patser-v3e ";                                                      
my $path_alphabet = "/ref/gzlab/software/VRMOD_V0.1/alphabet";

 input position weigth matrix                                                                                       
my $dir = "/ref/gzlab/data/Mouse_phylonet/Mouse_phylonet_L30_Short_combined_Strategy2_final";  

4. example command line
   perl  VRMOD_V0.1.pl "input dir path" "fasta input sequence file name"
   
5. output file will be automatically generated in the input directory named .pred_mod.txt . Each input file will have a corresponding output file. If the command was successfully completed a line of "# program finished ." will be printed at the end of the output file.
   
