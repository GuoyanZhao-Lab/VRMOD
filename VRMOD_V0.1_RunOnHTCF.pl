#!/usr/bin/perl
use strict;

my $usage = '
This script will run VRMOD_V0.1.pl on htcf.

perl VRMOD_V0.1_RunOnHTCF.pl <input dir> <sequence file>

<input dir> = directory of input .fa file
<sequence file> = fasta format sequence file name in the input dir

This script will call VRMOD_V0.1.pl

===============================================================================

';
die $usage unless scalar @ARGV == 2;
my ( $input_dir, $input_file ) = @ARGV;

my $HOME = $ENV{HOME};

# may need to modify to adapt the your operating system
#my $run_script_path = `dirname $0`;
#chomp $run_script_path;
#my $run_script_path = "/usr/bin/perl /ref/gzlab/software/PhyloNet_tools_Mammal/";
my $run_script_path = "/usr/bin/perl /ref/gzlab/software/VRMOD_V0.1/";

# no need to change
my $job_file = $input_dir."/PredictCRM_job".$$.".sh";
#my $job_file = "PredictCRM_job".$$.".sh";
open(STCH, ">$job_file") or die "can not open file $job_file!\n";
chmod 0755, $job_file;
print STCH "#!/bin/bash\n";
#print STCH "#SBATCH --mem=240G\n";
print STCH "#SBATCH --mem=100G\n";
print STCH "#SBATCH --error=".$job_file.".error \n";
print STCH "#SBATCH --output=".$job_file.".out \n";
print STCH "INDIR=$input_dir \n";
print STCH "INFILE=$input_file \n";
print STCH "TEMP=/scratch/gzlab/temp/temp_gzhao/ \n";
#print STCH "set -x\n";

#print STCH $run_script_path."/"."ModulePred_S1_NoCompositionMotifs_Mammal_ForOneSeq.pl \${INDIR}  \${INFILE}    \n\n";
print STCH $run_script_path."/"."VRMOD_V0.1.pl \${INDIR}  \${INFILE}    \n\n";
close STCH;

`sbatch $job_file `;

exit;


