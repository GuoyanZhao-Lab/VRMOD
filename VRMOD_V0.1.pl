#!/usr/bin/perl -w
use strict;

my $usage = ' This script is designed to find out how the 5000
orso motifs are distributed in each given promoter sequence and predicts
modules based on the distribution.
This is the final version of VRMOD without the need to change any parameters.
It can be used for both single sequence and multiple sequence input prediction. 

perl script <input dir> <sequence file name>
<input dir> = path to the dir
<sequence file> = fasta format sequence file name

output prediction to automatically generated file named .pred_mod.txt

';

die $usage unless @ARGV == 2;
my ($input_dir, $input_file) = @ARGV;

my $temp_dir = "/scratch/gzlab/temp/";
###########################################################################
my $HOME = $ENV{HOME};
my $seq_f = $input_dir."/".$input_file;
my $out_file = $input_file;
$out_file =~ s/\.fa|\.fasta/.pred_mod.txt/;
$out_file = $input_dir."/".$out_file;

# software path
my $path_patser = "/ref/gzlab/software/VRMOD_V0.1/patser-v3e ";
my $path_alphabet = "/ref/gzlab/software/VRMOD_V0.1/alphabet";

# input position weigth matrix
my $dir = "/ref/gzlab/data/Mouse_phylonet/Mouse_phylonet_L30_Short_combined_Strategy2_final";

############################################################################
# significance cutoff 
my $Z_score_cutoff = 2.3263; # corresponding to p vlaue = 0.01 of one-tailed
#my $Z_score_cutoff = 3.09; # corresponding to p vlaue = 1.0E-3
#my $Z_score_cutoff = 3.719; # corresponding to p vlaue = 1.0E-4
#my $Z_score_cutoff = 4.265; # corresponding to p vlaue = 1.0E-5
my $distance_cutoff = 30; # distance between the nearby peakes
my $distance_cutoff_CRM = 30; # distance between the nearby CRMs to merge

##########################################################################
# Read matrices and descriptions in dir
my @names = (); # Names of the matrix file and desc file
my %width = (); # matrix name => width of matrix
&read_names( \@names, $dir ); # Read list of matrix names in dir
&get_width($dir, \@names, \%width);

##################################################################
# read input query sequence. 
my %seq = (); # gene name => sequence
&getGeneSeq($seq_f, \%seq);
#print "finished geting gene seq\n";

##################################################################
# get all motif sites for all the sequences
my %sites = (); #  [all the sites for all matrices in all genes]}
&get_allMotifSites_for_multiple_seq(\%sites, \@names, \%seq);

#################################################################
# for each gene, find modules in it
# open output file
open(OUT, ">$out_file" ) or die "Can not open $out_file!\n";	

foreach my $gene (keys %sites) {
	my $seq_length = length($seq{$gene});

	my $modules_ref = &find_module($gene, $seq_length, \@{$sites{$gene}});
	foreach my $mod (@{$modules_ref}) {
		print OUT  $mod, "\n";
	}
}	
print OUT  "# program finished .\n";
close(OUT);

exit;

################################################################
sub find_module() {
	my ($gene, $length, $sites_ref) = @_;

	# bin the positions, each position is a bin  
	my $binNum = $length;
	my @bins = ();
	for (my $i = 0; $i < $binNum; $i++) {
		$bins[$i] = 0;
	}

	foreach my $pos (@{$sites_ref}) {
		$bins[$pos]++;
	}

	# calculate average, SD
	my $sum = 0;
	for (my $i = 0; $i < $binNum; $i++) {
		$sum += $bins[$i];
	}
	my $ave = $sum/$binNum;
	my $variance = 0;
	for (my $i = 0; $i < $binNum; $i++) {
		$variance += ($bins[$i] - $ave)*($bins[$i] -  $ave);
	}
	$variance = $variance/$binNum;
	my $SD = sqrt($variance);
#	print "ave = $ave, SD = $SD\n";
	if ($SD == 0) {
		return;
	}

	my @Z_score = ();
	for (my $i = 0; $i < $binNum; $i++) {
		$Z_score[$i] = ($bins[$i] - $ave)/$SD;
		#print "value = $bins[$i], z score = $Z_score[$i]\n";
	
#		print out Z score distribution
#		print $i+1, "\t", $Z_score[$i], "\n";
	}

	# first find all the peaks that have z score >= cutoff
	my @peaks = ();
	for (my $i = 0; $i < (scalar @Z_score); $i++) {
		if ($Z_score[$i] >= $Z_score_cutoff) {
			push @peaks, $i;
		}
	}
	my %peak_use = (); # position => used 1, not used 0
	foreach my $peak_pos (@peaks) {
		$peak_use{$peak_pos} = 0;
	}

	# extend from each peak
	my $startPos = 0;
	my $lastPos = 0;
	my $distance = 0; # distance from the last peak being included to the next negative z score position 
	my @module = ();
	foreach my $peak_pos (@peaks) {
		my $startPos = $peak_pos;
		my $lastPos = $peak_pos;
		if (!$peak_use{$peak_pos}) {
			# extend backward 5'
			my $distance_start = $peak_pos; # start for calculate distance
			for (my $i = $peak_pos - 1; $i > 0; $i--) {
				if ($Z_score[$i] > 0) { # position has positive Z score
					$startPos = $i + 1; # sequence start position
					if ($Z_score[$i] >= $Z_score_cutoff) { # is a peak
						$peak_use{$i} = 1;
						$distance_start = $i;
					}
				}
				else { # position has negative Z score
					$distance = $distance_start - $i + 1;
					if ($distance > $distance_cutoff) {
						last;
					}
				}
			}

			# extend forward 3'
			for (my $i = $peak_pos + 1; $i < $binNum; $i++) {
				if ($Z_score[$i] > 0) { # position has positive Z score
					$lastPos = $i + 1; # sequence end position
					if ($Z_score[$i] >= $Z_score_cutoff) { # is a peak
						$peak_use{$i} = 1;
						$distance_start = $i;
					}
				}
				else { # position has negative positive Z score
					$distance = $i - $distance_start + 1;
					if ($distance > $distance_cutoff) {
						last;
					}
				}
			}
	
#			my $mod_start = $startPos - 15; # subtract 15 bp to cover the possible longest motif
			my $mod_start = $startPos - 30; # subtract 30 bp to cover the possible longest motif
			if ($mod_start < 0) { $mod_start = 0;}
#			my $mod_end = $lastPos + 15;  # add 15 bp to cover the possible longest motif
			my $mod_end = $lastPos + 30;  # add 30 bp to cover the possible longest motif
			if ($mod_end > $length) { $mod_end = $length;}

			my $mod = $gene."_".$mod_start."_".$mod_end;

#			print "predicted module: $mod\n";
			push @module, $mod;
		}
	}
	undef @bins;
	undef @Z_score;
=head1
	print "before procesing:\n";
	foreach my $mod (@module) {
		print $mod, "\n";
	}
=cut

	# merge modules if their ends are within $distance_cutoff_CRM bp
	my @post_processed_modules = ();
	if (scalar @module) {
#		print "at merging modules: $module[0] \n";

		my ($name11,  $start_i, $end_i) = split("_", $module[0]);
		my $new_start = $start_i;
		my $new_end = $end_i;
		my $i = 1;
		while( $i <= (scalar @module - 1)) {
			my ($name21,  $start_i, $end_i) = split ("_", $module[$i]);
			if ($start_i - $new_end <= $distance_cutoff_CRM) { 
				$new_end = $end_i;
				$i++;
			}
			else {
				my $p_modules = $gene."_".$new_start."_".$new_end;
				push @post_processed_modules, $p_modules;
	
				$new_start = $start_i;
				$new_end = $end_i;
				$i++;
			}
		}

		# for the last module
		my $p_module = $gene."_".$new_start."_".$new_end;
		push @post_processed_modules, $p_module;
	}

=head1	
	print "after procesing:\n";
	foreach my $mod (@post_processed_modules) {
		print $mod, "\n";
	}
	print "*****************\n\n";
=cut

	undef (@module);
	return \@post_processed_modules;
}


#####################################################################
# This function accept a hash that stores sequence information, generate 
# a file in consensus format, scan more all the motif sites using patser.
sub get_allMotifSites_for_multiple_seq() {
	my ($sites_ref, $names_ref, $seq_ref) = @_;

	my $pat_f = $temp_dir."/temp_pat".$$;
	my $con_f = $temp_dir."/temp_con".$$;
	open( Con,">$con_f") or die "can notp open $con_f.\n";
	foreach my $gene (keys %{$seq_ref}) {
		print Con "$gene \\ ";
		print Con $seq_ref->{$gene}, "\n\\\n";
	}
	close Con;

	foreach my $mat (@{$names_ref}) {
		# Run patser default cutoff and get sites for each motif
		my $command = $path_patser." -a ".$path_alphabet." -m ".$dir."/".$mat.".matrix -f $con_f -c -d2 -li > $pat_f ";
		#print $command, "\n";
		system( $command );
	
		open ( PAT, $pat_f ) || die "Cannot open $pat_f.";
		while ( <PAT> ) {
			if ( /\s+(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
#				print "name = $1, pos = $2, score = $3\n";
				my $name = $1;
				my $pos = $2;
				if ($pos =~ /C/) {
					$pos =~ s/C//;
				}
		
				# use the center of motif
				my $cen = $pos + int($width{$mat}*0.5);
				push @{$sites{$name}}, $cen;	
#				print "name = $name, pos = $cen\n";
			}
		}
		close ( PAT );
	}
	unlink $pat_f;
	unlink $con_f;
}

################################################################3
# This function read all the matrices and description files in dir
sub read_names()
{
    my ( $name_r, $dir ) = @_;
    opendir (DH, $dir) or die "Cannot open $dir";
    foreach my $file (readdir DH) {
		if ($file =~ /matrix/) {
		    my $name = $file;
	    	$name =~ s/\.matrix//;
#	   	 	print "matrix name = \\$name\\\n";
	    	push @$name_r, $name;
		}
    }
    close DH;
}


#############################################
# This function get the width of all the  matrices in dir
sub get_width {
    my ($dir, $name_ref, $width_ref) = @_;

    my @line;
    my $desc;

    for ( my $i=0; $i < scalar @{$name_ref}; $i++ )
    {
        my $alignMatrix = "";
        my ($A, $C, $G, $T) = 0;

        $desc = $dir."/".$name_ref->[$i].".desc";
        open ( IN, $desc ) || die "Cannot open $desc.";
        while ( <IN> ) {
            if ( /width\s*=\s*(\d+)/ )   {
                $width_ref->{$name_ref->[$i]} = $1;
                last;
            }
        }
        close ( IN );
    }
}


#############################################################
# Check whether an array has an element
sub array_has
{
  my ( $array_r, $ele ) = @_;
  for ( my $i=0; $i<=$#$array_r; $i++ )
  {
    if ( $ele eq $array_r->[$i] )
    {
      return 1;
    }
  }
  return 0;
}

#######################################################################
# read fasta file, For multiple gene input
sub getGeneSeq() {
    my ($fileName, $hashRef) =  @_;
    open (File, $fileName) or die "Can not open file $fileName !\n";

    # get all the gene Name and seq name=> seq
    my $oldSeperator = $/;
    $/ = ">"; # seperates two records
    my $line = <File>;
    while ($line = <File>){
		#print $line;
		my @fields = split ("\\n", $line);	
		my $First_line = shift @fields;

		my @temp = split (" ", $First_line);
		my $name = shift @temp;

		my $seq = join ("", @fields);
		$seq =~ s/>//g;
		$seq =~ s/\s//g;
		#lc $seq;
 
		print "name is \\$name\\\n";
		print "seq is \n";
		print "\\$seq\\\n";
		$hashRef->{$name} = $seq;
    }
    $/ = $oldSeperator;
	close File;
}

