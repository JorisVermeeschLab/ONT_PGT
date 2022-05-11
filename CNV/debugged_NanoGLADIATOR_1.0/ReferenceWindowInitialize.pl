#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use Cwd 'abs_path';
use File::Path;
  
######################################################################
#
# Variables initialization with default values
#
######################################################################

our ($record,@words,@unlinkfiles);

our ($Assembly,$Window,$Target_File_Out,$mode);

our ($Program_Folder_Path,$Target_Name,$Source_Data);

our ($Path_2_Wig,$Path_2_fasta);

our ($verbose,$help,$man);

our $Target_Filt_Path;

######################################################################
#
# Defining options
#
######################################################################



GetOptions( 'M=s' => \$mode,
            'S=s' => \$Source_Data,
            'W=s' => \$Window,
            'A=s' => \$Assembly,
            'T=s' => \$Target_Name,
) or die "Invalid arguments!";

die "Missing -M " unless $mode;
die "Missing -S " unless $Source_Data;
die "Missing -W " unless $Window;
die "Missing -A " unless $Assembly;
die "Missing -T " unless $Target_Name;

######################################################################
#
# Defining system variables
#
######################################################################



my ($myscriptname,$myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);

$Program_Folder_Path="$workingfolder";

print "Program folder is: $Program_Folder_Path\n";


$Target_Filt_Path="$Program_Folder_Path/data/targets/$Assembly/$Target_Name";




######################################################################
#
# Removing empty lines from source file
# and creating support temporary files
#
######################################################################

my $range=99999;
my $ID=int(rand($range));
my $filename="source.$ID";

if(-e $filename){ 
  $ID=$ID+100000;
}

system qq(awk NF $Source_Data > source.$ID);

######################################################################
#
# Reading source file
#
######################################################################

open(CHECKBOOK,"source.$ID") || die "Couldn't open the source file!";

while($record=<CHECKBOOK>){
  chomp($record);
  @words=split(' ',$record);
  
  $Path_2_Wig=$words[0];
  $Path_2_fasta=$words[1];
      

######################################################################
#
# Target initialization
#
######################################################################

if ($mode eq 'genome'){
print "Creating Consecutive and non-overlapping windows of $Window bp Size \n";
system qq(R --slave --args $Program_Folder_Path,$Target_Name,$Assembly,$Window,$Path_2_fasta < $Program_Folder_Path/lib/R/CreateBedGenome.R);
}


$Target_File_Out="$Program_Folder_Path/data/targets/$Assembly/$Target_Name/TargetWindow.bed";
print "Calculating Mappability and GC content...\n";
system qq($Program_Folder_Path/lib/bash/./TargetCreate.sh $Path_2_Wig $Target_File_Out $Program_Folder_Path $Assembly $Target_Name $Path_2_fasta);
print "...done!\n";

}


close(CHECKBOOK);
@unlinkfiles=("source.$ID");
unlink @unlinkfiles;

######################################################################
#
# Documentation
#
######################################################################

=head1 SYNOPSIS 

 perl ReferenceWindowInitialize.pl [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.
       -M                           Set mode (at present only the genome option is allowed) to inizitialize windows information.
       -S                           The source file that contain the complete path to the mappability and reference file.
       -T                           The name of the windows information
       -A                           The assembly name (hg19 or hg38)
       -W                           The window size in bp (10000)
 Function:
 
ReferenceWindowInitialize.pl initialises windows informations for further data processing with the Nano-GLADIATOR package. 
It requires 4 arguments (one source files - with space-delimited paths to source data for mappability and GC-content calculations),
 a label name, window size and assembly to run properly. 
 A sub-folder with the specified target name will be created under "NanoGLADIATOR/data/targets/hgXX".

 Example: 
 
 NanoGLADIATOR> perl ReferenceWindowInitialize.pl -M genome -S SourceTarget.txt -T LabelName -W 50000 -A hg19
 
=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.


=back

=head1 DESCRIPTION

ReferenceWindowInitialize.pl is a Perl script which is part of the NanoGLADIATOR package. It includes all of the first step operations of the NanoGLADIATOR package. It calculates GC-content and Mappability for consecutive and non-overlapping windows of size W.

It requires, as arguments, the path to a source file (the default source file is "SourceTarget.txt" which is placed in the main NanoGLADIATOR folder) containing the paths to source data (for the calculations of mappability and GC-content), the path to the target input file, a "target name", the window size and the assembly. Setting the label name as "LabelName", all data calculated will be saved in the "LabelName" folder in (if you are using the hg19 assembly) NanoGLADIATOR/data/targets/hg19/LabelName.

The allowed assemblies are hg19 and hg38.

NanoGLADIATOR is freely available to the community for non-commercial use. For questions or comments, please contact "albertomagi@gmail.com".

=cut

