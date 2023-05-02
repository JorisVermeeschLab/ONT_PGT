#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use File::Path;
use Cwd 'abs_path';

######################################################################
#
# Variables initialization with default values
#
######################################################################

my ($record,@words,@unlinkfiles,@SamplesVect);

my ($Program_Folder_Path,$Main_Output_Folder_Path,$Sample_Out_Folder_Path,$Sample_Label,$Analysis_Label,$Assembly,$Target_Name,$Target_File_Path,$RC_Folder_Path,$RCNorm_Folder_Path,$Images_Folder_Path,$R_Target_Folder,$R_Target_Path,$R_Norm_Path,$Input_File_Path,$Mode,$Sample_Res_Folder_Path,$Sample_Plot_Folder_Path,$Sample_Data_Folder_Path,$processors);

my ($verbose,$help,$man);

$Mode = "antani";
$Target_Name = "supercazzola";


######################################################################
#
#  Defining system variables
#
######################################################################

my ($myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);


$Program_Folder_Path="$workingfolder";


######################################################################
#
# Reading user's options
#
######################################################################

GetOptions('processors=s'=>\$processors,'assembly|a=s'=>\$Assembly,'target|t=s'=>\$Target_Name,'output|o=s'=>\$Main_Output_Folder_Path,'mode|e:s'=>\$Mode,'verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error: the number of arguments found at command line is incorrect.");

$Input_File_Path=$ARGV[0];




die ("User ERROR: Please set option --mode with a valid string ('pooling' or 'paired' or 'nocontrol').\n") if( $Mode ne "pooling" && $Mode ne "paired" && $Mode ne "nocontrol" );


#die ("User ERROR: the target specified with --target wasn't found for assembly $Assembly. \nPlease check --target and --assembly. \nError found") if (! -d "$Program_Folder_Path/data/targets/$Assembly/$Target_Name" );

######################################################################
#
#  Checking system folders
#
######################################################################

print "Checking output folders...\n";

if(-e $Main_Output_Folder_Path){ 
  print "'$Main_Output_Folder_Path' folder ready!\n";
}
else{
  mkpath($Main_Output_Folder_Path);
  if(-e $Main_Output_Folder_Path){
    print "'$Main_Output_Folder_Path' folder created!\n";
  }  
}

print "Checking output subfolders...\n";


if(-e "$Main_Output_Folder_Path/Results"){
  print "'$Main_Output_Folder_Path/Results' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/Results");
  print "Creating '$Main_Output_Folder_Path/Results' folder...\n";
}

if(-e "$Main_Output_Folder_Path/Plots"){
  print "'$Main_Output_Folder_Path/Plots' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/Plots");
  print "Creating '$Main_Output_Folder_Path/Plots' folder...\n";
}

if(-e "$Main_Output_Folder_Path/.tmp"){
  print "'$Main_Output_Folder_Path/.tmp' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/.tmp");
  print "Creating '$Main_Output_Folder_Path/.tmp' folder...\n";
}


######################################################################
#
#  When mode is pooling generates the Read Count for Controls
#
######################################################################


if($Mode eq "pooling"){ 
  print "Creating Pooling Control!\n";
  system qq(R --slave --args $Program_Folder_Path,$Main_Output_Folder_Path,$Program_Folder_Path/data/targets/$Assembly/$Target_Name,$Input_File_Path,$Mode,$Target_Name < $Program_Folder_Path/lib/R/PoolingCreateControl.R);
  print "Pooling Control Created!\n";
}




######################################################################
#
#  Creating Multi Processor Analysis
#
######################################################################


print "Preparing Multiprocessor Analysis...\n";
system qq(R --slave --args $Program_Folder_Path,$Input_File_Path,$Main_Output_Folder_Path,$Target_Name,$Assembly,$processors,$Mode < $Program_Folder_Path/lib/R/DataAnalysisParallel.R);
print "Multiprocessor Analysis Complete\n";
print "Starting Multiprocessor Analysis!\n";


######################################################################
#
#  Running Multi Processor Analysis
#
######################################################################

my $Input_File_Parallel="$Main_Output_Folder_Path/.tmp/ParallelRCAnalysis.sh";


open(CHECKBOOK,"$Input_File_Parallel") || die "Couldn't open the input file $Input_File_Parallel.";
my @pids;
while($record=<CHECKBOOK>){
  my $childpid = fork() or exec($record);
  push(@pids,$childpid);
}

print "My Children: ", join(' ',@pids), "\n";
waitpid($_,0) for @pids;

print "Multiprocessor Analysis Complete!\n";




######################################################################
#
#  Documentation
#
######################################################################

=head1 SYNOPSIS 

 =head1 SYNOPSIS 

 NanoGLADIATORAnalysis.pl [arguments] [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.
           --output                 Select mapping quality for .bam file filtering; if omitted default value is 0.
           --processors             Select the number of processors for parallel sample analysis.
           --assembly               Select assembly name. Available options: "hg19" or "hg38".
           --target                 Select Target window already initialized.
           --mode                   Select mode for experimental design between paired and nocontrol.

 Function:
 
 NanoGLADIATORAnalysis.pl performs segmentation of the normalized RC and predict the state (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification) and the allelic fraction of each segmented region.

 Example: 
 
 NanoGLADIATOR> perl NanoGLADIATORAnalysis.pl ExperimentalFileAnalysis.txt --processors 6 --target w1000 --assembly hg19 --output /.../OutXCAVATOR/Results_MyProject_w50K --mode paired/nocontrol

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.

=item B<--processors>

The number of thread to use for the analysis.

=item B<--output>

The output folder for resulting files.

=item B<--assembly>

The assembly exploited for read mapping and target initialization. 

=item B<--target>

The "target name" used for target initialization with ReferenceWindowInitialize.pl.

=item B<--mode>

The experimental design mode to use. The possible options are "pooling", "paired" or "nocontrol".

=back

=head1 DESCRIPTION

NanoGLADIATORAnalysis.pl performs the segmentation of the RC by means of the Shifting Level Model algorithm and exploits FractionCall algorithm to classify each segmented region as one of the five possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification) and predict its allelic fraction.


Nano-GLADIATOR is freely available to the community for non-commercial use. For questions or comments, please contact "albertomagi@gmail.com".
=cut


