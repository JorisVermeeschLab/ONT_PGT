#!/usr/bin/perl

use warnings;
use Getopt::Long;
use strict;
use File::Path;
use Cwd 'abs_path';
use File::Find;
#use Forks::Super;

my ($target,$assembly,$MAPQ,$Program_Folder_Path);
my ($record,@words,@unlinkfiles,@SamplesVect);
my ($verbose,$help,$man,$mode);
my ($PathIn,$PathOut,$NReads,$LabelOut,$GenomeRef,$processors,$fout1);

$MAPQ=20;
######################################################################
#
#  Reading user's options
#
######################################################################

GetOptions( 'label=s' => \$LabelOut,
            'ref=s' => \$GenomeRef,
            'nreads=s' => \$NReads,
            'pathout=s' => \$PathOut,
            'pathin=s' => \$PathIn,
            'processors=s' => \$processors,
            'target=s' => \$target,
            'assembly=s' => \$assembly,
            'mode=s' => \$mode,
            'mapq=i'=>\$MAPQ
) or die "Invalid arguments!";

#die "Missing --label " unless $VCFIn;
#die "Missing -O " unless $OutFolder;
#die "Missing -L " unless $OutName;



######################################################################
#
# Defining system variables
#
######################################################################

my ($myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);


$Program_Folder_Path="$workingfolder";

#########################################
#                                       #
#  Checking Fast5 folder for Files      #
#                                       #
#########################################

sub count_file{
my $source = shift;
my $files = 0;
File::Find::find( sub {
    if (-f $File::Find::name) {
        $files++
    } 
}, $source );
return $files;
}

my $NumFile = 0;

while ($NumFile < $NReads) {
$NumFile = count_file($PathIn);
}

print "$NumFile\n";

mkpath ("$PathOut");

#########################################
#                                       #
#  Creating FastQ from Fast5 Files      #
#                                       #
#########################################

my $FastqOut = "$PathOut/$LabelOut.fastq";


print "Creating Fastq File with Poretools...\n";
system qq(poretools fastq $PathIn > $FastqOut);
print "Fastq file created!!\n";

#########################################
#                                       #
#  Aligning FastQ with Minimap2         #
#                                       #
#########################################

my $SamOut = "$PathOut/$LabelOut.sam";

print "Starting minimap2 alignment...\n";
system qq(minimap2 --MD -ax map-ont -t $processors $GenomeRef $FastqOut > $SamOut);
print "Sam file created!!\n";


my $BamOut = "$PathOut/$LabelOut.bam";

print "Converting sam to bam ...\n";
system qq(samtools view -b $SamOut > $BamOut);
print "Bam file created!!\n";


my $SortBamOut = "$PathOut/$LabelOut.srt.bam";

print "Sorting bam file ...\n";
system qq(samtools sort -@ $processors $BamOut > $SortBamOut);
print "Bam file sorted!!\n";

print "Indexing sorted bam file ...\n";
system qq(samtools index $SortBamOut);
print "Sorted Bam file indexed!!\n";


print "Removing support files ...\n";
system qq(rm $SamOut $BamOut $FastqOut);


###################################################
#                                                 #
#  Create FilePrepare and File Analysis......     #
#                                                 #
###################################################
my $FilePrepare = "$PathOut/$LabelOut.prepare.txt";
my $FileAnalysis = "$PathOut/$LabelOut.analysis.txt";
my $PathPrepare = "$PathOut/Prepare/$LabelOut";
my $PathAnalysis = "$PathOut/Analysis";

mkpath ("$PathOut/Prepare");
mkpath ("$PathAnalysis");


open $fout1, '>', $FilePrepare or die $!;
print $fout1 $SortBamOut;
print $fout1 " ";
print $fout1 $PathPrepare;
print $fout1 " ";
print $fout1 "$LabelOut\n";
close $fout1;

open $fout1, '>', $FileAnalysis or die $!;
print $fout1 "T1";
print $fout1 " ";
print $fout1 $PathPrepare;
print $fout1 " ";
print $fout1 "$LabelOut\n";
close $fout1;


#########################################
#                                       #
#  Running Data Prepare........         #
#                                       #
#########################################


system qq(perl $Program_Folder_Path/lib/perl/ReadPerla.pl $FilePrepare --assembly $assembly --target $target);

#########################################
#                                       #
#  Running Data Analysis.......         #
#                                       #
#########################################


system qq(perl $Program_Folder_Path/lib/perl/AnalyzePerla.pl $FileAnalysis --assembly $assembly --output $PathAnalysis --target $target --mode nocontrol)



