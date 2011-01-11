#!/usr/bin/perl
# Example taken from Chapters 8 and 12 of "Beginning perl for bioinformatics"
# Modified by Ioannis Filippis

# Examples 12-1 and 12-2   Parse alignments from BLAST output file

use strict;
use warnings;
#use BeginPerlBioinfo;     # see Chapter 6 about this module
use Getopt::Long;

# declare and initialize variables
my $beginning_annotation = '';
my $ending_annotation = '';
my %alignments = (  );
my $alignment = '';
my @HSPs = (  );
my($expect, $query, $query_start, $query_end, $subject, $subject_start, $subject_end) = ('','','','','','','');
my %args = ();
my $residues;
my @residues;
my $origQueryResidues = '';
my $psipredResidues = '';
my $origQueryRes;
my $alignQueryRes;
my $psipred;
my $help;
my $queryFile = '';
my $queryFasta = '';
my $queryHeader = '';
my $queryId = '';
my $subjectGiven = '';
my $subjectFasta = '';
my $finalQueryFasta = '';
my $finalSubjectFasta = '';
my $queryOrig2Align;
my $queryAlign2Orig;
my $subjectOrig2Align;
my $subjectAlign2Orig;
my $psipredPred;
my $psipredConf;
my $dssp;

# prints usage 
# - if not all necessary arguments are passed or
# - help option is passed
GetOptions("help|?" =>\$help,
	   "input=s"=>\$args{input},
	   "pdbcode=s"=>\$args{pdbcode},
	   "chaincode=s"=>\$args{chaincode},
	   "alignment=s"=>\$args{alignment},
	   "residues:s"=>\$residues,
	   "secstruc:s"=>\$psipred);
if (defined $help) {
    usage();
    die;   
}
foreach (keys %args) {
    if (!(defined $args{$_} && exists $args{$_})) {
	usage();
	die;
    }
}

$queryFile = join( '', get_file_data($args{input}));
($queryHeader, $queryFasta) = ($queryFile =~ />(.*?)\n(.*)/ms);
$queryFasta =~ s/\n//g;
($queryId) = ($queryHeader =~ /^([^ \t]+)/);

$subjectGiven = $args{pdbcode}."|".$args{chaincode};
$subjectGiven =~ tr/a-z/A-Z/;

$subjectFasta = `dumpseq -p $args{pdbcode}$args{chaincode} -s -N`;
$subjectFasta =~ s/\n//;

if (defined $psipred) {
    ($psipredPred, $psipredConf) = parse_psipred($psipred);
    $psipredPred =~ tr/C/ /;
    $psipredConf =~ tr/C/ /;

    #This works only in shell : DO NOT DELETE
    #$dssp= `dumppdb -p $args{pdbcode}$args{chaincode} -s | dssp -- 2>/dev/null | perl -pe 's/^.*\.$//' | sed -e '/^$/d' | sed -e '1d' | cut -c 17 | perl -pe 's/\n//'`; 
    `dumppdb -p $args{pdbcode}$args{chaincode} -s | dssp -- 2>/dev/null 1>$args{pdbcode}$args{chaincode}.dssp`;
    $dssp = parse_dssp("$args{pdbcode}$args{chaincode}.dssp");
    unlink("$args{pdbcode}$args{chaincode}.dssp");
    $dssp =~ tr/[HGI]/H/;
    $dssp =~ tr/[EB]/E/;
    $dssp =~ tr/[TS ]/ /;
}

parse_blast(\$beginning_annotation, \$ending_annotation, \%alignments, $args{alignment});
if ($alignments{$subjectGiven}) {
    $alignment = $alignments{$subjectGiven};
    @HSPs = parse_blast_alignment_HSP($alignment);
    ($expect, $query, $query_start, $query_end, $subject, $subject_start, $subject_end) = extract_HSP_information($HSPs[1]);

    $finalQueryFasta = "-" x ($subject_start - 1) . substr($queryFasta, 0, $query_start-1).$query.substr($queryFasta, $query_end) . "-" x (length($subjectFasta) - $subject_end);
    $finalSubjectFasta = substr($subjectFasta, 0, $subject_start-1) . "-" x ($query_start - 1) . $subject. "-" x (length($queryFasta) - $query_end) . substr($subjectFasta, $subject_end);
    
    ($queryOrig2Align, $queryAlign2Orig) = get_mapping($finalQueryFasta);
    ($subjectOrig2Align, $subjectAlign2Orig) = get_mapping($finalSubjectFasta);
    
    if(defined $psipred) {
    	$psipredPred = get_align_seq($psipredPred, $queryAlign2Orig, '-');
    	$psipredConf = get_align_seq($psipredConf, $queryAlign2Orig, '-');
    	$dssp = get_align_seq($dssp, $subjectAlign2Orig, '-');    	
    }
    
    if(defined $residues) {
	@residues = split(/,/, $residues);
	for my $res (@residues) {
	    $finalSubjectFasta = lowercase($finalSubjectFasta, $subjectOrig2Align, $res);
	    $alignQueryRes = $subjectOrig2Align->[$res-1];
	    $origQueryRes = $queryAlign2Orig->[$alignQueryRes-1];
	    $origQueryResidues .= $origQueryRes.",";
	    $finalQueryFasta = lowercase($finalQueryFasta, $queryOrig2Align, $origQueryRes);
	    if(defined $psipred) {
	    	$dssp = lowercase($dssp, $subjectOrig2Align, $res);
	    	$psipredPred = lowercase($psipredPred, $queryOrig2Align, $origQueryRes);
	    	if (substr($psipredPred, $alignQueryRes-1, 1) eq substr($dssp, $alignQueryRes-1, 1)) {
	    		$psipredResidues .= $origQueryRes.",";	
	    	} else {
	    		$psipredResidues .= "-1,";
	    	}
	    }
	}
	$origQueryResidues =~ s/,$//;
	$psipredResidues =~ s/,$//;
    }

    if (length($finalQueryFasta) != length($finalSubjectFasta)) {
	die "Calculated sequences differ in length!!!!!!!!!!!\n";
    } else {
	print ">".$queryHeader." | $query_start-$query_end | ".($subject_start-1+$query_start)."-".($subject_start-1+$query_start-1+length($query))."\n".$finalQueryFasta."\n";
	print ">".$args{pdbcode}.$args{chaincode}." | $subject_start-$subject_end | ".($query_start-1+$subject_start)."-".($query_start-1+$subject_start-1+length($subject))."\n".$finalSubjectFasta."\n";
    	if(defined $residues) {
    		print "\n\nOriginal Subject Residues : $residues\n";
    		print "Original Query Residues   : $origQueryResidues\n";
    	}
    	if(defined $psipred) {	
    		print "Psipred Agreement Residues: $psipredResidues\n";
    		print "\n\n$queryId: $finalQueryFasta\n$queryId: $psipredConf\n$queryId: $psipredPred\n";
    		print "$args{pdbcode}$args{chaincode}: $finalSubjectFasta\n$args{pdbcode}$args{chaincode}: $dssp\n";
    	}
    }
} else {
    die "No such pdb found in the blast file!!!!!!!!!!!!\n";
}

exit;

sub usage {
    print "Usage:\n perl parse_blast.pl\n";
    print "\t  -i  :  file with input target sequence in FASTA format\n";
    print "\t  -p  :  pdb code of the subject/template sequence\n";
    print "\t  -c  :  pdb chain code of the subject/template sequence\n";
    print "\t  -a  :  psi-blact classic output file for input target sequence\n";
    print "\t [-r] :  comma-separated list of residues to be examined. Optional\n";
    print "\t [-s] :  psipred output file. Optional\n";
    exit;
}

sub get_mapping {
	my $seq = shift;
	my @chars = split(//,$seq);
	my $origCharPos = 0;
	my $alignCharPos = 0;
       	my @orig2align = ();
       	my @align2orig = ();
	
	for my $char (@chars) {
	    $alignCharPos++;
	    if (!($char eq "-")) {
		$origCharPos++;
		push(@orig2align, $alignCharPos);
		push(@align2orig, $origCharPos);
	    } else {
	    	push(@align2orig, -1);
	    }
	}
	
	return (\@orig2align, \@align2orig);
}

sub lowercase {
	my ($seq, $orig2align, $origCharPos) = @_;
	my $alignCharPos = $orig2align->[$origCharPos-1];
	return (substr($seq, 0, $alignCharPos-1).lcfirst(substr($seq, $alignCharPos-1)));
}

sub get_align_seq {
	my ($seq, $align2orig, $gapChar) = @_;
	my @map = @{$align2orig};
	my $alignSeq = '';
	
	for my $pos (@map) {
		if ($pos == -1) {
			$alignSeq .= $gapChar;		
		} else {
			$alignSeq .= substr($seq, ($pos-1), 1);
		}
	}
	
	return $alignSeq;
}

sub parse_psipred {

    my $filename = shift;
    my $psipred_output_file = '';
    my $pred = '';
    my $conf = '';
    
    $psipred_output_file = join( '', get_file_data($filename));
    $pred = join ( '' , ($psipred_output_file =~ /^Pred: (.*)\n/gm) );
    $conf = join ( '' , ($psipred_output_file =~ /^Conf: (.*)\n/gm) );
    
    return ($pred, $conf);
}

sub parse_dssp {

    my $filename = shift;
    my $dssp_output_file = '';
    my $dssp = '';
    my $mainPart = '';

    $dssp_output_file = join( '', get_file_data($filename));
    ($mainPart) = ($dssp_output_file =~ /.*CA\s+(^.*)/ms);
    $dssp = join ( '' , ($mainPart =~ /^.{16}(.{1}).*\n/gm) );
    $dssp =~ s/\n+$//;

    return ($dssp);
}

# get_file_data
#
# A subroutine to get data from a file given its filename

################################################################################
# Subroutine from Chapter 8
################################################################################

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

################################################################################
# Subroutines for Example 12-1
################################################################################

# parse_blast
#
# -parse beginning and ending annotation, and alignments,
#     from BLAST output file

sub parse_blast {

    my($beginning_annotation, $ending_annotation, $alignments, $filename) = @_;

    # $beginning_annotation-reference to scalar
    # $ending_annotation   -reference to scalar
    # $alignments          -reference to hash
    # $filename            -scalar
    
    # declare and initialize variables
    my $blast_output_file = '';
    my $alignment_section = '';
    
    # Get the BLAST program output into an array from a file
    $blast_output_file = join( '', get_file_data($filename));

    # Extract the beginning annotation, alignments, and ending annotation
    ($$beginning_annotation, $alignment_section, $$ending_annotation)
    = ($blast_output_file =~ /(.*?)(^>.*)(^  Database:.*)/ms);

    # Populate %alignments hash
    # key = ID of hit
    # value = alignment section
    %$alignments = parse_blast_alignment($alignment_section);
}

# parse_blast_alignment
#
# -parse the alignments from a BLAST output file,
#       return hash with
#       key = ID
#       value = text of alignment

sub parse_blast_alignment {

    my($alignment_section) = @_;
    
    # declare and initialize variables
    my(%alignment_hash) = (  );

    # loop through the scalar containing the BLAST alignments,
    # extracting the ID and the alignment and storing in a hash
    #
    # The regular expression matches a line beginning with >,
    # and containing the ID between the first pair of | characters;
    # followed by any number of lines that don't begin with >

    while($alignment_section =~ /^>.*\n(^(?!>).*\n)+/gm) {
        my($value) = $&;
	my($key) = ($value =~ /pdb\|(.*?) /);
        $alignment_hash{$key} = $value;
    }

    return %alignment_hash;
}

################################################################################
# Subroutines for Example 12-2
################################################################################

# parse_blast_alignment_HSP
#
# -parse beginning annotation, and HSPs,
#     from BLAST alignment
#     Return an array with first element set to the beginning annotation,
#    and each successive element set to an HSP

sub parse_blast_alignment_HSP {

    my($alignment ) = @_;

    # declare and initialize variables
    my $beginning_annotation  = '';
    my $HSP_section  = '';
    my @HSPs = (  );
    
    # Extract the beginning annotation and HSPs
    ($beginning_annotation, $HSP_section )
        = ($alignment =~ /(.*?)(^ Score =.*)/ms);

    # Store the $beginning_annotation as the first entry in @HSPs
    push(@HSPs, $beginning_annotation);

    # Parse the HSPs, store each HSP as an element in @HSPs
    while($HSP_section =~ /(^ Score =.*\n)(^(?! Score =).*\n)+/gm) {
        push(@HSPs, $&);
    }

    # Return an array with first element = the beginning annotation,
    # and each successive element = an HSP
    return(@HSPs);
}

# extract_HSP_information
#
# -parse a HSP from a BLAST output alignment section
#        - return array with elements:
#    Expect value
#    Query string
#    Query range 
#    Subject string
#    Subject range

sub extract_HSP_information {

    my($HSP) = @_;
    
    # declare and initialize variables
    my($expect) = '';
    my($query) = '';
    my($query_start) = '';
    my($query_end) = '';
    my($subject) = '';
    my($subject_start) = '';
    my($subject_end) = '';

    ($expect) = ($HSP =~ /Expect = (\S+),/);

    $query = join ( '' , ($HSP =~ /^Query.*\n/gm) );

    $subject = join ( '' , ($HSP =~ /^Sbjct.*\n/gm) );

    ($query_start, $query_end) = ($query =~ /(\d+).*\D(\d+)/s);

    ($subject_start, $subject_end) = ($subject =~ /(\d+).*\D(\d+)/s);

    #$query =~ s/[^acgt]//g;

    #$subject =~ s/[^acgt]//g;

    $query = join ( '' , ($HSP =~ /^Query:[\s\d]+([^\s]+)[\s\d]+\n/gm) );

    $subject = join ( '' , ($HSP =~ /^Sbjct:[\s\d]+([^\s]+)[\s\d]+\n/gm) );
    
    return ($expect, $query, $query_start, $query_end, $subject, $subject_start, $subject_end);
}
