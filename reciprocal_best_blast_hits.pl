#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Bio::SearchIO;
use File::Basename;

BEGIN {
	*STDOUT->autoflush();
	*STDERR->autoflush();
}

my $options               = check_params();
my %hits_blast_a          = ();
my @queries_names_blast_a = ();
my %hits_blast_b          = ();
my @queries_names_blast_b = ();

print "*** Parsing  a ...\n";
parse_blast_report( $options->{'a'}, \%hits_blast_a, \@queries_names_blast_a );
print "    Parsing  a ... DONE\n";
print "*** Parsing  b ...\n";
parse_blast_report( $options->{'b'}, \%hits_blast_b, \@queries_names_blast_b );
print "    Parsing  b ... DONE\n";
print "*** Finding  RBB a ...\n";
find_reciprocal_best_hits( \%hits_blast_a, \@queries_names_blast_a,
	\%hits_blast_b, 'a', $options->{'a'} );
print "    Finding  RBB a ...DONE\n";
print "*** Finding  RBB b ...\n";
find_reciprocal_best_hits( \%hits_blast_b, \@queries_names_blast_b,
	\%hits_blast_a, 'b', $options->{'b'} );
print "    Finding  RBB b ...DONE\n";

#--------------------------------------------------------------------------------
sub parse_blast_report {
	my ( $blast_report, $hits_blast_ref, $queries_names_blast_ref ) = @_;

	my $searchio = new Bio::SearchIO(
		-format => 'blast',
		-file   => $blast_report
	);

  RESULT:
	while ( my $result = $searchio->next_result ) {
		my $query_name = $result->query_name();
		push @$queries_names_blast_ref, $query_name;
		my $hits = $result->num_hits();

		if ( !$hits ) {
			$hits_blast_ref->{$query_name} = 'NONE';
			next RESULT;
		}

		my $hit      = $result->next_hit();
		my $hit_name = $hit->name();
		$hits_blast_ref->{$query_name} = $hit_name;
	}
}

sub find_reciprocal_best_hits {
	my ( $hits_a_ref, $names_a_ref, $hits_b_ref, $file_type, $file_name ) = @_;

	my $basename     = basename( $file_name, '.blastplus_report' );
	my $RBB_file     = 'RBB_' . $file_type . '_' . $basename . '.csv';
	my $NO_HIT_file  = 'NO_HIT_' . $file_type . '_' . $basename . '.txt';
	my $NOT_RBB_file = 'NOT_RBB_' . $file_type . '_' . $basename . '.csv';

	open my $NO_HIT_fh, '>', $NO_HIT_file
	  or croak("Can't open output file [$NO_HIT_file]: $!.\n");
	open my $RBB_fh, '>', $RBB_file
	  or croak("Can't open output file [$RBB_file]: $!.\n");
	open my $NOT_RBB_fh, '>', $NOT_RBB_file
	  or croak("Can't open output file [$NOT_RBB_file]: $!.\n");

  QUERY:
	foreach my $query_name_a (@$names_a_ref) {
		my $hit_name_a = $hits_a_ref->{$query_name_a};

		if ( $hit_name_a eq 'NONE' ) {
			print $NO_HIT_fh $query_name_a . "\n";
			next QUERY;
		}

		if ( $hits_b_ref->{$hit_name_a} eq $query_name_a ) {
			print $RBB_fh $query_name_a . ',' . "$hit_name_a" . "\n";
		}
		else {
			print $NOT_RBB_fh $query_name_a . ',' . "$hit_name_a" . "\n";
		}
	}

	close $NO_HIT_fh  or croak("Failed to close file [$NO_HIT_file]: $!\n");
	close $RBB_fh     or croak("Failed to close file [$RBB_file]: $!\n");
	close $NOT_RBB_fh or croak("Failed to close file [$NOT_RBB_file]: $!\n");
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "a:s", "b:s" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0") if !( $options{'a'} && $options{'b'} );

	return \%options;
}

__DATA__

=head1 NAME

   reciprocal_best_blast_hits.pl

=head1 COPYRIGHT

   Copyright (c) 2016, Jiri Stiller. All rights reserved.

=head1 DESCRIPTION

   This script parses 2 input files - two reciprocal blast reports. 
   In the first step it retrieves hit names for each gene in each gene set and stores them in separated hashes.
   In the second step for each gene set 3 files, NO_HIT, RBB (reciprocal best blast) and NOT_RBB are generated containing the corresponding genes.

=head1 SYNOPSIS

   reciprocal_best_blast_hits.pl -a <blast_report_1> -b <blast_report_2> [-help] [-man]

      [-help]     Displays basic usage information
      [-man]      Displays more detailed information
      -a          blast report input 1
      -b          blast report input 2
       
   Example:   
      reciprocal_best_blast_hits.pl -a sb_ncbi_vs_sb_phytozome.blastplus_report -b sb_phytozome_vs_sb_ncbi.blastplus_report
      
   Output:
      NO_HIT_a_sb_ncbi_vs_sb_phytozome.txt
      NO_HIT_b_sb_phytozome_vs_sb_ncbi.txt
      NOT_RBB_a_sb_ncbi_vs_sb_phytozome.csv
      NOT_RBB_b_sb_phytozome_vs_sb_ncbi.csv
      RBB_a_sb_ncbi_vs_sb_phytozome.csv
      RBB_b_sb_phytozome_vs_sb_ncbi.csv
   
=cut
