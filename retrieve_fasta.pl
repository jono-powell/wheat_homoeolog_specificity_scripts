#!/usr/bin/perl

use strict;
use warnings;

use IO::Handle;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;
use File::Find;
use Carp;

BEGIN {
	*STDOUT->autoflush();
	*STDERR->autoflush();
}
my $options = check_params();

my $found_ids = 0;
my %ids = ();

retrieve_ids_and_create_outputfile();
my $total_ids = keys %ids;
my $output_file;
retrive_fastas();

print "*** Did not find any ids!!!" if !$found_ids;
print "*** Total ids to be retrieved [$total_ids]\n";
print "*** Found                     [$found_ids]\n";

#--------------------------------------------------------

sub retrieve_ids_and_create_outputfile {

	open my $ids_fh, '<', $options->{'ids'}
	  or croak("Can't open ids input file [$options->{'ids'}]: $!.\n");

	while ( my $line = <$ids_fh> ) {
		chomp $line;
		$ids{$line}++;
	}

	close $ids_fh or croak("Failed to close [$options->{'ids'}]: $!\n");

	my $basename = basename( $options->{'ids'}, '.txt' );
	$output_file = $basename . '.fa';

}

sub retrive_fastas {
	my $seqOUT =
	  Bio::SeqIO->new( '-format' => 'fasta', '-file' => ">$output_file" );

	my $seqIN = Bio::SeqIO->new(
		'-format' => 'fasta',
		'-file'   => $options->{'f'}
	);

  FASTA:
	while ( my $fasta_obj = $seqIN->next_seq() ) {
		my $id = $fasta_obj->display_id();

		foreach my $id_to_find ( keys %ids ) {

			if ( $id =~ m/\A$id_to_find\z/ ) {
				$found_ids++;
				$seqOUT->write_seq($fasta_obj);
				last FASTA if ( $found_ids == $total_ids );
				next FASTA;
			}
		}
	}
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "re:s", "f:s", "ids:s" );

	exec("pod2usage $0") if $options{'help'};
	exec("perldoc $0")   if $options{'man'};
	exec("pod2usage $0") if !( $options{'f'} && $options{'ids'} );

	return \%options;
}

__DATA__

=head1 NAME

   retrieve_fastas.pl

=head1 COPYRIGHT

   Copyright (c) 2016, Jiri Stiller. All rights reserved.

=head1 DESCRIPTION

   This scripts reads in fasta ids from txt file (one id per line) 
   and then retrieves correspondings fastas from mutiple fasta file
   and writes them to a new fasta file.

=head1 SYNOPSIS

   retrieve_fastas.pl -f <file_name>  -ids <file name> [-help] [-man]

      -f         Input multiple fasta file from which fastas will be retrieved based on ids
      -ids       txt file containing fasta ids we want to retrieve                     
                      
      Example:
      retreive_fastas.pl -f multiple.fa -ids ids.txt    
         
      Output file:
      ids.fa
         
=cut

