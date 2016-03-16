#!/usr/bin/perl

use strict;
use warnings;

use File::Find;

use lib '/home/sti089/code/perl/lib';

use Getopt::Long;
use Carp;
use File::Basename;
use File::Spec;
use File::Path;
use GetConfig;
use Bio::SeqIO;

BEGIN {
	select(STDERR);
	$| = 1;
	select(STDOUT);
	$| = 1;
}
my $config_file = '/home/sti089/code/perl/homeo/allhits.config';
my $config      = GetConfig::get_config($config_file);

my $options = check_params();

chdir( $options->{'dir'} )
  or croak "Can't change to dir [$options->{'dir'}]: $!.\n";

my @seq_files = glob "$options->{'dir'}/*.fa";

my $count = 0;
print "\n**************************************************\n";
print "*** TOTAL NUMBER OF SEQUENCES TO BE PROCESSED: [", scalar(@seq_files),
  "] ***\n";
print "****************************************************\n";

my $out_file = 'summary.csv';
open my $out_fh, '>', $out_file
  or croak("Can't open output file [$out_file]: $!.\n");
print $out_fh
"gene_id,total_seqs,cons_min,total_good_seqs,total_bad_seqs,good_ids,bad_ids\n";

SEQ:
foreach my $seq_file (@seq_files) {
	$count++;

	print "\n********************************************************\n";
	print "*** PROCESSING SEQUENCE No: [$count] FILE: [$seq_file] ***\n";
	print "**********************************************************\n";

	my $file_basename = fileparse( $seq_file, '.fa' );
	my $total_seqs = `grep -c '>' $seq_file`;
	chomp $total_seqs;
	print "****** TOTAL SEQS [$total_seqs]\n";

	next SEQ if !( $total_seqs > 1 );

	my %bad_seqs = ();

	for my $pos ( 2 .. $total_seqs ) {
		print "*** POS [$pos]\n";
		my $input_clustalw = $file_basename . '_' . $pos . '.fa';
		my $seqIN =
		  Bio::SeqIO->new( '-format' => 'fasta', '-file' => $seq_file );
		my $seqOUT = Bio::SeqIO->new(
			'-format' => 'fasta',
			'-file'   => ">$input_clustalw"
		);

		my $count_seq = 0;

		while ( my $fasta_obj = $seqIN->next_seq() ) {
			$count_seq++;
			last if ( $count_seq > $pos );
			next if $bad_seqs{$count_seq};
			$seqOUT->write_seq($fasta_obj);
		}

		my $total_seqs_pos = `grep -c '>' $input_clustalw`;
		chomp $total_seqs_pos;
		my $clustalw_basename = fileparse( $input_clustalw, '.fa' );
		execute_clustalw($clustalw_basename);
		execute_gblocks( $clustalw_basename, $total_seqs_pos );
		convert_pir2fa( $clustalw_basename, \%bad_seqs, $pos );

		if ( $options->{'del'} ) {
			my @seq_files_to_del =
			  glob "$options->{'dir'}/${clustalw_basename}*";

			foreach my $file (@seq_files_to_del) {
				unlink($file) || carp("Can't unlink $file: $!");
			}
		}

	}

	########## Final run with good seqs
	my $input_clustalw = $file_basename . '_' . 'final' . '.fa';
	my $seqIN  = Bio::SeqIO->new( '-format' => 'fasta', '-file' => $seq_file );
	my $seqOUT = Bio::SeqIO->new(
		'-format' => 'fasta',
		'-file'   => ">$input_clustalw"
	);

	my $count_seq = 0;
	my @good_ids  = ();
	my @bad_ids   = ();

	while ( my $fasta_obj = $seqIN->next_seq() ) {
		$count_seq++;
		my $id = $fasta_obj->display_id();

		if ( $bad_seqs{$count_seq} ) {
			push @bad_ids, $id;
			next;
		}
		else {
			push @good_ids, $id;
			$seqOUT->write_seq($fasta_obj);
		}
	}

	my $total_seqs_pos = `grep -c '>' $input_clustalw`;
	chomp $total_seqs_pos;

	if ( $total_seqs_pos > 1 ) {
		my $clustalw_basename = fileparse( $input_clustalw, '.fa' );
		execute_clustalw($clustalw_basename);
		execute_gblocks( $clustalw_basename, $total_seqs_pos );
		convert_pir2fa( $clustalw_basename, \%bad_seqs );
	}

	write_summary( $file_basename, $total_seqs, \@good_ids, \@bad_ids );

	print "*** TOTAL NUMBER OF PROCESSED SEQUENCES: [$count] ***\n";
	print "**********************************************************\n";
}

close $out_fh or croak("Failed to close [$out_file]: $!\n");

print "*** Concatenating conserved sequences for all genes...\n";
system("cat *cons.fasta > ../concat_cons.fa");
print "\tDone.\n";

#----------------------------------------------------------------------------------
sub write_summary {
	my ( $file_basename, $total_seqs, $good_ids_ref, $bad_ids_ref ) = @_;

	my $total_good = scalar @$good_ids_ref;
	my $total_bad  = scalar @$bad_ids_ref;
	my $good_ids   = join( ' ', @$good_ids_ref );
	my $bad_ids    = join( ' ', @$bad_ids_ref );
	my $cons       = $options->{'cons_min'} ? $options->{'cons_min'} : 'ND';
	print $out_fh
"$file_basename,$total_seqs,$cons,$total_good,$total_bad,$good_ids,$bad_ids\n";

}

sub convert_pir2fa {
	my ( $basename, $bad_seqs_ref, $pos ) = @_;
	my $input_file = $basename . '.pir-gb1';

	my $type = $options->{'cons_min'} ? '_' . $options->{'cons_min'} : '';
	my $output_file = $basename . $type . '_cons.fasta';
	my $seqIN = Bio::SeqIO->new( '-format' => 'pir', '-file' => $input_file );
	my $seqOUT =
	  Bio::SeqIO->new( '-format' => 'fasta', '-file' => ">$output_file" );

	my $delete_file = 0;

  FASTA:
	while ( my $fasta_obj = $seqIN->next_seq() ) {
		my $seq    = $fasta_obj->seq();
		my $length = length($seq);

		if ( $options->{'cons_min'} ) {

			if ( $length < $options->{'cons_min'} ) {
				print "**** Too short cons region!!! [$length]\n";
				$bad_seqs_ref->{$pos}++ if $pos;
				last FASTA;
			}
		}

		$seq =~ s/-//gxms;
		$fasta_obj->seq($seq);
		$seqOUT->write_seq($fasta_obj);
	}
}

sub execute_clustalw {
	my ($basename) = @_;

	my $clustalw_cmd =
	    "$config->{'clustalw_exec'} "
	  . "$config->{'clustalw_opts'} "
	  . '-infile='
	  . "$basename" . '.fa';

	my $search_type = 'clustalw';
	print_start_execution( $search_type, $clustalw_cmd );

	system($clustalw_cmd) == 0
	  or croak "Failed  command [$clustalw_cmd]: $!.\n";
	print_end_execution( $search_type, $clustalw_cmd );
}

sub execute_gblocks {
	my ( $basename, $total_seqs_pos ) = @_;

	my $gblocks_cmd =
	    "$config->{'gblocks_exec'} "
	  . "$basename" . '.pir'
	  . ' -t=d -e=-gb1 -b2='
	  . $total_seqs_pos . ' -b1='
	  . $total_seqs_pos
	  . ' -b4=5 -d=y -b5=a';
	my $search_type = 'gblocks';

	print_start_execution( $search_type, $gblocks_cmd );
	system($gblocks_cmd);
	print_end_execution( $search_type, $gblocks_cmd );
}

sub print_start_execution {
	my ( $search_type, $cmd ) = @_;

	print "\n*** Executing $search_type command:\n[$cmd]\n";
}

sub print_end_execution {
	my ( $search_type, $cmd ) = @_;
	print "*** Successfully finished $search_type command\n";
}

sub check_params {
	my @standard_options = ( "help+", "man+" );
	my %options;

	# Add any other command line options, and the code to handle them
	GetOptions( \%options, @standard_options, "dir:s", "cons_min:i", "del+" );

	# If the -help option is set, print the usage and exit
	exec("pod2usage $0") if $options{'help'};

	# If the -man option is set, run perldoc for this
	exec("perldoc $0") if $options{'man'};

	exec("pod2usage $0") if !( $options{'dir'} );

	if ( $options{'dir'} ) {
		if ( $options{'dir'} !~ m#\A \/ #xms ) {
			croak "dir option must be a full path\n";
		}
	}

	return \%options;
}

__DATA__

=head1 NAME

   pipeline_allhits.pl

=head1 COPYRIGHT

   Copyright (c) 2016, Jiri Stiller. All rights reserved.
   
=head1 DESCRIPTION

   This script processes all sequences for each gene stored in <gene>.fa file. 
   The retrieved sequences in <gene>.fa file are ordered according to the blast result starting with the best hit.
   The script runs clustalw, then gblocks and finally converts .pir to fasta format.
   Option -cons_min specifies the minimum size of the conserved region.
   The comparison starts with the first 2 sequences. It then loops through remaining sequences, 
   adding sequences one by one and recording sequences which don't satisfy minimum conserved region.
   In the end the final comparison is run with only sequences that share conserved region (good sequences). 
   summary.csv files provides details about every gene - total of sequences to compare, good and bad sequences.
   -dir option specifies where the files to be processed are stored. Full path must be specified.  
   The final result and intermediate files are stored in the same directory as -dir option. 
   -del options deletes files of intermediate comparison steps.
   The script does cd to directory specified by -dir option so it can be run from anywhere.
   Approx 500 files are processed per hour (each file has usually 6-10 seqs when blast is run with 10 hits setting).
   If there are only 2 sequences in the <gene>.fa file and no conserved region is found between them then only <gene>_final.fa 
   file is genereated containing only the first sequence.
   As a last step the script concatenates conserved sequences for all genes into a new file 'concat_cons.fa' 
   located in a one directory above directory specified by -dir option (../dir).
   This script requires a config file allhits.config and library GetConfig.pm 
   allhits.config and pipeline_allhits.pl must be modified according to GetConfig.pn, gblocks and clustalw locations.    
 
=head1 SYNOPSIS

   pipeline_allhits.pl -dir <dir_name> [-cons_min <value>] [-del] [-help] [-man]
   
   -dir  <dir name>      Dir containing sequence files to be analysed
   -cons_min             Min size of conserved region
   -del                  Delete temporary files                           
   -help                 Displays basic usage information
   -man                  Displays more detailed information
   
   Example:
       pipeline_allhits.pl -dir /data/tim/DE_genes_Chara_E86_Fp_vs_Mock_Bowtie2_178K_Unigenes/allhits_exon_300 -cons_min 200 -del
   
   Results:
       in directory: 
            allhits_exon_300/
       files inside above directory for each gene, e.g., Ta#S12890914.fa:
            Ta#S12890914.fa
            Ta#S12890914_final_200_cons.fasta
            Ta#S12890914_final.dnd
            Ta#S12890914_final.fa
            Ta#S12890914_final.pir
            Ta#S12890914_final.pir-gb1
            Ta#S12890914_final.pir-gb1.htm
            Ta#S12890914_final.pir-gb1PS

       plus summary file:
            summary.csv
            
=cut




