package GetConfig;
{
	use Exporter;
	@EXPORT_OK = qw(get_config);

	use strict;
	use Memoize;
	use Carp;
	use Safe;

	memoize('GetConfig::get_config');

	sub get_config {

		my ($config_file) = shift @_;

		my $compartment = new Safe;
		$compartment->permit(
			qw(:base_core concat padany anonlist anonhash refgen));

		if ( !-f $config_file ) {
			croak
			  "Error parsing [$config_file]: file does not exist.";
		}

		if ( !$compartment->rdo($config_file) ) {
			croak "Error parsing '$config_file': $@.";
		}

		my %config = $compartment->reval('return(%CONFIG)');
		if ($@) {
			croak "Error getting config file [$config_file]: $@.";
		}

		return ( \%config );
	}

	sub get_var {
		my ( $config_ref, $key ) = @_;

		if ( !defined $config_ref->{$key} ) {
			croak "Variable [$key] not defined.";
		}
		return $config_ref->{$key};
	}
}

__DATA__

=head1 NAME

   GetConfig.pm

=head1 COPYRIGHT

   Copyright (c) 2016, Jiri Stiller. All rights reserved.
   
=head1 DESCRIPTION

=head1 SYNOPSIS

=cut
