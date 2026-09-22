package Mods::Meth2Rep;

use strict;
use warnings;
use Exporter qw(import);
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;
use IO::Uncompress::Gunzip qw($GunzipError);

use Mods::GenoMetaAss qw(getAssemblPath parseSupportReads);

our @EXPORT_OK = qw(
	read_target_mgs read_modbam_manifest read_mgs_report
	build_meth2rep_plan representative_fasta_path
);

sub _trim {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/^\s+|\s+$//g;
	return $value;
}

sub _add_target {
	my ($targets, $value, $origin) = @_;
	$value = _trim($value);
	return if $value eq '' || $value =~ /^#/;
	die "Invalid MGS identifier '$value' in $origin; identifiers may contain only letters, numbers, dots, underscores and hyphens\n"
		unless $value =~ /\A[A-Za-z0-9_.-]+\z/;
	$targets->{$value} = 1;
}

sub read_target_mgs {
	my (%options) = @_;
	my %targets;
	for my $value (split /,/, ($options{mgs} // '')) {
		_add_target(\%targets, $value, '--mgs');
	}
	if (defined($options{mgs_file}) && $options{mgs_file} ne '') {
		open my $fh, '<', $options{mgs_file}
			or die "Cannot read MGS selection file $options{mgs_file}: $!\n";
		while (my $line = <$fh>) {
			$line =~ s/[\r\n]+$//;
			next if $line =~ /^\s*#/;
			$line =~ s/\s+#.*$//;
			for my $value (split /[\t,\s]+/, $line) {
				_add_target(\%targets, $value, $options{mgs_file});
			}
		}
		close $fh or die "Cannot close MGS selection file $options{mgs_file}: $!\n";
	}
	die "Select at least one explicit MGS with --mgs and/or --mgs-file\n"
		unless keys %targets;
	return \%targets;
}

sub read_modbam_manifest {
	my ($file) = @_;
	die "A non-empty --modbam-manifest is required\n"
		unless defined($file) && -s $file;
	open my $fh, '<', $file or die "Cannot read modBAM manifest $file: $!\n";
	my (@header, %column, %entries);
	my $line_number = 0;
	while (my $line = <$fh>) {
		$line_number++;
		$line =~ s/[\r\n]+$//;
		next if $line =~ /^\s*$/;
		if (!@header) {
			$line =~ s/^#//;
			@header = map { lc _trim($_) } split /\t/, $line, -1;
			$column{$header[$_]} = $_ for 0 .. $#header;
			for my $required (qw(sample scope technology modbam)) {
				die "modBAM manifest $file is missing required column '$required'\n"
					unless exists $column{$required};
			}
			next;
		}
		next if $line =~ /^\s*#/;
		my @fields = split /\t/, $line, -1;
		my %row = map {
			$_ => _trim($fields[$column{$_}] // '')
		} qw(sample scope technology modbam);
		die "modBAM manifest $file line $line_number has an empty sample\n"
			if $row{sample} eq '';
		die "modBAM manifest $file line $line_number has invalid scope '$row{scope}' (expected primary or support)\n"
			unless $row{scope} eq 'primary' || $row{scope} eq 'support';
		$row{technology} = uc $row{technology};
		die "modBAM manifest $file line $line_number has technology '$row{technology}'; only ONT and PB carry supported methylation tags\n"
			unless $row{technology} eq 'ONT' || $row{technology} eq 'PB';
		die "modBAM manifest $file line $line_number has an empty modbam path\n"
			if $row{modbam} eq '';
		$row{modbam} = File::Spec->rel2abs($row{modbam}, dirname(File::Spec->rel2abs($file)));
		my $key = "$row{sample}\t$row{scope}";
		$row{line} = $line_number;
		$entries{$key} ||= [];
		for my $existing (@{$entries{$key}}) {
			die "modBAM manifest $file repeats donor '$row{modbam}' for sample '$row{sample}' scope '$row{scope}'\n"
				if $existing->{modbam} eq $row{modbam};
			die "modBAM manifest $file mixes technologies for sample '$row{sample}' scope '$row{scope}'\n"
				if $existing->{technology} ne $row{technology};
		}
		push @{$entries{$key}}, \%row;
	}
	close $fh or die "Cannot close modBAM manifest $file: $!\n";
	die "modBAM manifest $file contains no data rows\n" unless keys %entries;
	return \%entries;
}

sub read_mgs_report {
	my ($file, $selected) = @_;
	die "MGS MAG report is missing or empty: $file\n" unless defined($file) && -s $file;
	my $fh = IO::Uncompress::Gunzip->new($file, Transparent => 1, MultiStream => 1)
		or die "Cannot read MGS MAG report $file: $GunzipError\n";
	my (@header, %column, %result, %seen_mag);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '';
		my @fields = split /\t/, $line, -1;
		if (!@header) {
			@header = @fields;
			$column{$header[$_]} = $_ for 0 .. $#header;
			for my $required (qw(MAG MGS Representative4MGS)) {
				die "MGS MAG report $file is missing required column '$required'\n"
					unless exists $column{$required};
			}
			next;
		}
		my $mgs = $fields[$column{MGS}] // '';
		next unless exists $selected->{$mgs};
		my $mag = $fields[$column{MAG}] // '';
		die "MGS MAG report $file has an empty MAG for selected $mgs\n" if $mag eq '';
		die "MGS MAG report $file repeats MAG '$mag' in selected $mgs\n"
			if $seen_mag{"$mgs\t$mag"}++;
		my $is_canopy = $mag =~ /^Cano__/ ? 1 : 0;
		push @{$result{$mgs}{rows}}, {
			mag => $mag,
			is_canopy => $is_canopy,
			representative => (($fields[$column{Representative4MGS}] // '') eq '*') ? 1 : 0,
		};
	}
	close $fh or die "Cannot close MGS MAG report $file: $!\n";

	for my $mgs (sort keys %{$selected}) {
		die "Requested MGS '$mgs' was not found in $file\n"
			unless exists $result{$mgs};
		my @mags = grep { !$_->{is_canopy} } @{$result{$mgs}{rows}};
		if (!@mags) {
			$result{$mgs}{available} = 0;
			$result{$mgs}{reason} = 'canopy_only_no_mag_reference';
			next;
		}
		my @representatives = grep { $_->{representative} } @mags;
		die "Selected MGS '$mgs' has " . scalar(@representatives)
			. " non-Canopy representatives in $file; exactly one is required\n"
			unless @representatives == 1;
		$result{$mgs}{available} = 1;
		$result{$mgs}{representative_mag} = $representatives[0]{mag};
		$result{$mgs}{member_mags} = [map { $_->{mag} } @mags];
	}
	return \%result;
}

sub _mag_parts {
	my ($mag) = @_;
	die "MAG identifier '$mag' does not have the required sample__bin form\n"
		unless defined($mag) && $mag =~ /\A(.+)__(.+)\z/;
	return ($1, $2);
}

sub representative_fasta_path {
	my ($directory, $mgs, $mag) = @_;
	my $name = "$mgs.ctgs.$mag";
	$name =~ s/\.gz\z//i;
	$name .= '.fna' unless $name =~ /\.(?:fa|fna|fasta)\z/i;
	$name .= '.gz';
	return File::Spec->catfile($directory, $name);
}

sub _read_bin_assignments {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read bin assignments $file: $!\n";
	my (%bins, %seen_contig);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line =~ /^\s*$/;
		my @fields = split /\t/, $line, -1;
		die "Malformed bin assignment in $file at line $.\n"
			unless @fields >= 2 && $fields[0] ne '' && $fields[1] ne '';
		next if $fields[0] eq 'Sequence ID';
		if (exists($seen_contig{$fields[0]}) && $seen_contig{$fields[0]} ne $fields[1]) {
			die "Contig '$fields[0]' has conflicting bins '$seen_contig{$fields[0]}' and '$fields[1]' in $file\n";
		}
		next if exists $seen_contig{$fields[0]};
		$seen_contig{$fields[0]} = $fields[1];
		push @{$bins{$fields[1]}}, $fields[0];
	}
	close $fh or die "Cannot close bin assignments $file: $!\n";
	return \%bins;
}

sub _canonical_sample {
	my ($map, $sample) = @_;
	return $sample if exists $map->{$sample};
	return $map->{altNms}{$sample}
		if ref($map->{altNms}) eq 'HASH' && exists($map->{altNms}{$sample});
	die "MAG sample '$sample' is not present in the resolved catalogue mapping metadata\n";
}

sub _scope_technology {
	my ($record, $scope, $sample) = @_;
	if ($scope eq 'primary') {
		return '' unless $record->{hasPrimaryRds};
		return uc($record->{SeqTech} // '');
	}
	return '' unless defined($record->{SupportReads}) && $record->{SupportReads} ne '';
	my ($technology) = parseSupportReads($record->{SupportReads});
	return uc($technology // '');
}

sub _cram_for_scope {
	my ($record, $sample, $scope) = @_;
	my $suffix = $scope eq 'support' ? '.sup-smd.cram' : '-smd.cram';
	return File::Spec->catfile($record->{wrdir}, 'mapping', "$sample$suffix");
}

sub build_meth2rep_plan {
	my (%options) = @_;
	my $map = $options{map};
	my $groups = $options{assembly_groups};
	my $manifest = $options{manifest};
	my $selected = $options{selected};
	my $modes = $options{modes};
	die "Mapping metadata are required\n" unless ref($map) eq 'HASH';
	die "Assembly-group metadata are required\n" unless ref($groups) eq 'HASH';
	die "modBAM manifest data are required\n" unless ref($manifest) eq 'HASH';
	die "At least one meth2rep mode is required\n" unless ref($modes) eq 'ARRAY' && @{$modes};
	my %valid_mode = map { $_ => 1 } qw(mgs2rep rep2rep);
	die "Unknown meth2rep mode in " . join(',', @{$modes}) . "\n"
		if grep { !$valid_mode{$_} } @{$modes};

	my $report = read_mgs_report($options{mgs_report}, $selected);
	my (%requested_groups, %resolved_groups, %inputs);
	$inputs{$options{mgs_report}} = 1;
	$inputs{$options{manifest_file}} = 1 if defined $options{manifest_file};
	my %requested_mode = map { $_ => 1 } @{$modes};
	for my $mgs (sort keys %{$selected}) {
		next unless $report->{$mgs}{available};
		my @source_mags = $requested_mode{mgs2rep}
			? @{$report->{$mgs}{member_mags}}
			: ($report->{$mgs}{representative_mag});
		for my $mag (@source_mags) {
			my ($mag_sample) = _mag_parts($mag);
			$mag_sample = _canonical_sample($map, $mag_sample);
			my $group = $map->{$mag_sample}{AssGroup};
			die "MAG '$mag' has no assembly group in mapping metadata\n"
				unless defined($group) && $group ne '';
			die "MAG '$mag' refers to assembly group '$group', which is absent from assembly-group metadata\n"
				unless exists $groups->{$group};
			$requested_groups{$group} = 1;
		}
	}
	for my $group (sort keys %requested_groups) {
		my @samples = @{$groups->{$group}{SmplID} || []};
		my @workdirs = @{$groups->{$group}{wrdir} || []};
		die "Assembly group '$group' has inconsistent sample/work-directory metadata\n"
			unless @samples == @workdirs && @samples;
		my ($assembly, @eligible_samples);
		for my $index (0 .. $#samples) {
			my $sample = $samples[$index];
			next if -e File::Spec->catfile($workdirs[$index], 'SMPL.empty');
			my $candidate = getAssemblPath($workdirs[$index]);
			die "Assembly group '$group' sample '$sample' has no resolved assembly directory\n"
				unless defined($candidate) && -d $candidate;
			$candidate = abs_path($candidate);
			if (defined($assembly) && $candidate ne $assembly) {
				die "Assembly group '$group' resolves to multiple assemblies: '$assembly' and '$candidate'\n";
			}
			$assembly = $candidate;
			push @eligible_samples, $sample;
		}
		next unless @eligible_samples;
		my $reference = File::Spec->catfile($assembly, 'scaffolds.fasta.filt');
		die "Assembly reference is missing or empty for group '$group': $reference\n" unless -s $reference;
		$inputs{$reference} = 1;
		$resolved_groups{$group} = {
			assembly => $assembly,
			reference => $reference,
			samples => \@eligible_samples,
		};
	}

	my (%assignment_cache, %targets, %contig_targets, @unavailable, %needed_groups);
	for my $mgs (sort keys %{$selected}) {
		if (!$report->{$mgs}{available}) {
			for my $mode (@{$modes}) {
				push @unavailable, { mgs => $mgs, mode => $mode, reason => $report->{$mgs}{reason} };
			}
			next;
		}
		my $representative = $report->{$mgs}{representative_mag};
		my $representative_fasta = representative_fasta_path(
			$options{representatives_dir}, $mgs, $representative,
		);
		die "Representative contig FASTA is missing for selected $mgs ($representative): $representative_fasta\n"
			unless -s $representative_fasta;
		$representative_fasta = abs_path($representative_fasta);
		$inputs{$representative_fasta} = 1;
		for my $mode (@{$modes}) {
			my $target_key = "$mode\t$mgs";
			my @source_mags = $mode eq 'rep2rep'
				? ($representative) : @{$report->{$mgs}{member_mags}};
			my $target = $targets{$target_key} = {
				key => $target_key, mode => $mode, mgs => $mgs,
				representative_mag => $representative,
				reference => $representative_fasta,
				sources => [], groups => {}, scope_keys => {}, available => 1,
			};
			for my $mag (@source_mags) {
				my ($mag_sample, $bin) = _mag_parts($mag);
				$mag_sample = _canonical_sample($map, $mag_sample);
				my $group = $map->{$mag_sample}{AssGroup};
				die "MAG '$mag' has no eligible resolved assembly for group '$group'\n"
					unless defined($group) && exists $resolved_groups{$group};
				my $assembly = $resolved_groups{$group}{assembly};
				my $assignment = File::Spec->catfile(
					$assembly, 'Binning', $options{binner}, $mag_sample,
				);
				die "Bin assignment file is missing for MAG '$mag': $assignment\n"
					unless -s $assignment;
				$assignment = abs_path($assignment);
				$inputs{$assignment} = 1;
				$assignment_cache{$assignment} ||= _read_bin_assignments($assignment);
				my $contigs = $assignment_cache{$assignment}{$bin};
				die "MAG '$mag' bin '$bin' has no contigs in $assignment\n"
					unless ref($contigs) eq 'ARRAY' && @{$contigs};
				push @{$target->{sources}}, {
					mag => $mag, group => $group, assembly => $assembly,
					assignment => $assignment, bin => $bin,
					contigs => [@{$contigs}],
				};
				$target->{groups}{$group} = 1;
				$needed_groups{$group} = 1;
				for my $contig (@{$contigs}) {
					$contig_targets{$assembly}{$contig}{$target_key} = 1;
				}
			}
		}
	}

	my (%scopes, %manifest_used);
	for my $group (sort keys %needed_groups) {
		for my $sample (@{$resolved_groups{$group}{samples}}) {
			die "Assembly-group sample '$sample' is absent from mapping metadata\n"
				unless exists $map->{$sample};
			for my $scope (qw(primary support)) {
				my $technology = _scope_technology($map->{$sample}, $scope, $sample);
				next unless $technology eq 'ONT' || $technology eq 'PB';
				my $manifest_key = "$sample\t$scope";
				die "No original modBAM is declared for long-read sample '$sample' scope '$scope' (technology $technology)\n"
					unless exists $manifest->{$manifest_key};
				my (@modbams, %seen_modbam);
				for my $entry (@{$manifest->{$manifest_key}}) {
					die "modBAM manifest technology '$entry->{technology}' does not match mapping technology '$technology' for sample '$sample' scope '$scope'\n"
						unless $entry->{technology} eq $technology;
					die "Original modBAM is missing or empty for selected sample '$sample' scope '$scope': $entry->{modbam}\n"
						unless -s $entry->{modbam};
					my $canonical_modbam = abs_path($entry->{modbam});
					die "Original modBAM '$entry->{modbam}' repeats the same physical donor for sample '$sample' scope '$scope'\n"
						if $seen_modbam{$canonical_modbam}++;
					push @modbams, $canonical_modbam;
				}
				my $cram = _cram_for_scope($map->{$sample}, $sample, $scope);
				die "Required assembly backmapping CRAM is missing or empty for sample '$sample' scope '$scope': $cram\n"
					unless -s $cram;
				die "Assembly backmapping CRAM has no completion stone and may be partial: $cram.sto\n"
					unless -e "$cram.sto";
				$cram = abs_path($cram);
				my $scope_key = "$sample\t$scope";
				$scopes{$scope_key} = {
					key => $scope_key, sample => $sample, scope => $scope,
					technology => $technology, modbams => \@modbams,
					cram => $cram, assembly => $resolved_groups{$group}{assembly},
					assembly_reference => $resolved_groups{$group}{reference},
					group => $group,
				};
				$manifest_used{$manifest_key} = 1;
				$inputs{$cram} = 1;
				$inputs{"$cram.sto"} = 1;
				(my $reference_stat = $cram) =~ s/\.cram\z/.reference.stat/;
				$inputs{$reference_stat} = 1 if -e $reference_stat;
				$inputs{$_} = 1 for @modbams;
				for my $target (values %targets) {
					next unless $target->{groups}{$group};
					$target->{scope_keys}{$scope_key} = 1;
				}
			}
		}
	}

	for my $target_key (sort keys %targets) {
		my $target = $targets{$target_key};
		next if keys %{$target->{scope_keys}};
		$target->{available} = 0;
		$target->{reason} = 'no_eligible_ont_or_pb_scope';
		push @unavailable, {
			mgs => $target->{mgs}, mode => $target->{mode}, reason => $target->{reason},
		};
	}

	return {
		report => $report,
		targets => \%targets,
		scopes => \%scopes,
		contig_targets => \%contig_targets,
		unavailable => \@unavailable,
		inputs => [sort grep { defined($_) && $_ ne '' } keys %inputs],
		manifest_unused => [sort grep { !$manifest_used{$_} } keys %{$manifest}],
	};
}

1;
