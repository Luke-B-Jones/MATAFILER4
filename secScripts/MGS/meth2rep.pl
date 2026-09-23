#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../..";
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Spec;
use File::Basename qw(dirname);
use File::Find qw(find);
use Cwd qw(abs_path);
use File::Temp qw(tempdir);
use Digest::SHA qw(sha1_hex sha256_hex);
use DB_File;
use Fcntl qw(O_CREAT O_RDWR LOCK_EX LOCK_NB);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP;
use Text::ParseWords qw(shellwords);
use Time::HiRes ();

use Mods::Checkpoint qw(write_checkpoint checkpoint_valid read_checkpoint);
use Mods::GenoMetaAss qw(getDirsPerAssmblGrp parseSupportReads);
use Mods::Meth2Rep qw(
	read_target_mgs read_modbam_manifest build_meth2rep_plan
);

my $VERSION = '0.7';
my $HISTORICAL_SOURCE_FILTER =
	'inherited_unknown: the retained assembly CRAM was filtered when MATAFILER created it, but legacy CRAM completion stones do not record the run-specific filter parameters';
my @publication_partials;
END {
	for my $file (@publication_partials) {
		unlink $file if defined($file) && -f $file;
	}
}

sub usage {
	return <<'USAGE';
Usage:
  meth2rep.pl --mgs-dir DIR --map FILE[,FILE...]
    --modbam-manifest FILE (--mgs MGS.1,MGS.2 | --mgs-file FILE)
    (--mgs2rep | --rep2rep | both) [-o DIR]

Modes (either or both):
  --mgs2rep  Candidate reads supporting any MAG in an MGS are freshly aligned
             to that MGS representative MAG.
  --rep2rep  Only candidate reads supporting the representative MAG are freshly
             aligned to that same representative MAG.

Required manifest header (tab separated):
  sample  scope  technology  modbam

scope must be primary or support; technology must be ONT or PB. Repeated
sample/scope rows declare multiple donor BAMs. A QNAME must be unique across
those BAMs for its sample/scope; ambiguous identities stop the run.

Controls:
  --source-min-mapq INT       Additional MAPQ floor on retained CRAM records (default 10)
  --source-min-coverage FLOAT Additional aligned-query floor on retained CRAM records (default 0.5)
  --target-min-mapq INT       Minimum representative MAPQ for both technologies
  --target-min-coverage FLOAT Minimum aligned query fraction for both technologies
  --target-max-edit-rate FLOAT Maximum representative edit rate for both technologies
  --target-min-end-clip INT   Minimum two-ended clipping rejected by bamFilter
  --mapper-filter-ont STRING  bamFilter arguments (default "0.15 0.5 10 0")
  --mapper-filter-pb STRING   bamFilter arguments (default "0.05 0.5 30 0")
  --minimap2-preset-ont NAME  ONT preset (default map-ont; modern option: lr:hq)
  --minimap2-preset-pb NAME   PacBio preset (default map-pb; HiFi option: map-hifi)
  --supplementary-alignments POLICY
                              drop (default) or keep non-overlapping records
  --allow-missing-mn          Permit legacy donors without MN (warning is recorded)
  --output-format FORMAT      bam (default) or self-contained cram
  --threads INT               Mapping/sorting threads (default 4)
  --memory-gb INT             Sort-memory sizing budget, not a process cap (default 32)
  --tmp DIR                   Parent for disposable intermediates
  --keep-read-ids             Retain gzip-compressed candidate-name audit files
  --out-manifest FILE         Completed-unit/donor TSV (default OUT/manifest.tsv)
  --override                  Rebuild selected units even when their checkpoints validate
  --redo                      Alias of --override
  --plan-only                 Run input/tool/output preflight and write a preview plan
  --samtools COMMAND          Override samtools (otherwise resolved from PATH)
  --minimap2 COMMAND          Override minimap2 (otherwise resolved from PATH)
  --bam-filter COMMAND        Override bundled first-party bamFilter.pl
USAGE
}

sub shell_quote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub tool_command {
	my ($tool, @arguments) = @_;
	my @prefix = shellwords($tool // '');
	die "Empty external-tool command\n" unless @prefix;
	return join(' ', map { shell_quote($_) } (@prefix, @arguments));
}

sub find_on_path {
	my ($program) = @_;
	return unless defined($program) && $program ne '';
	for my $directory (split /:/, ($ENV{PATH} // '')) {
		$directory = '.' if $directory eq '';
		my $candidate = File::Spec->catfile($directory, $program);
		return abs_path($candidate) || File::Spec->rel2abs($candidate)
			if -f $candidate && -x $candidate;
	}
	return;
}

sub file_sha256 {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot hash $file: $!\n";
	binmode $fh;
	my $digest = Digest::SHA->new(256);
	$digest->addfile($fh);
	close $fh or die "Cannot close $file after hashing: $!\n";
	return $digest->hexdigest;
}

sub canonical_command {
	my ($command, $label) = @_;
	my @tokens = shellwords($command // '');
	die "$label command is empty\n" unless @tokens;
	my $executable;
	if (File::Spec->file_name_is_absolute($tokens[0]) || $tokens[0] =~ m{/}) {
		my $candidate = File::Spec->rel2abs($tokens[0]);
		$executable = abs_path($candidate) || $candidate;
	} else {
		$executable = find_on_path($tokens[0]);
	}
	die "$label executable '$tokens[0]' is missing or not executable\n"
		unless defined($executable) && -f $executable && -x $executable;
	$tokens[0] = $executable;
	for my $index (1 .. $#tokens) {
		next unless -f $tokens[$index];
		$tokens[$index] = abs_path($tokens[$index]) || File::Spec->rel2abs($tokens[$index]);
	}
	return join(' ', map { shell_quote($_) } @tokens);
}

sub command_identity {
	my ($command) = @_;
	my @tokens = shellwords($command);
	my @parts = ('argv=' . join("\x1f", @tokens));
	for my $token (@tokens) {
		next unless -f $token;
		my @stat = stat($token);
		push @parts, join(':', $token, $stat[7], $stat[9], file_sha256($token));
	}
	return sha256_hex(join("\n", @parts));
}

sub run_shell {
	my ($description, $command) = @_;
	my $status = system('bash', '-o', 'pipefail', '-c', $command);
	return if $status == 0;
	my $detail = $status == -1 ? "could not execute: $!"
		: ($status & 127) ? 'terminated by signal ' . ($status & 127)
		: 'exit code ' . ($status >> 8);
	die "$description failed ($detail)\n";
}

sub open_command {
	my ($description, $command) = @_;
	open my $fh, '-|', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub open_sink_command {
	my ($description, $command) = @_;
	open my $fh, '|-', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub close_command {
	my ($fh, $description) = @_;
	return if close $fh;
	my $status = $?;
	my $detail = $status == -1 ? "could not execute: $!"
		: ($status & 127) ? 'terminated by signal ' . ($status & 127)
		: 'exit code ' . ($status >> 8);
	die "$description failed ($detail)\n";
}

sub capture_command {
	my ($description, $command) = @_;
	my $fh = open_command($description, $command);
	local $/;
	my $output = <$fh> // '';
	close_command($fh, $description);
	$output =~ s/[\r\n]+$//;
	return $output;
}

sub sort_unique {
	my ($input, $output, $tmp_parent) = @_;
	local $ENV{LC_ALL} = 'C';
	my @command = ('sort', '-u');
	push @command, ('-T', $tmp_parent) if defined($tmp_parent) && -d $tmp_parent;
	push @command, ('-o', $output, $input);
	my $status = system @command;
	die "Sorting read names failed (exit " . ($status == -1 ? -1 : $status >> 8) . ")\n"
		if $status != 0;
}

sub count_lines {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot count $file: $!\n";
	my $count = 0;
	$count++ while <$fh>;
	close $fh or die "Cannot close $file: $!\n";
	return $count;
}

sub concatenate_files {
	my ($output, @inputs) = @_;
	open my $out, '>', $output or die "Cannot create $output: $!\n";
	for my $input (@inputs) {
		open my $in, '<', $input or die "Cannot read $input: $!\n";
		while (read($in, my $buffer, 1024 * 1024)) {
			print {$out} $buffer or die "Cannot write $output: $!\n";
		}
		close $in or die "Cannot close $input: $!\n";
	}
	close $out or die "Cannot close $output: $!\n";
}

sub remove_temp_files {
	for my $file (@_) {
		next unless defined($file) && $file ne '' && -e $file;
		unlink $file or die "Cannot remove temporary file $file: $!\n";
	}
}

sub filter_bam_by_names {
	my (%options) = @_;
	my $database = "$options{workdir}/names." . sha1_hex(
		join("\0", $options{names}, $options{input}, $options{output}),
	) . '.db';
	my %wanted;
	tie %wanted, 'DB_File', $database, O_RDWR | O_CREAT, 0600, $DB_HASH
		or die "Cannot create disk-backed read-name index $database: $!\n";
	open my $names, '<', $options{names} or die "Cannot read $options{names}: $!\n";
	while (my $name = <$names>) {
		$name =~ s/[\r\n]+$//;
		$wanted{$name} = 1 if $name ne '';
	}
	close $names or die "Cannot close $options{names}: $!\n";

	my $source_description = "reading $options{input} for name filtering";
	my $source = open_command(
		$source_description,
		tool_command($options{samtools}, 'view', '-h', $options{input}),
	);
	my $sink_description = "writing filtered BAM $options{output}";
	my $sink = open_sink_command(
		$sink_description,
		tool_command($options{samtools}, 'view', '-b', '-o', $options{output}, '-'),
	);
	while (my $line = <$source>) {
		if ($line =~ /^\@/) {
			print {$sink} $line or die "Cannot stream header to $options{output}: $!\n";
			next;
		}
		my ($name, $flag) = split /\t/, $line, 3;
		next unless defined($name) && exists $wanted{$name};
		next if ($options{exclude_flags} || 0) && ($flag & $options{exclude_flags});
		print {$sink} $line or die "Cannot stream record to $options{output}: $!\n";
	}
	close_command($source, $source_description);
	close_command($sink, $sink_description);
	untie %wanted or die "Cannot close disk-backed read-name index $database: $!\n";
	unlink $database or die "Cannot remove temporary read-name index $database: $!\n" if -e $database;
}

sub assert_name_subset {
	my ($wanted_file, $found_file, $context) = @_;
	open my $wanted, '<', $wanted_file or die "Cannot read $wanted_file: $!\n";
	open my $found, '<', $found_file or die "Cannot read $found_file: $!\n";
	my $wanted_name = <$wanted>;
	my $found_name = <$found>;
	chomp $wanted_name if defined $wanted_name;
	chomp $found_name if defined $found_name;
	my (@missing, $missing_count);
	while (defined $wanted_name) {
		while (defined($found_name) && $found_name lt $wanted_name) {
			$found_name = <$found>;
			chomp $found_name if defined $found_name;
		}
		if (!defined($found_name) || $found_name ne $wanted_name) {
			$missing_count++;
			push @missing, $wanted_name if @missing < 5;
		}
		$wanted_name = <$wanted>;
		chomp $wanted_name if defined $wanted_name;
	}
	close $wanted or die "Cannot close $wanted_file: $!\n";
	close $found or die "Cannot close $found_file: $!\n";
	if ($missing_count) {
		die "$context: $missing_count candidate read name(s) are absent from the declared original modBAM set"
			. (@missing ? ' (examples: ' . join(', ', @missing) . ')' : '')
			. ". Refusing a silently incomplete methylation result.\n";
	}
}

sub file_component {
	my ($value) = @_;
	$value =~ s/([^A-Za-z0-9_.-])/sprintf('_%02X', ord($1))/ge;
	return $value;
}

sub native_sequence {
	my ($fields, $origin) = @_;
	my ($name, $flag, $cigar, $sequence) = @{$fields}[0, 1, 5, 9];
	die "$origin record '$name' has no complete SEQ\n"
		if !defined($sequence) || $sequence eq '' || $sequence eq '*';
	die "$origin record '$name' is hard-clipped ($cigar); its full native sequence is unavailable\n"
		if defined($cigar) && $cigar =~ /H/;
	if ($flag & 0x10) {
		$sequence = reverse $sequence;
		$sequence =~ tr/ACGTRYKMSWBDHVNacgtrykmswbdhvn/TGCAYRMKSWVHDBNtgcayrmkswvhdbn/;
	}
	return uc $sequence;
}

sub sort_identity_records {
	my ($input, $output, $tmp_parent) = @_;
	local $ENV{LC_ALL} = 'C';
	my @command = ('sort', '-k1,1');
	push @command, ('-T', $tmp_parent) if defined($tmp_parent) && -d $tmp_parent;
	push @command, ('-o', $output, $input);
	my $status = system @command;
	die "Sorting read-identity records failed (exit " . ($status == -1 ? -1 : $status >> 8) . ")\n"
		if $status != 0;
}

sub identity_record {
	my ($fh, $file) = @_;
	my $line = <$fh>;
	return unless defined $line;
	$line =~ s/[\r\n]+\z//;
	my ($name, $digest, @extra) = split /\t/, $line, -1;
	die "Malformed read-identity record in $file\n"
		unless defined($name) && $name ne '' && defined($digest)
		&& $digest =~ /\A[0-9a-f]{64}\z/ && !@extra;
	return [$name, $digest];
}

sub assert_source_donor_identity {
	my ($source_file, $donor_file, $context) = @_;
	open my $source_fh, '<', $source_file or die "Cannot read $source_file: $!\n";
	open my $donor_fh, '<', $donor_file or die "Cannot read $donor_file: $!\n";
	my $source = identity_record($source_fh, $source_file);
	my $donor = identity_record($donor_fh, $donor_file);
	my ($last_source, $last_donor);
	while (defined $source) {
		die "$context: source CRAM contains more than one primary record named '$source->[0]'\n"
			if defined($last_source) && $last_source eq $source->[0];
		while (defined($donor) && $donor->[0] lt $source->[0]) {
			die "$context: declared donor set contains more than one record named '$donor->[0]'\n"
				if defined($last_donor) && $last_donor eq $donor->[0];
			$last_donor = $donor->[0];
			$donor = identity_record($donor_fh, $donor_file);
		}
		die "$context: candidate read '$source->[0]' is absent from the declared original modBAM set\n"
			unless defined($donor) && $donor->[0] eq $source->[0];
		die "$context: native SEQ for '$source->[0]' differs between the assembly CRAM and declared original modBAM; QNAME alone is not sufficient provenance\n"
			unless $donor->[1] eq $source->[1];
		$last_source = $source->[0];
		$last_donor = $donor->[0];
		$source = identity_record($source_fh, $source_file);
		$donor = identity_record($donor_fh, $donor_file);
	}
	close $source_fh or die "Cannot close $source_file: $!\n";
	close $donor_fh or die "Cannot close $donor_file: $!\n";
}

sub bam_count {
	my ($samtools, $bam) = @_;
	my $count = capture_command(
		"counting records in $bam",
		tool_command($samtools, 'view', '-c', $bam),
	);
	die "samtools returned a non-integer record count for $bam: '$count'\n"
		unless $count =~ /^\d+$/;
	return 0 + $count;
}

sub bam_tag_counts {
	my ($samtools, $bam) = @_;
	my $fh = open_command("reading tags from $bam", tool_command($samtools, 'view', $bam));
	my ($records, $mm, $ml) = (0, 0, 0);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		my @fields = split /\t/, $line, -1;
		die "Malformed transferred SAM record from $bam\n" unless @fields >= 11;
		my (@mm_tags, @ml_tags, @mn_tags);
		for my $field (@fields[11 .. $#fields]) {
			push @mm_tags, $field if $field =~ /\AM[Mm]:/;
			push @ml_tags, $field if $field =~ /\AM[Ll]:/;
			push @mn_tags, $field if $field =~ /\AMN:/;
		}
		die "Transferred record '$fields[0]' in $bam does not contain exactly one MM:Z and one ML:B:C tag\n"
			unless @mm_tags == 1 && @ml_tags == 1
			&& $mm_tags[0] =~ /\AM[Mm]:Z:/ && $ml_tags[0] =~ /\AM[Ll]:B:C(?:,\d+)*\z/;
		die "Transferred record '$fields[0]' in $bam has missing or stale MN for SEQ length "
			. length($fields[9]) . "\n"
			unless @mn_tags == 1 && $mn_tags[0] =~ /\AMN:i:(\d+)\z/
			&& $1 == length($fields[9]);
		$records++;
		$mm++;
		$ml++;
	}
	close_command($fh, "reading tags from $bam");
	return ($records, $mm, $ml);
}

sub validate_donor_subset {
	my ($samtools, $bam, $found_raw, $identity_raw, $allow_missing_mn) = @_;
	my $fh = open_command("validating donor records in $bam", tool_command($samtools, 'view', $bam));
	open my $names, '>', $found_raw or die "Cannot create $found_raw: $!\n";
	my $identities;
	if (defined($identity_raw) && $identity_raw ne '') {
		open $identities, '>', $identity_raw or die "Cannot create $identity_raw: $!\n";
	}
	my ($previous, $records, $mm, $ml, $missing_mn);
	$records = $mm = $ml = $missing_mn = 0;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		my @fields = split /\t/, $line, -1;
		die "Malformed donor SAM record from $bam\n" unless @fields >= 11;
		my ($name, $flag) = @fields[0, 1];
		die "Donor modBAM $bam contains a paired record '$name'; meth2rep only supports singleton ONT/PB records\n"
			if $flag & 0x1;
		die "Donor modBAM $bam contains more than one primary record named '$name'; CRAM-to-donor identity is ambiguous\n"
			if defined($previous) && $previous eq $name;
		my $sequence = native_sequence(\@fields, 'Donor');
		my (@mm_fields, @ml_fields, @mn_fields);
		for my $field (@fields[11 .. $#fields]) {
			push @mm_fields, $field if $field =~ /\AM[Mm]:/;
			push @ml_fields, $field if $field =~ /\AM[Ll]:/;
			push @mn_fields, $field if $field =~ /\AMN:/;
		}
		die "Donor record '$name' must contain exactly one paired MM:Z/ML:B:C tag set\n"
			unless @mm_fields == 1 && @ml_fields == 1
			&& $mm_fields[0] =~ /\AM[Mm]:Z:(?:[ACGTUN][+-](?:[A-Za-z]+|\d+)[.?]?(?:,\d+)*;)*\z/
			&& $ml_fields[0] =~ /\AM[Ll]:B:C(?:,\d+)*\z/;
		die "Donor record '$name' contains duplicate MN tags\n" if @mn_fields > 1;
		if (!@mn_fields) {
			die "Donor record '$name' lacks MN; use --allow-missing-mn only for explicitly reviewed legacy input\n"
				unless $allow_missing_mn;
			$missing_mn++;
		} else {
			die "Donor record '$name' has malformed or stale MN\n"
				unless $mn_fields[0] =~ /\AMN:i:(\d+)\z/ && $1 == length($fields[9]);
		}
		$previous = $name;
		print {$names} "$name\n" or die "Cannot write $found_raw: $!\n";
		if ($identities) {
			print {$identities} "$name\t", sha256_hex($sequence), "\n"
				or die "Cannot write $identity_raw: $!\n";
		}
		$records++;
		$mm++;
		$ml++;
	}
	close $names or die "Cannot close $found_raw: $!\n";
	close $identities or die "Cannot close $identity_raw: $!\n" if $identities;
	close_command($fh, "validating donor records in $bam");
	return ($records, $mm, $ml, $missing_mn);
}

sub validate_filter {
	my ($value, $name) = @_;
	my @parts = split /\s+/, $value;
	die "$name must contain four bamFilter values: max_edit_rate min_query_coverage min_mapq min_end_clip\n"
		unless @parts == 4
		&& $parts[0] =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $parts[0] >= 0 && $parts[0] <= 1
		&& $parts[1] =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $parts[1] >= 0 && $parts[1] <= 1
		&& $parts[2] =~ /^\d+$/ && $parts[2] <= 255
		&& $parts[3] =~ /^\d+$/;
	return @parts;
}

sub cigar_query_coverage {
	my ($cigar, $context) = @_;
	my ($aligned, $length) = (0, 0);
	my $reconstructed = '';
	while ($cigar =~ /([1-9]\d*)([MIDNSHP=X])/g) {
		my ($span, $op) = ($1, $2);
		$reconstructed .= "$span$op";
		if ($op =~ /[MIS=XH]/) { $length += $span; }
		if ($op =~ /[MI=X]/) { $aligned += $span; }
	}
	die "Malformed source CIGAR '$cigar' for $context\n"
		unless $cigar ne '*' && $reconstructed eq $cigar && $length > 0;
	return $aligned / $length;
}

sub fingerprint_inputs {
	my ($inputs, $parameters) = @_;
	my @records;
	for my $file (@{$inputs}) {
		my @stat = stat($file);
		die "Cannot fingerprint meth2rep input $file\n" unless @stat;
		push @records, join("\t", $file, map { 0 + $stat[$_] } (0, 1, 7, 9, 10));
	}
	push @records, map { "option\t$_\t$parameters->{$_}" } sort keys %{$parameters};
	return sha256_hex(join("\n", @records));
}

sub require_mgs_readiness {
	my ($directory) = @_;
	die "MGS directory is missing: $directory\n" unless -d $directory;
	my %stages = (
		Stage1 => 'stage-1',
		BinExtr => 'extract-bin-contigs',
	);
	my %manifests;
	for my $name (sort keys %stages) {
		my $stone = File::Spec->catfile($directory, 'LOGandSUB', 'checkpoints', "$name.stone");
		die "Required MGS progress checkpoint is missing or legacy-empty: $stone\n" unless -s $stone;
		my $manifest = read_checkpoint($stone)
			or die "MGS progress checkpoint is malformed: $stone\n";
		die "MGS progress checkpoint has wrong stage: $stone\n"
			unless ($manifest->{parameters}{stage} // '') eq $stages{$name};
		die "MGS progress checkpoint has stale/missing outputs: $stone\n"
			unless checkpoint_valid($stone, parameters => { stage => $stages{$name} });
		$manifests{$name} = $manifest;
	}
	return \%manifests;
}

sub read_json_file {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read $file: $!\n";
	local $/;
	my $data = eval { JSON::PP->new->decode(<$fh> // '') };
	close $fh or die "Cannot close $file: $!\n";
	return $data if !$@ && ref($data) eq 'HASH';
	die "Invalid meth2rep unit JSON: $file\n";
}

sub write_json_file {
	my ($file, $data) = @_;
	my $partial = "$file.part.$$";
	push @publication_partials, $partial;
	open my $fh, '>', $partial or die "Cannot create $partial: $!\n";
	print {$fh} JSON::PP->new->ascii->canonical->pretty->encode($data)
		or die "Cannot write $partial: $!\n";
	close $fh or die "Cannot close $partial: $!\n";
	rename $partial, $file or die "Cannot publish $file: $!\n";
}

sub intersection_count {
	my ($left_file, $right_file) = @_;
	open my $left, '<', $left_file or die "Cannot read $left_file: $!\n";
	open my $right, '<', $right_file or die "Cannot read $right_file: $!\n";
	my ($a, $b, $count) = (scalar(<$left>), scalar(<$right>), 0);
	while (defined($a) && defined($b)) {
		chomp($a, $b);
		if ($a eq $b) { $count++; $a = <$left>; $b = <$right>; }
		elsif ($a lt $b) { $a = <$left>; }
		else { $b = <$right>; }
	}
	close $left or die "Cannot close $left_file: $!\n";
	close $right or die "Cannot close $right_file: $!\n";
	return $count;
}

sub write_intersection_origin_rows {
	my ($candidate_file, $donor_file, $sink, $sample, $scope, $modbam) = @_;
	open my $candidate, '<', $candidate_file or die "Cannot read $candidate_file: $!\n";
	open my $donor, '<', $donor_file or die "Cannot read $donor_file: $!\n";
	my ($a, $b) = (scalar(<$candidate>), scalar(<$donor>));
	while (defined($a) && defined($b)) {
		chomp($a, $b);
		if ($a eq $b) {
			print {$sink} join("\t", $sample, $scope, $a, $modbam), "\n"
				or die "Cannot write compressed read-origin rows: $GzipError\n";
			$a = <$candidate>; $b = <$donor>;
		} elsif ($a lt $b) { $a = <$candidate>; }
		else { $b = <$donor>; }
	}
	close $candidate or die "Cannot close $candidate_file: $!\n";
	close $donor or die "Cannot close $donor_file: $!\n";
}

sub unit_paths {
	my ($state_dir, $target, $sample) = @_;
	my $base = File::Spec->catfile($state_dir, 'units', $target->{mode}, $target->{mgs}, file_component($sample));
	return ("$base.json", "$base.stone");
}

sub output_alignment_path {
	my ($out_dir, $target, $sample, $format) = @_;
	my $name = join('__', file_component($sample), $target->{mode},
		file_component($target->{representative_mag})) . ".mod.$format";
	return File::Spec->catfile($out_dir, $target->{mgs}, $name);
}

sub alignment_index_path {
	my ($alignment, $format) = @_;
	return $alignment . ($format eq 'cram' ? '.crai' : '.bai');
}

sub reference_lengths {
	my ($fasta) = @_;
	my $fh = IO::Uncompress::Gunzip->new($fasta, Transparent => 1, MultiStream => 1)
		or die "Cannot read representative FASTA $fasta: $GunzipError\n";
	my (%lengths, $name);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		if ($line =~ /^>(\S+)/) {
			$name = $1;
			die "Representative FASTA $fasta repeats contig '$name'\n" if exists $lengths{$name};
			$lengths{$name} = 0;
		} elsif ($line ne '') {
			die "Representative FASTA $fasta has sequence before a header\n" unless defined $name;
			$line =~ s/\s+//g;
			$lengths{$name} += length($line);
		}
	}
	close $fh or die "Cannot close representative FASTA $fasta: $!\n";
	die "Representative FASTA $fasta contains no nonempty contigs\n"
		unless keys(%lengths) && !grep { $_ <= 0 } values %lengths;
	return \%lengths;
}

sub bam_coverage {
	my ($samtools, $bam, $lengths) = @_;
	my $total_bases = 0;
	$total_bases += $_ for values %{$lengths};
	my ($covered_bases, $depth_sum) = (0, 0);
	my $command = tool_command($samtools, 'depth', '-d', 0, '-q', 0, '-Q', 0, $bam);
	my $fh = open_command("computing reference coverage for $bam", $command);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		my ($contig, $position, $depth) = split /\t/, $line, -1;
		die "Malformed samtools depth result for $bam: $line\n"
			unless defined($depth) && exists($lengths->{$contig})
			&& $position =~ /^\d+$/ && $position >= 1 && $position <= $lengths->{$contig}
			&& $depth =~ /^\d+$/;
		$covered_bases++ if $depth > 0;
		$depth_sum += $depth;
	}
	close_command($fh, "computing reference coverage for $bam");
	die "Coverage calculation for $bam reported $covered_bases positions across only $total_bases reference bases\n"
		if $covered_bases > $total_bases;
	return {
		reference_bases => $total_bases, covered_bases => $covered_bases,
		breadth_fraction => $total_bases ? $covered_bases / $total_bases : 0,
		mean_depth => $total_bases ? $depth_sum / $total_bases : 0,
		depth_sum => $depth_sum,
		method => 'samtools depth -d 0 -q 0 -Q 0; reference denominator includes every representative contig',
	};
}

sub read_filter_report {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read alignment-filter report $file: $!\n";
	my %patterns = (
		input_records => qr/^Input records:\s*(\d+)/,
		retained_mapped => qr/^Retained mapped records:\s*(\d+)/,
		newly_filtered => qr/^Newly filtered records:\s*(\d+)/,
		malformed => qr/^Malformed SAM records skipped:\s*(\d+)/,
		mapq_rejected => qr/^\s*Mapping quality \([^)]*\):\s*(\d+)/,
		coverage_rejected => qr/^\s*Query coverage \([^)]*\):\s*(\d+)/,
		edit_rate_rejected => qr/^\s*Edit rate \([^)]*\):\s*(\d+)/,
		end_clip_rejected => qr/^\s*Both ends clipped \([^)]*\):\s*(\d+)/,
	);
	my (%values, @warnings);
	while (my $line = <$fh>) {
		for my $key (keys %patterns) {
			$values{$key} = 0 + $1 if $line =~ $patterns{$key};
		}
		push @warnings, $line if $line =~ /\b(?:warning|error|malformed)\b/i
			&& $line !~ /^Malformed SAM records skipped:/;
	}
	close $fh or die "Cannot close $file: $!\n";
	for my $required (qw(input_records retained_mapped newly_filtered malformed)) {
		die "Alignment-filter report $file lacks '$required' statistics\n"
			unless exists $values{$required};
	}
	die "Alignment filter skipped $values{malformed} malformed SAM record(s); refusing an incomplete output. "
		. join(' | ', @warnings[0 .. ($#warnings < 4 ? $#warnings : 4)]) . "\n"
		if $values{malformed};
	return (\%values, \@warnings);
}

sub read_transfer_stats {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read tag-transfer statistics $file: $!\n";
	my %values;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		my ($key, $value, @extra) = split /\t/, $line, -1;
		die "Malformed tag-transfer statistics in $file\n"
			unless defined($key) && $key =~ /\A[a-z_]+\z/
			&& defined($value) && $value =~ /\A\d+\z/ && !@extra
			&& !exists($values{$key});
		$values{$key} = 0 + $value;
	}
	close $fh or die "Cannot close $file: $!\n";
	for my $required (qw(input_alignments output_alignments output_read_names
		dropped_no_primary_groups dropped_no_primary_alignments
		dropped_supplementary_alignments stripped_sa_tags legacy_missing_mn_reads
		legacy_tag_names_normalized implicit_mode_groups_normalized)) {
		die "Tag-transfer statistics $file lacks '$required'\n"
			unless exists $values{$required};
	}
	return \%values;
}

sub materialize_cram_reference {
	my ($source, $workdir, $samtools, $cache) = @_;
	return $cache->{$source} if exists $cache->{$source};
	my $destination = File::Spec->catfile($workdir,
		'reference.' . sha1_hex($source) . '.fa');
	my $input = IO::Uncompress::Gunzip->new($source, Transparent => 1, MultiStream => 1)
		or die "Cannot read representative FASTA $source for CRAM encoding: $GunzipError\n";
	open my $output, '>', $destination or die "Cannot create $destination: $!\n";
	while (my $line = <$input>) {
		print {$output} $line or die "Cannot write $destination: $!\n";
	}
	close $input or die "Cannot finish reading representative FASTA $source: $GunzipError\n";
	close $output or die "Cannot close $destination: $!\n";
	run_shell("indexing temporary CRAM reference", tool_command($samtools, 'faidx', $destination));
	$cache->{$source} = $destination;
	return $destination;
}

sub recorded_source_paths {
	my ($sample, $record) = @_;
	my $per_sample = File::Spec->catfile($record->{wrdir}, 'input_raw.txt');
	my @paths;
	my $source = '';
	if (-s $per_sample) {
		open my $fh, '<', $per_sample or die "Cannot read $per_sample: $!\n";
		local $/;
		@paths = split /;/, (<$fh> // '');
		close $fh or die "Cannot close $per_sample: $!\n";
		$source = $per_sample;
	} else {
		my $cohort = File::Spec->catfile(dirname($record->{wrdir}), 'Input_raw.txt');
		if (-s $cohort) {
			open my $fh, '<', $cohort or die "Cannot read $cohort: $!\n";
			while (my $line = <$fh>) {
				my ($row_sample, $value) = split /\t/, $line, 2;
				next unless defined($value) && $row_sample eq $sample;
				@paths = split /;/, $value;
				$source = $cohort;
				last;
			}
			close $fh or die "Cannot close $cohort: $!\n";
		}
	}
	for (@paths) { s/^\s+|\s+$//g; }
	@paths = grep { $_ ne '' } @paths;
	return ($source, \@paths);
}

sub write_plan {
	my ($file, $plan) = @_;
	my $partial = "$file.part.$$";
	push @publication_partials, $partial;
	open my $fh, '>', $partial or die "Cannot create $partial: $!\n";
	print {$fh} join("\t", qw(mgs mode representative_mag representative_fasta source_mag assembly_group assembly bin_assignment source_contigs eligible_scopes eligible_sample_scopes status)), "\n";
	for my $key (sort keys %{$plan->{targets}}) {
		my $target = $plan->{targets}{$key};
		next unless $target->{available};
		for my $source (@{$target->{sources}}) {
			my @source_scope_keys = sort grep {
				$plan->{scopes}{$_}{group} eq $source->{group}
			} keys %{$target->{scope_keys}};
			my @source_scopes = map {
				"$plan->{scopes}{$_}{sample}:$plan->{scopes}{$_}{scope}"
			} @source_scope_keys;
			print {$fh} join("\t",
				$target->{mgs}, $target->{mode}, $target->{representative_mag}, $target->{reference},
				$source->{mag}, $source->{group}, $source->{assembly}, $source->{assignment},
				scalar(@{$source->{contigs}}), scalar(@source_scopes), join(',', @source_scopes), 'ready',
			), "\n";
		}
	}
	for my $row (@{$plan->{unavailable}}) {
		print {$fh} join("\t", $row->{mgs}, $row->{mode}, '-', '-', '-', '-', '-', '-', 0, 0, '-', $row->{reason}), "\n";
	}
	close $fh or die "Cannot close $partial: $!\n";
	rename $partial, $file or die "Cannot publish $file: $!\n";
}

sub validate_cram_reference {
	my ($scope) = @_;
	(my $stat_file = $scope->{cram}) =~ s/\.cram\z/.reference.stat/;
	return unless -s $stat_file;
	open my $fh, '<', $stat_file or die "Cannot read $stat_file: $!\n";
	my $line = <$fh> // '';
	close $fh or die "Cannot close $stat_file: $!\n";
	$line =~ s/[\r\n]+$//;
	my ($recorded_size, $recorded_mtime) = split /\s+/, $line;
	my @stat = stat($scope->{assembly_reference});
	die "Assembly reference fingerprint in $stat_file does not match $scope->{assembly_reference}; the CRAM may target another assembly generation\n"
		unless defined($recorded_size) && defined($recorded_mtime)
		&& $recorded_size =~ /^\d+$/ && $recorded_mtime =~ /^\d+$/
		&& $recorded_size == $stat[7] && $recorded_mtime == $stat[9];
}

sub scan_scope_candidates {
	my (%options) = @_;
	my $scope = $options{scope};
	my $plan = $options{plan};
	my $assembly_targets = $plan->{contig_targets}{$scope->{assembly}} || {};
	my %target_files;
	for my $target_key (sort keys %{$plan->{targets}}) {
		my $target = $plan->{targets}{$target_key};
		next unless $target->{available} && $target->{scope_keys}{$scope->{key}}
			&& $options{active_target_scopes}{$target_key}{$scope->{key}};
		my $hash = sha1_hex($scope->{key} . "\0" . $target_key);
		$target_files{$target_key} = "$options{workdir}/candidate.$hash.raw";
	}
	validate_cram_reference($scope);
	run_shell("checking CRAM $scope->{cram}", tool_command($options{samtools}, 'quickcheck', $scope->{cram}));
	my $command = tool_command(
		$options{samtools}, 'view', '-@', $options{threads},
		'-T', $scope->{assembly_reference}, '-F', 3844,
		'-q', $options{source_min_mapq}, $scope->{cram},
	);
	my $fh = open_command("scanning candidate CRAM $scope->{cram}", $command);
	open my $identity, '>', $options{source_identity_raw}
		or die "Cannot create $options{source_identity_raw}: $!\n";
	my (%handles, %last_used, $clock);
	while (my $line = <$fh>) {
		my @fields = split /\t/, $line, 12;
		die "Malformed candidate CRAM alignment in $scope->{cram}\n" unless @fields >= 11;
		my ($name, $contig) = @fields[0, 2];
		next if $fields[4] == 255;
		next if cigar_query_coverage($fields[5], "$name in $scope->{cram}") < $options{source_min_coverage};
		next unless exists $assembly_targets->{$contig};
		my @matching_targets = grep { exists $target_files{$_} }
			keys %{$assembly_targets->{$contig}};
		next unless @matching_targets;
		my $source_native = native_sequence(\@fields, "Assembly CRAM $scope->{cram}");
		print {$identity} "$name\t", sha256_hex($source_native), "\n"
			or die "Cannot write $options{source_identity_raw}: $!\n";
		for my $target_key (@matching_targets) {
			$clock++;
			if (!exists $handles{$target_key}) {
				if (keys(%handles) >= 64) {
					my ($oldest) = sort { $last_used{$a} <=> $last_used{$b} } keys %handles;
					close $handles{$oldest} or die "Cannot close candidate-name partition: $!\n";
					delete $handles{$oldest};
					delete $last_used{$oldest};
				}
				open my $out, '>>', $target_files{$target_key}
					or die "Cannot append $target_files{$target_key}: $!\n";
				$handles{$target_key} = $out;
			}
			$last_used{$target_key} = $clock;
			print {$handles{$target_key}} "$name\n"
				or die "Cannot write $target_files{$target_key}: $!\n";
		}
	}
	close_command($fh, "scanning candidate CRAM $scope->{cram}");
	close $identity or die "Cannot close $options{source_identity_raw}: $!\n";
	close $handles{$_} or die "Cannot close candidate-name partition: $!\n" for keys %handles;
	my %result;
	for my $target_key (sort keys %target_files) {
		my $raw = $target_files{$target_key};
		if (!-e $raw) {
			open my $empty, '>', $raw or die "Cannot create $raw: $!\n";
			close $empty or die "Cannot close $raw: $!\n";
		}
		my $sorted = "$raw.names";
		sort_unique($raw, $sorted, $options{workdir});
		$result{$target_key} = { file => $sorted, count => count_lines($sorted) };
	}
	return \%result;
}

sub align_and_transfer {
	my (%options) = @_;
	my $started = Time::HiRes::time();
	my $hash = sha1_hex(join("\0", $options{scope}{key}, $options{target}{key}));
	my $prefix = "$options{workdir}/align.$hash";
	my $donor_target = "$prefix.donor.name.bam";
	filter_bam_by_names(
		samtools => $options{samtools}, names => $options{names},
		input => $options{donor_union}, output => $donor_target,
		workdir => $options{workdir}, exclude_flags => 0,
	);
	my $donor_count = bam_count($options{samtools}, $donor_target);
	die "Target donor extraction returned $donor_count records for $options{candidate_count} candidate names\n"
		unless $donor_count == $options{candidate_count};

	my $acceptor_unsorted = "$prefix.acceptor.unsorted.bam";
	my $preset = $options{preset};
	my @filter = $options{scope}{technology} eq 'ONT'
		? @{$options{filter_ont}} : @{$options{filter_pb}};
	my $rg_id = file_component("$options{scope}{sample}.$options{scope}{scope}");
	my $rg = "\@RG\\tID:$rg_id\\tSM:$options{scope}{sample}\\tPL:"
		. ($options{scope}{technology} eq 'ONT' ? 'ONT' : 'PACBIO')
		. "\\tLB:meth2rep.$options{scope}{scope}";
	my $filter_report = "$prefix.bamFilter.log";
	my $mapper_report = "$prefix.minimap2.log";
	my $mapping_command = tool_command(
		$options{samtools}, 'fastq', '-n', $donor_target,
	) . ' | ' . tool_command(
		$options{minimap2}, '-2', '-a', '-Y', '-t', $options{threads}, '--secondary=no',
		'-x', $preset, '-R', $rg, $options{reference_index}, '-',
	) . ' 2> ' . shell_quote($mapper_report)
		. ' | ' . tool_command($options{bam_filter}, @filter)
		. ' 2> ' . shell_quote($filter_report)
		. ' | ' . tool_command(
			$options{samtools}, 'view', '-b', '-F', 4,
			'-o', $acceptor_unsorted, '-',
		);
	my $mapping_ok = eval {
		run_shell("fresh representative mapping for $options{target}{mgs}", $mapping_command);
		1;
	};
	if (!$mapping_ok) {
		my $reason = $@ || "representative mapping failed\n";
		for my $log ($mapper_report, $filter_report) {
			next unless -s $log;
			open my $fh, '<', $log or next;
			my @lines = <$fh>;
			close $fh;
			@lines = @lines[-10 .. -1] if @lines > 10;
			$reason .= "\n$log:\n" . join('', @lines);
		}
		die $reason;
	}
	my ($filter_stats, $filter_warnings) = read_filter_report($filter_report);
	open my $mapper_fh, '<', $mapper_report or die "Cannot read $mapper_report: $!\n";
	my @mapper_warnings = grep { /\b(?:warning|error)\b/i } <$mapper_fh>;
	close $mapper_fh or die "Cannot close $mapper_report: $!\n";
	my @warnings = (@{$filter_warnings}, @mapper_warnings);
	s/\s+\z// for @warnings;
	my $accepted_mapping_count = bam_count($options{samtools}, $acceptor_unsorted);
	if (!$accepted_mapping_count) {
		remove_temp_files($donor_target, $acceptor_unsorted, $filter_report, $mapper_report);
		return {
			aligned => 0, donor => $donor_count, filter_stats => $filter_stats,
			warnings => \@warnings,
			elapsed_seconds => Time::HiRes::time() - $started,
		};
	}

	my $acceptor_name = "$prefix.acceptor.name.bam";
	run_shell(
		"name-sorting fresh representative mappings",
		tool_command($options{samtools}, 'sort', '-n', '-@', $options{threads},
			'-m', $options{sort_memory}, '-o', $acceptor_name, $acceptor_unsorted),
	);
	my $transferred_name = "$prefix.transferred.name.bam";
	my $transfer_stats_file = "$prefix.transfer.tsv";
	my @transfer_options = (
		'--samtools', $options{samtools}, '--donor', $donor_target,
		'--acceptor', $acceptor_name, '--output', $transferred_name,
		'--stats', $transfer_stats_file,
		'--supplementary', $options{supplementary_alignments},
	);
	push @transfer_options, '--allow-missing-mn' if $options{allow_missing_mn};
	run_shell(
		"transferring MM/ML tags for full-sequence representative alignments",
		tool_command($^X, "$Bin/transfer_mod_tags.pl", @transfer_options),
	);
	my $transfer_stats = read_transfer_stats($transfer_stats_file);
	die "Tag-transfer input count changed from $accepted_mapping_count to $transfer_stats->{input_alignments} for $options{scope}{sample} $options{target}{mgs}\n"
		unless $transfer_stats->{input_alignments} == $accepted_mapping_count;
	my $transferred_count = bam_count($options{samtools}, $transferred_name);
	die "Tag-transfer BAM contains $transferred_count records but reports $transfer_stats->{output_alignments} for $options{scope}{sample} $options{target}{mgs}\n"
		unless $transferred_count == $transfer_stats->{output_alignments};
	push @warnings, "$transfer_stats->{dropped_no_primary_alignments} orphan supplementary alignment(s) were dropped after their primary failed filtering"
		if $transfer_stats->{dropped_no_primary_alignments};
	push @warnings, "$transfer_stats->{dropped_supplementary_alignments} supplementary alignment(s) were removed by the '$options{supplementary_alignments}' output policy"
		if $transfer_stats->{dropped_supplementary_alignments};
	push @warnings, "$transfer_stats->{legacy_missing_mn_reads} donor read(s) lacked MN and were accepted under --allow-missing-mn"
		if $transfer_stats->{legacy_missing_mn_reads};
	push @warnings, "$transfer_stats->{legacy_tag_names_normalized} legacy Mm/Ml tag name(s) were normalized to standard MM/ML"
		if $transfer_stats->{legacy_tag_names_normalized};
	push @warnings, "$transfer_stats->{implicit_mode_groups_normalized} MM group(s) without an explicit mode were normalized to the SAM-equivalent '.' mode for default modkit compatibility"
		if $transfer_stats->{implicit_mode_groups_normalized};
	remove_temp_files($donor_target, $acceptor_unsorted, $acceptor_name,
		$filter_report, $mapper_report, $transfer_stats_file);
	if (!$transferred_count) {
		remove_temp_files($transferred_name);
		return {
			donor => $donor_count, aligned => 0,
			filter_stats => $filter_stats, transfer_stats => $transfer_stats,
			warnings => \@warnings,
			elapsed_seconds => Time::HiRes::time() - $started,
		};
	}
	return {
		donor => $donor_count, aligned => $transferred_count,
		transferred_name => $transferred_name,
		filter_stats => $filter_stats, transfer_stats => $transfer_stats,
		warnings => \@warnings,
		elapsed_seconds => Time::HiRes::time() - $started,
	};
}

my ($mgs_dir, $map_file, $mgs_report, $representatives_dir, $binner, $manifest_file, $out_dir, $out_manifest);
my ($mgs_values, $mgs_file) = ('', '');
my ($mgs2rep, $rep2rep) = (0, 0);
my ($source_min_mapq, $source_min_coverage, $threads, $memory_gb) = (10, 0.5, 4, 32);
my ($target_min_mapq, $target_min_coverage, $target_max_edit_rate, $target_min_end_clip);
my ($filter_ont, $filter_pb) = ('0.15 0.5 10 0', '0.05 0.5 30 0');
my ($preset_ont, $preset_pb) = ('map-ont', 'map-pb');
my ($supplementary_alignments, $allow_missing_mn, $output_format) = ('drop', 0, 'bam');
my ($tmp_parent, $keep_read_ids, $redo, $plan_only, $help) = ('', 0, 0, 0, 0);
my ($samtools, $minimap2, $bam_filter) = ('', '', '');

GetOptions(
	'mgs-dir=s' => \$mgs_dir,
	'map=s' => \$map_file,
	'mgs-report=s' => \$mgs_report,
	'representatives-dir=s' => \$representatives_dir,
	'binner=s' => \$binner,
	'modbam-manifest=s' => \$manifest_file,
	'out|o=s' => \$out_dir,
	'out-manifest=s' => \$out_manifest,
	'mgs=s' => \$mgs_values,
	'mgs-file=s' => \$mgs_file,
	'mgs2rep!' => \$mgs2rep,
	'rep2rep!' => \$rep2rep,
	'source-min-mapq=i' => \$source_min_mapq,
	'source-min-coverage=f' => \$source_min_coverage,
	'target-min-mapq=i' => \$target_min_mapq,
	'target-min-coverage=f' => \$target_min_coverage,
	'target-max-edit-rate=f' => \$target_max_edit_rate,
	'target-min-end-clip=i' => \$target_min_end_clip,
	'mapper-filter-ont=s' => \$filter_ont,
	'mapper-filter-pb=s' => \$filter_pb,
	'minimap2-preset-ont|mapper-preset-ont=s' => \$preset_ont,
	'minimap2-preset-pb|mapper-preset-pb=s' => \$preset_pb,
	'supplementary-alignments=s' => \$supplementary_alignments,
	'allow-missing-mn!' => \$allow_missing_mn,
	'output-format=s' => \$output_format,
	'threads=i' => \$threads,
	'memory-gb=i' => \$memory_gb,
	'tmp=s' => \$tmp_parent,
	'keep-read-ids!' => \$keep_read_ids,
	'redo!' => \$redo,
	'override!' => \$redo,
	'plan-only!' => \$plan_only,
	'samtools=s' => \$samtools,
	'minimap2=s' => \$minimap2,
	'bam-filter=s' => \$bam_filter,
	'help|h' => \$help,
) or die usage();
if ($help) {
	print usage();
	exit 0;
}
die usage() if @ARGV;
die "--mgs-dir is required\n" unless defined($mgs_dir) && $mgs_dir ne '';
$mgs_dir = File::Spec->rel2abs($mgs_dir);
$mgs_dir =~ s{/\z}{} unless $mgs_dir eq '/';
$mgs_dir = abs_path($mgs_dir) || $mgs_dir;
my $mgs_progress = require_mgs_readiness($mgs_dir);
$mgs_report ||= File::Spec->catfile($mgs_dir, 'MAGvsGC.txt.gz');
$representatives_dir ||= File::Spec->catdir($mgs_dir, 'Genomes', 'MGS_ctg');
if (-e $mgs_report) {
	my $stage1_stone = File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone');
	die "MGS membership report was modified after the Stage1 completion checkpoint: $mgs_report\n"
		if (stat($mgs_report))[9] > (stat($stage1_stone))[9];
}
if (!defined($binner) || $binner eq '') {
	($binner) = $mgs_dir =~ m{(?:^|/)Bin_([A-Za-z0-9_.-]+)/?\z};
}
$out_dir ||= File::Spec->catdir($mgs_dir, 'Meth2Rep');
for my $required (
	['--map', $map_file], ['--mgs-report', $mgs_report],
	['--representatives-dir', $representatives_dir], ['--binner', $binner],
	['--modbam-manifest', $manifest_file], ['--out', $out_dir],
) {
	die "$required->[0] is required\n" unless defined($required->[1]) && $required->[1] ne '';
}
die "Enable --mgs2rep, --rep2rep, or both\n" unless $mgs2rep || $rep2rep;
die "--threads must be a positive integer\n" unless $threads > 0;
die "--memory-gb must be a positive integer\n" unless $memory_gb > 0;
die "--source-min-mapq must be between 0 and 254 (255 means unknown)\n"
	unless $source_min_mapq >= 0 && $source_min_mapq <= 254;
die "--source-min-coverage must be between 0 and 1\n"
	unless $source_min_coverage >= 0 && $source_min_coverage <= 1;
die "--output-format must be 'bam' or 'cram'\n"
	unless $output_format eq 'bam' || $output_format eq 'cram';
die "--supplementary-alignments must be 'drop' or 'keep'\n"
	unless $supplementary_alignments eq 'drop' || $supplementary_alignments eq 'keep';
for my $preset_spec (['--minimap2-preset-ont', $preset_ont], ['--minimap2-preset-pb', $preset_pb]) {
	die "$preset_spec->[0] contains unsafe or unsupported characters\n"
		unless $preset_spec->[1] =~ /\A[A-Za-z0-9_.:+-]+\z/;
}
my @filter_ont = validate_filter($filter_ont, '--mapper-filter-ont');
my @filter_pb = validate_filter($filter_pb, '--mapper-filter-pb');
if (defined $target_min_mapq) {
	die "--target-min-mapq must be between 0 and 254\n"
		unless $target_min_mapq >= 0 && $target_min_mapq <= 254;
	$filter_ont[2] = $filter_pb[2] = $target_min_mapq;
}
if (defined $target_min_coverage) {
	die "--target-min-coverage must be between 0 and 1\n"
		unless $target_min_coverage >= 0 && $target_min_coverage <= 1;
	$filter_ont[1] = $filter_pb[1] = $target_min_coverage;
}
if (defined $target_max_edit_rate) {
	die "--target-max-edit-rate must be between 0 and 1\n"
		unless $target_max_edit_rate >= 0 && $target_max_edit_rate <= 1;
	$filter_ont[0] = $filter_pb[0] = $target_max_edit_rate;
}
if (defined $target_min_end_clip) {
	die "--target-min-end-clip must be a nonnegative integer\n"
		unless $target_min_end_clip >= 0;
	$filter_ont[3] = $filter_pb[3] = $target_min_end_clip;
}
my $sort_memory_mb = int(($memory_gb * 1024 * 0.5) / $threads);
$sort_memory_mb = 1 if $sort_memory_mb < 1;
$sort_memory_mb = 1024 if $sort_memory_mb > 1024;
my $sort_memory = $sort_memory_mb . 'M';

my $selected = read_target_mgs(mgs => $mgs_values, mgs_file => $mgs_file);
my $manifest = read_modbam_manifest($manifest_file);
my ($assembly_groups, $map) = getDirsPerAssmblGrp($map_file);
my @modes = grep { $_->[1] } (['mgs2rep', $mgs2rep], ['rep2rep', $rep2rep]);
@modes = map { $_->[0] } @modes;
my $plan = build_meth2rep_plan(
	map => $map, assembly_groups => $assembly_groups, manifest => $manifest,
	manifest_file => File::Spec->rel2abs($manifest_file), selected => $selected,
	mgs_report => File::Spec->rel2abs($mgs_report),
	representatives_dir => File::Spec->rel2abs($representatives_dir),
	binner => $binner, modes => \@modes,
);
my %checkpoint_representatives = map {
	my $path = $_->{path} // '';
	my $canonical = $path ne '' ? abs_path($path) : undef;
	defined($canonical) ? ($canonical => 1) : ()
}
	@{$mgs_progress->{BinExtr}{outputs}};
for my $target (values %{$plan->{targets}}) {
	next unless $target->{available};
	die "Selected representative $target->{reference} is not recorded by the completed MGS bin-extraction checkpoint\n"
		unless $checkpoint_representatives{$target->{reference}};
}
my %encoded_samples;
for my $scope (values %{$plan->{scopes}}) {
	my $encoded = file_component($scope->{sample});
	die "Sample names '$encoded_samples{$encoded}' and '$scope->{sample}' collide in meth2rep output filenames\n"
		if exists($encoded_samples{$encoded}) && $encoded_samples{$encoded} ne $scope->{sample};
	$encoded_samples{$encoded} = $scope->{sample};
}
my %plan_inputs = map { $_ => 1 } @{$plan->{inputs}};
$plan_inputs{File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone')} = 1;
$plan_inputs{File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone')} = 1;
for my $map_path (split /,/, $map_file) {
	$map_path = File::Spec->rel2abs($map_path);
	die "Resolved mapping input is missing: $map_path\n" unless -e $map_path;
	$plan_inputs{$map_path} = 1;
}
if ($mgs_file ne '') {
	my $selection_path = File::Spec->rel2abs($mgs_file);
	$plan_inputs{$selection_path} = 1;
}
$plan->{inputs} = [sort keys %plan_inputs];

$out_dir = File::Spec->rel2abs($out_dir);
$out_dir =~ s{/\z}{} unless $out_dir eq '/';
$out_dir = abs_path($out_dir) if -d $out_dir;
die "--out cannot be the MGS directory or one of its ancestors\n"
	if index("$mgs_dir/", "$out_dir/") == 0;
for my $input (@{$plan->{inputs}}) {
	die "meth2rep output directory contains a required input ($input); choose an isolated --out\n"
		if index($input, "$out_dir/") == 0;
}
make_path($out_dir);
$out_dir = abs_path($out_dir) || $out_dir;
die "--out resolves to the MGS directory or one of its ancestors\n"
	if index("$mgs_dir/", "$out_dir/") == 0;
for my $input (@{$plan->{inputs}}) {
	my $canonical_input = abs_path($input) || $input;
	die "meth2rep output directory contains a required input ($canonical_input); choose an isolated --out\n"
		if index($canonical_input, "$out_dir/") == 0;
}
my $state_dir = File::Spec->catdir($out_dir, '.meth2rep');
die "Legacy meth2rep mode-first directories are present in $out_dir; choose a fresh -o directory to avoid mixing output layouts\n"
	if -d File::Spec->catdir($out_dir, 'mgs2rep')
	|| -d File::Spec->catdir($out_dir, 'rep2rep');
make_path($state_dir);
$out_manifest ||= File::Spec->catfile($out_dir, 'manifest.tsv');
$out_manifest = File::Spec->rel2abs($out_manifest);
die "--out-manifest collides with a required input file: $out_manifest\n"
	if grep { $_ eq $out_manifest } @{$plan->{inputs}};
die "--out-manifest collides with a meth2rep control file: $out_manifest\n"
	if grep { $out_manifest eq File::Spec->catfile($state_dir, $_) }
		qw(plan.tsv plan.preview.tsv summary.tsv provenance.json complete.stone .meth2rep.lock);
die "--out-manifest cannot be placed inside a meth2rep data subdirectory: $out_manifest\n"
	if grep { index($out_manifest, File::Spec->catdir($out_dir, $_) . '/') == 0 }
		qw(.meth2rep read_ids);
die "--out-manifest cannot be placed inside a per-MGS output directory or another nested meth2rep data directory: $out_manifest\n"
	if index($out_manifest, "$out_dir/") == 0
	&& dirname($out_manifest) ne $out_dir;
die "--out-manifest names an existing directory rather than a TSV file: $out_manifest\n"
	if -d $out_manifest;
make_path(dirname($out_manifest)) unless -d dirname($out_manifest);
make_path(File::Spec->catdir($state_dir, 'units'));
my $lock_path = File::Spec->catfile($state_dir, '.meth2rep.lock');
open my $lock_fh, '>>', $lock_path or die "Cannot open meth2rep lock $lock_path: $!\n";
flock($lock_fh, LOCK_EX | LOCK_NB) or die "Another meth2rep run is using $out_dir\n";
my $plan_file = File::Spec->catfile($state_dir, $plan_only ? 'plan.preview.tsv' : 'plan.tsv');
my $summary_file = File::Spec->catfile($state_dir, 'summary.tsv');
my $provenance_file = File::Spec->catfile($state_dir, 'provenance.json');
my $checkpoint_file = File::Spec->catfile($state_dir, 'complete.stone');
my (%tool_versions, %tool_identities);
my $has_runnable_target = grep { $_->{available} } values %{$plan->{targets}};
if ($has_runnable_target) {
	$samtools ||= find_on_path('samtools')
		or die "samtools was not found on PATH; provide --samtools COMMAND\n";
	$minimap2 ||= find_on_path('minimap2')
		or die "minimap2 was not found on PATH; provide --minimap2 COMMAND\n";
	$bam_filter ||= join(' ', shell_quote($^X),
		shell_quote(File::Spec->catfile($Bin, '..', 'assemblies', 'bamFilter.pl')));
	$samtools = canonical_command($samtools, 'samtools');
	$minimap2 = canonical_command($minimap2, 'minimap2');
	$bam_filter = canonical_command($bam_filter, 'bamFilter');
	for my $tool_spec (
		['samtools', $samtools, '--version'], ['minimap2', $minimap2, '--version'],
	) {
		my $version = capture_command(
			"checking $tool_spec->[0]", tool_command($tool_spec->[1], $tool_spec->[2]),
		);
		($version) = split /\n/, $version, 2;
		$tool_versions{$tool_spec->[0]} = $version;
	}
	my %required_presets = map {
		my $technology = $_->{technology};
		($technology => ($technology eq 'ONT' ? $preset_ont : $preset_pb));
	} values %{$plan->{scopes}};
	for my $technology (sort keys %required_presets) {
		capture_command(
			"validating minimap2 preset $required_presets{$technology} for $technology",
			tool_command($minimap2, '-x', $required_presets{$technology}, '--version'),
		);
	}
	$tool_identities{samtools} = command_identity($samtools);
	$tool_identities{minimap2} = command_identity($minimap2);
	$tool_identities{bam_filter} = command_identity($bam_filter);
} else {
	($samtools, $minimap2, $bam_filter) = ('not_used_this_invocation') x 3;
	%tool_versions = map { $_ => 'not_used_this_invocation' } qw(samtools minimap2);
	%tool_identities = map { $_ => 'not_used_this_invocation' } qw(samtools minimap2 bam_filter);
}
my %checkpoint_parameters = (
	component_version => $VERSION,
	modes => join(',', @modes),
	mgs => join(',', sort keys %{$selected}),
	binner => $binner,
	source_min_mapq => $source_min_mapq,
	source_min_coverage => $source_min_coverage,
	mapper_filter_ont => join(' ', @filter_ont),
	mapper_filter_pb => join(' ', @filter_pb),
	minimap2_preset_ont => $preset_ont,
	minimap2_preset_pb => $preset_pb,
	supplementary_alignments => $supplementary_alignments,
	allow_missing_mn => $allow_missing_mn ? 1 : 0,
	output_format => $output_format,
	samtools_identity => $tool_identities{samtools},
	minimap2_identity => $tool_identities{minimap2},
	bam_filter_identity => $tool_identities{bam_filter},
	samtools_version => $tool_versions{samtools},
	minimap2_version => $tool_versions{minimap2},
	keep_read_ids => $keep_read_ids ? 1 : 0,
);
$checkpoint_parameters{input_fingerprint} = fingerprint_inputs($plan->{inputs}, \%checkpoint_parameters);

my (%active_target_scopes, %cached_units, %unit_context);
for my $target_key (sort keys %{$plan->{targets}}) {
	my $target = $plan->{targets}{$target_key};
	next unless $target->{available};
	my %sample_scopes;
	for my $scope_key (keys %{$target->{scope_keys}}) {
		my $sample = $plan->{scopes}{$scope_key}{sample};
		push @{$sample_scopes{$sample}}, $scope_key;
	}
	for my $sample (sort keys %sample_scopes) {
		my @scope_keys = sort @{$sample_scopes{$sample}};
		my %unit_groups = map { $plan->{scopes}{$_}{group} => 1 } @scope_keys;
		my %files = map { $_ => 1 } (
			File::Spec->rel2abs(__FILE__),
			File::Spec->catfile($Bin, 'transfer_mod_tags.pl'),
			File::Spec->catfile($Bin, '..', 'assemblies', 'bamFilter.pl'),
			File::Spec->catfile($Bin, '..', '..', 'Mods', 'Meth2Rep.pm'),
			File::Spec->catfile($Bin, '..', '..', 'Mods', 'Checkpoint.pm'),
			File::Spec->rel2abs($mgs_report),
			File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone'),
			File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone'),
			$target->{reference},
			(map { $_->{assignment} }
				grep { $unit_groups{$_->{group}} } @{$target->{sources}}),
			(map { File::Spec->rel2abs($_) } split /,/, $map_file),
		);
		for my $scope_key (@scope_keys) {
			my $scope = $plan->{scopes}{$scope_key};
			$files{$_} = 1 for ($scope->{cram}, "$scope->{cram}.sto",
				$scope->{assembly_reference}, @{$scope->{modbams}});
			(my $reference_stat = $scope->{cram}) =~ s/\.cram\z/.reference.stat/;
			$files{$reference_stat} = 1 if -e $reference_stat;
		}
		my %unit_parameters = (
			component_version => $VERSION, mode => $target->{mode},
			mgs => $target->{mgs}, sample => $sample,
			binner => $binner,
			source_min_mapq => $source_min_mapq,
			source_min_coverage => $source_min_coverage,
			mapper_filter_ont => join(' ', @filter_ont),
			mapper_filter_pb => join(' ', @filter_pb),
			minimap2_preset_ont => $preset_ont,
			minimap2_preset_pb => $preset_pb,
			supplementary_alignments => $supplementary_alignments,
			allow_missing_mn => $allow_missing_mn ? 1 : 0,
			output_format => $output_format,
			samtools_identity => $tool_identities{samtools},
			minimap2_identity => $tool_identities{minimap2},
			bam_filter_identity => $tool_identities{bam_filter},
			samtools_version => $tool_versions{samtools},
			minimap2_version => $tool_versions{minimap2},
			keep_read_ids => $keep_read_ids ? 1 : 0,
		);
		$unit_parameters{input_fingerprint} = fingerprint_inputs([sort keys %files], \%unit_parameters);
		my ($json, $stone) = unit_paths($state_dir, $target, $sample);
		my $unit_key = "$target_key\t$sample";
		my $expected_alignment = output_alignment_path($out_dir, $target, $sample, $output_format);
		my $expected_index = alignment_index_path($expected_alignment, $output_format);
		my $alternate_format = $output_format eq 'bam' ? 'cram' : 'bam';
		my $alternate_alignment = output_alignment_path($out_dir, $target, $sample, $alternate_format);
		my $alternate_index = alignment_index_path($alternate_alignment, $alternate_format);
		die "Untracked output alignment or index already exists: $expected_alignment; inspect it and use --override only if replacement is intended\n"
			if (-e $expected_alignment || -e $expected_index) && !-e $stone && !$redo;
		die "An alternate-format meth2rep output already exists ($alternate_alignment); use --override to replace it without leaving ambiguous BAM/CRAM siblings\n"
			if (-e $alternate_alignment || -e $alternate_index) && !$redo;
		$unit_context{$unit_key} = {
			json => $json, stone => $stone, parameters => \%unit_parameters,
			scope_keys => \@scope_keys,
			expected_alignment => $expected_alignment,
			expected_index => $expected_index,
			alternate_alignment => $alternate_alignment,
			alternate_index => $alternate_index,
		};
		if (!$redo && -s $stone && checkpoint_valid($stone, parameters => \%unit_parameters) && -s $json) {
			my $record = read_json_file($json);
			die "Unit checkpoint metadata mismatch: $json\n"
				unless ($record->{mgs} // '') eq $target->{mgs}
				&& ($record->{mode} // '') eq $target->{mode}
				&& ($record->{sample} // '') eq $sample;
			$cached_units{$unit_key} = $record;
			next;
		}
		$active_target_scopes{$target_key}{$_} = 1 for @scope_keys;
	}
}

if ($plan_only) {
	write_plan($plan_file, $plan);
	print "meth2rep v$VERSION: input, tool, cache, and output preflight passed; preview written to $plan_file\n";
	exit 0;
}
unlink $checkpoint_file or die "Cannot invalidate stale meth2rep checkpoint $checkpoint_file: $!\n"
	if -e $checkpoint_file;
write_plan($plan_file, $plan);

my $workdir = '';
if (keys %active_target_scopes) {
	$tmp_parent = $ENV{SLURM_TMPDIR}
		if $tmp_parent eq '' && defined($ENV{SLURM_TMPDIR}) && -d $ENV{SLURM_TMPDIR};
	$tmp_parent = File::Spec->rel2abs($tmp_parent) if $tmp_parent ne '';
	make_path($tmp_parent) if $tmp_parent ne '' && !-d $tmp_parent;
	my %temp_options = (CLEANUP => 1);
	$temp_options{DIR} = $tmp_parent if $tmp_parent ne '';
	$workdir = tempdir('meth2rep.XXXXXX', %temp_options);
}

print "meth2rep v$VERSION: " . scalar(keys %{$selected}) . " selected MGS; modes "
	. join(', ', @modes) . "; " . scalar(keys %{$plan->{scopes}}) . " eligible sample/scope CRAM(s)\n";
warn "Manifest row is not required by the selected MGS plan: $_\n"
	for @{$plan->{manifest_unused}};

my (%candidate, %scope_results, %reference_indexes, %nonempty_targets, %donor_scopes);
my (%source_identity, %cram_reference_cache);
my (%source_scan_seconds, %donor_stream_seconds);
for my $scope_key (sort keys %{$plan->{scopes}}) {
	my $scope = $plan->{scopes}{$scope_key};
	next unless grep { $active_target_scopes{$_}{$scope_key} } keys %active_target_scopes;
	my $scan_started = Time::HiRes::time();
	my $scope_hash = sha1_hex($scope_key);
	my $source_identity_raw = "$workdir/source.$scope_hash.identity.raw";
	$candidate{$scope_key} = scan_scope_candidates(
		scope => $scope, plan => $plan, workdir => $workdir,
		samtools => $samtools, threads => $threads, source_min_mapq => $source_min_mapq,
		source_min_coverage => $source_min_coverage,
		active_target_scopes => \%active_target_scopes,
		source_identity_raw => $source_identity_raw,
	);
	my $source_identity_sorted = "$source_identity_raw.sorted";
	sort_identity_records($source_identity_raw, $source_identity_sorted, $workdir);
	$source_identity{$scope_key} = $source_identity_sorted;
	remove_temp_files($source_identity_raw);
	$source_scan_seconds{$scope_key} = Time::HiRes::time() - $scan_started;
	my @nonempty_targets = grep { $candidate{$scope_key}{$_}{count} > 0 }
		sort keys %{$candidate{$scope_key}};
	next unless @nonempty_targets;
	$nonempty_targets{$scope_key} = \@nonempty_targets;
	push @{$donor_scopes{$_}}, $scope_key for @{$scope->{modbams}};
}

# Keep candidate membership scope-specific, but defer donor access until every
# CRAM has been scanned. This guarantees one record-streaming pass through each
# distinct physical modBAM even when a manifest explicitly reuses it.
my (%donor_name_bam, %scope_donor_names);
for my $modbam (sort keys %donor_scopes) {
	my $donor_started = Time::HiRes::time();
	my @scope_keys = sort @{$donor_scopes{$modbam}};
	my %technologies = map { $plan->{scopes}{$_}{technology} => 1 } @scope_keys;
	die "Original modBAM '$modbam' is declared for more than one sequencing technology ("
		. join(', ', sort keys %technologies) . "); split it into unambiguous manifest donors\n"
		if keys(%technologies) > 1;
	my $donor_hash = sha1_hex($modbam);
	my $union_raw = "$workdir/donor.$donor_hash.union.raw";
	concatenate_files(
		$union_raw,
		map {
			my $scope_key = $_;
			map { $candidate{$scope_key}{$_}{file} } @{$nonempty_targets{$scope_key}}
		} @scope_keys,
	);
	my $union_names = "$workdir/donor.$donor_hash.union.names";
	sort_unique($union_raw, $union_names, $workdir);
	# Raw, unmapped modBAMs legitimately have no @SQ records, which makes older
	# samtools quickcheck versions reject them.  Parsing the full selected stream
	# below catches truncation; this header read supplies the early format check.
	run_shell("checking donor modBAM header $modbam",
		tool_command($samtools, 'view', '-H', $modbam) . ' >/dev/null');
	my $donor_unsorted = "$workdir/donor.$donor_hash.primary.bam";
	filter_bam_by_names(
		samtools => $samtools, names => $union_names,
		input => $modbam, output => $donor_unsorted,
		workdir => $workdir, exclude_flags => 3840,
	);
	my $donor_name = "$workdir/donor.$donor_hash.name.bam";
	run_shell(
		"name-sorting donor subset $modbam",
		tool_command($samtools, 'sort', '-n', '-@', $threads,
			'-m', $sort_memory, '-o', $donor_name, $donor_unsorted),
	);
	my $found_raw = "$workdir/donor.$donor_hash.found.raw";
	my ($donor_count) = validate_donor_subset(
		$samtools, $donor_name, $found_raw, undef, $allow_missing_mn,
	);
	my $found_names = "$workdir/donor.$donor_hash.found.names";
	sort_unique($found_raw, $found_names, $workdir);
	$donor_name_bam{$modbam} = $donor_name;
	$donor_stream_seconds{$modbam} = Time::HiRes::time() - $donor_started;
	remove_temp_files($union_raw, $union_names, $donor_unsorted, $found_raw, $found_names);
}

for my $scope_key (sort keys %nonempty_targets) {
		my $scope = $plan->{scopes}{$scope_key};
		my $scope_hash = sha1_hex($scope_key);
		my $scope_union_raw = "$workdir/scope.$scope_hash.union.raw";
		concatenate_files($scope_union_raw,
			map { $candidate{$scope_key}{$_}{file} } @{$nonempty_targets{$scope_key}});
		my $scope_union_names = "$workdir/scope.$scope_hash.union.names";
		sort_unique($scope_union_raw, $scope_union_names, $workdir);
		my @scope_donor_bams;
		for my $modbam (@{$scope->{modbams}}) {
			my $donor_hash = sha1_hex("$scope_key\0$modbam");
			my $scope_donor_unsorted = "$workdir/scope.donor.$donor_hash.unsorted.bam";
			filter_bam_by_names(
				samtools => $samtools, names => $scope_union_names,
				input => $donor_name_bam{$modbam}, output => $scope_donor_unsorted,
				workdir => $workdir, exclude_flags => 0,
			);
			my $scope_donor_bam = "$workdir/scope.donor.$donor_hash.name.bam";
			run_shell("name-sorting scope donor subset",
				tool_command($samtools, 'sort', '-n', '-@', $threads,
					'-m', $sort_memory, '-o', $scope_donor_bam, $scope_donor_unsorted));
			my $raw_names = "$workdir/scope.donor.$donor_hash.names.raw";
			validate_donor_subset(
				$samtools, $scope_donor_bam, $raw_names, undef, $allow_missing_mn,
			);
			my $sorted_names = "$raw_names.sorted";
			sort_unique($raw_names, $sorted_names, $workdir);
			$scope_donor_names{$scope_key}{$modbam} = $sorted_names;
			push @scope_donor_bams, $scope_donor_bam;
			remove_temp_files($scope_donor_unsorted, $raw_names);
		}
		my $scope_donor_union = "$workdir/scope.$scope_hash.donors.name.bam";
		if (@scope_donor_bams == 1) {
			$scope_donor_union = $scope_donor_bams[0];
		} else {
			run_shell("merging selected original modBAM donors for $scope->{sample}:$scope->{scope}",
				tool_command($samtools, 'merge', '-n', '-@', $threads, '-f',
					$scope_donor_union, @scope_donor_bams));
		}
		my $combined_names_raw = "$workdir/scope.$scope_hash.found.raw";
		my $combined_identity_raw = "$workdir/scope.$scope_hash.identity.raw";
		my ($combined_count) = validate_donor_subset(
			$samtools, $scope_donor_union, $combined_names_raw,
			$combined_identity_raw, $allow_missing_mn,
		);
		my $combined_names = "$combined_names_raw.sorted";
		sort_unique($combined_names_raw, $combined_names, $workdir);
		my $combined_identity = "$combined_identity_raw.sorted";
		sort_identity_records($combined_identity_raw, $combined_identity, $workdir);
		assert_name_subset($scope_union_names, $combined_names,
			"sample '$scope->{sample}' scope '$scope->{scope}' across all declared modBAM donors");
		die "Donor identity collision in sample '$scope->{sample}' scope '$scope->{scope}'\n"
			unless $combined_count == count_lines($scope_union_names);
		assert_source_donor_identity(
			$source_identity{$scope_key}, $combined_identity,
			"sample '$scope->{sample}' scope '$scope->{scope}'",
		);
		remove_temp_files($scope_union_raw, $scope_union_names,
			$combined_names_raw, $combined_names,
			$combined_identity_raw, $combined_identity,
			$source_identity{$scope_key});
		my %mgs_for_scope;
		for my $target_key (@{$nonempty_targets{$scope_key}}) {
			my $target = $plan->{targets}{$target_key};
			$mgs_for_scope{$target->{mgs}}{$target->{mode}} = $target_key;
		}
		for my $mgs (sort keys %mgs_for_scope) {
			my $by_mode = $mgs_for_scope{$mgs};
			my $primary_key = exists($by_mode->{mgs2rep}) ? $by_mode->{mgs2rep} : $by_mode->{rep2rep};
			my $primary_target = $plan->{targets}{$primary_key};
			my $reference = $primary_target->{reference};
			my $preset = $scope->{technology} eq 'ONT' ? $preset_ont : $preset_pb;
			my $index_key = join("\0", $reference, $preset);
			if (!exists $reference_indexes{$index_key}) {
				my $index = "$workdir/reference." . sha1_hex($index_key) . '.mmi';
				run_shell(
					"indexing representative $mgs with preset $preset",
					tool_command($minimap2, '-x', $preset, '-t', $threads,
						'-d', $index, $reference),
				);
				$reference_indexes{$index_key} = $index;
			}
			my $primary_result = align_and_transfer(
				scope => $scope, target => $primary_target,
				names => $candidate{$scope_key}{$primary_key}{file},
				candidate_count => $candidate{$scope_key}{$primary_key}{count},
				donor_union => $scope_donor_union, reference_index => $reference_indexes{$index_key},
				samtools => $samtools, minimap2 => $minimap2,
				bam_filter => $bam_filter, filter_ont => \@filter_ont, filter_pb => \@filter_pb,
				preset => $preset, supplementary_alignments => $supplementary_alignments,
				allow_missing_mn => $allow_missing_mn,
				threads => $threads, sort_memory => $sort_memory, workdir => $workdir,
			);
			$scope_results{$primary_key}{$scope_key} = $primary_result;

			if (exists($by_mode->{rep2rep}) && $primary_key ne $by_mode->{rep2rep}) {
				my $rep_key = $by_mode->{rep2rep};
				my $rep_candidate = $candidate{$scope_key}{$rep_key};
				if (!$primary_result->{aligned}) {
					$scope_results{$rep_key}{$scope_key} = {
						donor => $rep_candidate->{count}, aligned => 0,
						reused_from => $primary_key, elapsed_seconds => 0,
					};
				} else {
					my $subset_started = Time::HiRes::time();
					my $subset = "$workdir/transferred." . sha1_hex("$scope_key\0$rep_key") . '.name.bam';
					filter_bam_by_names(
						samtools => $samtools, names => $rep_candidate->{file},
						input => $primary_result->{transferred_name}, output => $subset,
						workdir => $workdir, exclude_flags => 0,
					);
					my $subset_count = bam_count($samtools, $subset);
					$scope_results{$rep_key}{$scope_key} = {
						donor => $rep_candidate->{count}, aligned => $subset_count,
						($subset_count ? (transferred_name => $subset) : ()),
						reused_from => $primary_key,
						elapsed_seconds => Time::HiRes::time() - $subset_started,
					};
					remove_temp_files($subset) unless $subset_count;
				}
			}
			for my $completed_key (values %{$by_mode}) {
				my $result = $scope_results{$completed_key}{$scope_key};
				next unless $result && $result->{aligned};
				my $coordinate = "$workdir/coordinate."
					. sha1_hex("$scope_key\0$completed_key") . '.bam';
				run_shell(
					"coordinate-sorting transferred alignment",
					tool_command($samtools, 'sort', '-@', $threads,
						'-m', $sort_memory, '-o', $coordinate, $result->{transferred_name}),
				);
				$result->{coordinate_bam} = $coordinate;
				remove_temp_files($result->{transferred_name});
				delete $result->{transferred_name};
			}
		}
		my %scope_bams = map { $_ => 1 } (@scope_donor_bams, $scope_donor_union);
		remove_temp_files(sort keys %scope_bams);
}
remove_temp_files(values %donor_name_bam);

my @summary_rows;
push @summary_rows, {
	mgs => $_->{mgs}, mode => $_->{mode}, sample => '-', scopes => '-',
	status => $_->{reason}, candidates => 0, donor => 0, aligned => 0,
	mm => 0, ml => 0, alignment => '-', alignment_format => '-', index => '-',
} for @{$plan->{unavailable}};
my @checkpoint_outputs = ($plan_file);
my %reference_length_cache;
for my $target_key (sort keys %{$plan->{targets}}) {
	my $target = $plan->{targets}{$target_key};
	next unless $target->{available};
	my %sample_scopes;
	for my $scope_key (sort keys %{$target->{scope_keys}}) {
		my $scope = $plan->{scopes}{$scope_key};
		$sample_scopes{$scope->{sample}}{$scope->{scope}} = $scope_key;
	}
	for my $sample (sort keys %sample_scopes) {
		my $unit_key = "$target_key\t$sample";
		if (exists $cached_units{$unit_key}) {
			push @summary_rows, $cached_units{$unit_key}{summary};
			push @checkpoint_outputs, $unit_context{$unit_key}{json};
			my $cached_alignment = $cached_units{$unit_key}{summary}{alignment};
			my $cached_index = $cached_units{$unit_key}{summary}{index};
			push @checkpoint_outputs, $cached_alignment, $cached_index
				if defined($cached_alignment) && $cached_alignment ne '-';
			next;
		}
		my (@coordinate_bams, @scope_names, @name_files, @scope_reports, @unit_warnings);
		my ($candidate_count, $donor_count, $aligned_count, $processing_seconds) = (0, 0, 0, 0);
		my @unit_donors;
		for my $scope_name (sort keys %{$sample_scopes{$sample}}) {
			my $scope_key = $sample_scopes{$sample}{$scope_name};
			push @scope_names, $scope_name;
			my $candidate_record = $candidate{$scope_key}{$target_key};
			my $count = $candidate_record ? $candidate_record->{count} : 0;
			$candidate_count += $count;
			push @name_files, [$scope_name, $candidate_record->{file}] if $candidate_record && $count;
			my $scope_report = {
				scope => $scope_name,
				technology => $plan->{scopes}{$scope_key}{technology},
				assembly_cram => $plan->{scopes}{$scope_key}{cram},
				historical_source_filter_provenance => $HISTORICAL_SOURCE_FILTER,
				candidate_reads => $count,
				donor_reads => 0, accepted_alignment_records => 0,
				minimap2_preset => $plan->{scopes}{$scope_key}{technology} eq 'ONT'
					? $preset_ont : $preset_pb,
			};
			push @scope_reports, $scope_report;
			for my $modbam (@{$plan->{scopes}{$scope_key}{modbams}}) {
				my $donor_names = $scope_donor_names{$scope_key}{$modbam};
				my $donor_candidates = $count && $donor_names
					? intersection_count($candidate_record->{file}, $donor_names) : 0;
				push @unit_donors, {
					scope => $scope_name,
					technology => $plan->{scopes}{$scope_key}{technology},
					modbam => $modbam,
					candidate_reads => $donor_candidates,
				};
			}
			next unless $count;
			my $result = $scope_results{$target_key}{$scope_key}
				or die "Internal error: missing alignment result for $scope_key / $target_key\n";
			$scope_report->{donor_reads} = $result->{donor};
			$scope_report->{accepted_alignment_records} = $result->{aligned};
			$scope_report->{filter_stats} = $result->{filter_stats} if $result->{filter_stats};
			$scope_report->{transfer_stats} = $result->{transfer_stats} if $result->{transfer_stats};
			$scope_report->{reused_mapping_from} = 'mgs2rep' if $result->{reused_from};
			push @unit_warnings, @{$result->{warnings} || []};
			$processing_seconds += $result->{elapsed_seconds} || 0;
			$donor_count += $result->{donor};
			$aligned_count += $result->{aligned};
			next unless $result->{aligned};
			die "Internal error: missing coordinate-sorted transfer for $scope_key / $target_key\n"
				unless defined($result->{coordinate_bam}) && -s $result->{coordinate_bam};
			push @coordinate_bams, $result->{coordinate_bam};
		}
		my $attributed = 0;
		$attributed += $_->{candidate_reads} for @unit_donors;
		die "Donor attribution for $target->{mgs} $target->{mode} $sample accounts for $attributed of $candidate_count candidates\n"
			unless $attributed == $candidate_count;
		for my $left_index (0 .. $#name_files) {
			for my $right_index ($left_index + 1 .. $#name_files) {
				my $shared = intersection_count($name_files[$left_index][1], $name_files[$right_index][1]);
				die "Read identity collision for $target->{mgs} $target->{mode} $sample: $shared QNAME(s) appear in both $name_files[$left_index][0] and $name_files[$right_index][0] scopes\n"
					if $shared;
			}
		}

		my $target_dir = File::Spec->catdir($out_dir, $target->{mgs});
		my $unit = $unit_context{$unit_key};
		my $expected_alignment = $unit->{expected_alignment};
		my $expected_index = $unit->{expected_index};
		my ($status, $final_alignment, $final_index, $mm_count, $ml_count) =
			('no_candidates', '-', '-', 0, 0);
		my $lengths = $reference_length_cache{$target->{reference}}
			||= reference_lengths($target->{reference});
		my $reference_bases = 0;
		$reference_bases += $_ for values %{$lengths};
		my $coverage = {
			reference_bases => $reference_bases, covered_bases => 0,
			breadth_fraction => 0, mean_depth => 0, depth_sum => 0,
			method => 'no accepted representative alignments',
		};
		if ($candidate_count && !$aligned_count) {
			$status = 'no_target_alignment';
			push @unit_warnings, 'Candidate reads passed source filtering, but none passed representative-alignment filtering';
		} elsif ($aligned_count) {
			$status = 'complete';
			make_path($target_dir);
			$final_alignment = $expected_alignment;
			$final_index = $expected_index;
			my $partial_alignment = "$final_alignment.part.$$";
			my $partial_index = alignment_index_path($partial_alignment, $output_format);
			push @publication_partials, $partial_alignment, $partial_index;
			my $merged_bam = $output_format eq 'bam'
				? $partial_alignment
				: "$workdir/final." . sha1_hex("$target_key\0$sample") . '.bam';
			run_shell(
				"merging transferred sample scopes",
				tool_command($samtools, 'merge', '-@', $threads, '-f', $merged_bam, @coordinate_bams),
			);
			if ($output_format eq 'cram') {
				my $cram_reference = materialize_cram_reference(
					$target->{reference}, $workdir, $samtools, \%cram_reference_cache,
				);
				run_shell(
					"encoding self-contained final modCRAM",
					tool_command($samtools, 'view', '-@', $threads, '-C',
						'-T', $cram_reference, '--output-fmt-option', 'embed_ref=1',
						'-o', $partial_alignment, $merged_bam),
				);
			}
			run_shell("checking final mod$output_format",
				tool_command($samtools, 'quickcheck', $partial_alignment));
			my ($published_count, $published_mm, $published_ml) =
				bam_tag_counts($samtools, $partial_alignment);
			die "Published mod$output_format record count changed from $aligned_count to $published_count for $target->{mgs} $sample\n"
				unless $published_count == $aligned_count;
			($mm_count, $ml_count) = ($published_mm, $published_ml);
			run_shell("indexing final mod$output_format",
				tool_command($samtools, 'index', $partial_alignment, $partial_index));
			run_shell("validating final alignment index",
				tool_command($samtools, 'idxstats', $partial_alignment) . ' >/dev/null');
			$coverage = bam_coverage($samtools, $partial_alignment, $lengths);
			rename $partial_alignment, $final_alignment
				or die "Cannot publish $final_alignment: $!\n";
			rename $partial_index, $final_index
				or die "Cannot publish $final_index: $!\n";
			if ($redo) {
				unlink $unit->{alternate_alignment}
					or die "Cannot remove replaced alternate alignment $unit->{alternate_alignment}: $!\n"
					if -e $unit->{alternate_alignment};
				unlink $unit->{alternate_index}
					or die "Cannot remove replaced alternate index $unit->{alternate_index}: $!\n"
					if -e $unit->{alternate_index};
			}
			push @checkpoint_outputs, $final_alignment, $final_index;
		} else {
			push @unit_warnings, 'No source-MAG candidate reads passed the source filters';
			unlink $expected_alignment
				or die "Cannot remove stale derived alignment $expected_alignment: $!\n"
				if -e $expected_alignment;
			unlink $expected_index
				or die "Cannot remove stale derived index $expected_index: $!\n"
				if -e $expected_index;
			if ($redo) {
				unlink $unit->{alternate_alignment}
					or die "Cannot remove replaced alternate alignment $unit->{alternate_alignment}: $!\n"
					if -e $unit->{alternate_alignment};
				unlink $unit->{alternate_index}
					or die "Cannot remove replaced alternate index $unit->{alternate_index}: $!\n"
					if -e $unit->{alternate_index};
			}
		}
		remove_temp_files(@coordinate_bams);

		my $published_origin_file = '';
		if ($keep_read_ids && @name_files) {
			my $id_dir = File::Spec->catdir($state_dir, 'read_ids', $target->{mode}, $target->{mgs});
			make_path($id_dir);
			my $id_file = File::Spec->catfile($id_dir, file_component($sample) . '.read_origins.tsv.gz');
			my $id_partial = "$id_file.part.$$";
			push @publication_partials, $id_partial;
			my $gzip_fh = IO::Compress::Gzip->new($id_partial)
				or die "Cannot create $id_partial: $GzipError\n";
			print {$gzip_fh} "sample\tscope\tqname\toriginal_modbam\n"
				or die "Cannot write $id_partial: $GzipError\n";
			for my $spec (@name_files) {
				my $scope_key = $sample . "\t" . $spec->[0];
				for my $modbam (@{$plan->{scopes}{$scope_key}{modbams}}) {
					my $donor_names = $scope_donor_names{$scope_key}{$modbam};
					next unless $donor_names;
					write_intersection_origin_rows($spec->[1], $donor_names,
						$gzip_fh, $sample, $spec->[0], $modbam);
				}
			}
			close $gzip_fh or die "Cannot close $id_partial: $GzipError\n";
			rename $id_partial, $id_file or die "Cannot publish $id_file: $!\n";
			push @checkpoint_outputs, $id_file;
			$published_origin_file = $id_file;
		}
		my $row = {
			mgs => $target->{mgs}, mode => $target->{mode}, sample => $sample,
			scopes => join(',', @scope_names), status => $status,
			candidates => $candidate_count, donor => $donor_count,
			aligned => $aligned_count, mm => $mm_count, ml => $ml_count,
			alignment => $final_alignment,
			alignment_format => $final_alignment eq '-' ? '-' : $output_format,
			index => $final_index,
		};
		push @summary_rows, $row;
		make_path(dirname($unit->{json}));
		write_json_file($unit->{json}, {
			format => 'matafiler-meth2rep-unit-v2',
			mgs => $target->{mgs}, mode => $target->{mode}, sample => $sample,
			created_epoch => 0 + time,
			input_fingerprint => $unit->{parameters}{input_fingerprint},
			representative_mag => $target->{representative_mag},
			reference_fasta => $target->{reference},
			coverage => $coverage,
			processing => { alignment_transfer_wall_seconds => $processing_seconds },
			warnings => \@unit_warnings,
			scope_reports => \@scope_reports,
			tools => {
				samtools => { command => $samtools, version => $tool_versions{samtools}, identity => $tool_identities{samtools} },
				minimap2 => { command => $minimap2, version => $tool_versions{minimap2}, identity => $tool_identities{minimap2} },
				bam_filter => { command => $bam_filter, identity => $tool_identities{bam_filter} },
			},
			output => { alignment => $final_alignment, format => $row->{alignment_format}, index => $final_index },
			policies => {
				supplementary_alignments => $supplementary_alignments,
				allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
			},
			summary => $row, donors => \@unit_donors,
		});
		my @unit_outputs = ($unit->{json});
		push @unit_outputs, $final_alignment, $final_index if $final_alignment ne '-';
		push @unit_outputs, $published_origin_file if $published_origin_file ne '';
		write_checkpoint($unit->{stone},
			parameters => $unit->{parameters}, outputs => \@unit_outputs);
		push @checkpoint_outputs, $unit->{json};
	}
}

my $summary_partial = "$summary_file.part.$$";
push @publication_partials, $summary_partial;
open my $summary, '>', $summary_partial or die "Cannot create $summary_partial: $!\n";
print {$summary} join("\t", qw(mgs mode sample scopes status candidate_reads donor_reads aligned_records MM_records ML_records alignment_format alignment index)), "\n";
for my $row (sort {
	$a->{mgs} cmp $b->{mgs} || $a->{mode} cmp $b->{mode} || $a->{sample} cmp $b->{sample}
} @summary_rows) {
	print {$summary} join("\t", @{$row}{qw(mgs mode sample scopes status candidates donor aligned mm ml alignment_format alignment index)}), "\n";
}
close $summary or die "Cannot close $summary_partial: $!\n";
rename $summary_partial, $summary_file or die "Cannot publish $summary_file: $!\n";

my @manifest_units;
find({
	no_chdir => 1,
	wanted => sub {
		return unless -f $File::Find::name && $File::Find::name =~ /\.json\z/;
		my $json = $File::Find::name;
		(my $stone = $json) =~ s/\.json\z/.stone/;
		return unless -s $stone && checkpoint_valid($stone);
		my $unit = read_json_file($json);
		return unless ($unit->{format} // '') eq 'matafiler-meth2rep-unit-v2'
			&& ref($unit->{donors}) eq 'ARRAY' && ref($unit->{summary}) eq 'HASH';
		push @manifest_units, $unit;
	},
}, File::Spec->catdir($state_dir, 'units'));
my $manifest_partial = "$out_manifest.part.$$";
push @publication_partials, $manifest_partial;
open my $manifest_out, '>', $manifest_partial or die "Cannot create $manifest_partial: $!\n";
print {$manifest_out} join("\t", qw(mgs mode sample scope technology original_modbam donor_candidate_reads total_candidate_reads transferred_records MM_records ML_records status input_validation alignment_format alignment index input_fingerprint completed_epoch)), "\n";
for my $unit (sort {
	$a->{mgs} cmp $b->{mgs} || $a->{mode} cmp $b->{mode} || $a->{sample} cmp $b->{sample}
} @manifest_units) {
	for my $donor (sort {
		$a->{scope} cmp $b->{scope} || $a->{modbam} cmp $b->{modbam}
	} @{$unit->{donors}}) {
		my $row = $unit->{summary};
		my $unit_key = "$unit->{mode}\t$unit->{mgs}\t$unit->{sample}";
		my $input_validation = exists $unit_context{$unit_key}
			? 'current' : 'recorded_only';
		print {$manifest_out} join("\t",
			$unit->{mgs}, $unit->{mode}, $unit->{sample},
			$donor->{scope}, $donor->{technology}, $donor->{modbam},
			$donor->{candidate_reads}, $row->{candidates}, $row->{aligned},
			$row->{mm}, $row->{ml}, $row->{status}, $input_validation,
			$row->{alignment_format}, $row->{alignment}, $row->{index},
			$unit->{input_fingerprint}, $unit->{created_epoch},
		), "\n";
	}
}
close $manifest_out or die "Cannot close $manifest_partial: $!\n";
rename $manifest_partial, $out_manifest or die "Cannot publish $out_manifest: $!\n";

my %units_by_key = map { (join("\t", @{$_}{qw(mode mgs sample)}) => $_) } @manifest_units;
my %mgs_samples;
for my $unit_key (keys %unit_context) {
	my ($mode, $mgs, $sample) = split /\t/, $unit_key, 3;
	$mgs_samples{"$mgs\t$sample"} = 1;
}
for my $mgs_sample (sort keys %mgs_samples) {
	my ($mgs, $sample) = split /\t/, $mgs_sample, 2;
	my ($current_target) = grep { $_->{available} && $_->{mgs} eq $mgs }
		values %{$plan->{targets}};
	next unless $current_target;
	my ($recorded_source, $recorded_paths) = recorded_source_paths($sample, $map->{$sample});
	my @support_paths;
	if (($map->{$sample}{SupportReads} // '') ne '') {
		(undef, my $paths) = parseSupportReads($map->{$sample}{SupportReads});
		@support_paths = @{$paths};
	}
	my %recorded = map {
		my $canonical = -e $_ ? abs_path($_) : File::Spec->rel2abs($_);
		$canonical => 1
	} (@{$recorded_paths}, @support_paths);
	my (%modes_for_log, %declared_donors, @warnings);
	for my $mode (qw(mgs2rep rep2rep)) {
		my $unit_key = "$mode\t$mgs\t$sample";
		my $unit = $units_by_key{$unit_key} or next;
		if (($unit->{representative_mag} // '') ne $current_target->{representative_mag}) {
			push @warnings, "Historical $mode output targets a different representative and is not included in this log";
			next;
		}
		my $row = $unit->{summary};
		my $validation = exists($unit_context{$unit_key}) ? 'current' : 'recorded_only';
		my $alignment = $row->{alignment};
		$modes_for_log{$mode} = {
			status => $row->{status}, input_validation => $validation,
			checkpoint_reused_this_run => exists($cached_units{$unit_key}) ? JSON::PP::true : JSON::PP::false,
			candidate_reads => $row->{candidates}, donor_reads => $row->{donor},
			accepted_alignment_records => $row->{aligned},
			MM_records => $row->{mm}, ML_records => $row->{ml},
			alignment => $alignment, alignment_format => $row->{alignment_format},
			index => $row->{index},
			coverage => $unit->{coverage},
			processing => $unit->{processing},
			tools => $unit->{tools},
			scopes => $unit->{scope_reports},
			donors => $unit->{donors},
			completed_epoch => $unit->{created_epoch},
		};
		$declared_donors{$_->{modbam}} = 1 for @{$unit->{donors}};
		push @warnings, @{$unit->{warnings} || []};
	}
	for my $donor (sort keys %declared_donors) {
		push @warnings, "Declared original modBAM is absent from MATAFILER's recorded input paths (possibly relocated or generated separately): $donor"
			if (%recorded && !$recorded{$donor});
	}
	push @warnings, 'Neither per-sample input_raw.txt nor cohort Input_raw.txt was found; original donor paths could not be independently corroborated'
		if $recorded_source eq '';
	my %unique_warnings;
	@warnings = grep { !$unique_warnings{$_}++ } @warnings;
	my %sample_scopes;
	for my $target (values %{$plan->{targets}}) {
		next unless $target->{available} && $target->{mgs} eq $mgs;
		for my $scope_key (keys %{$target->{scope_keys}}) {
			$sample_scopes{$scope_key} = 1
				if $plan->{scopes}{$scope_key}{sample} eq $sample;
		}
	}
	my $log_file = File::Spec->catfile($out_dir, $mgs, file_component($sample) . '.meth2rep.json');
	make_path(dirname($log_file));
	write_json_file($log_file, {
		format => 'matafiler-meth2rep-sample-v2', component_version => $VERSION,
		mgs => $mgs, sample => $sample,
		representative_mag => $current_target->{representative_mag},
		representative_fasta => $current_target->{reference},
		generated_epoch => 0 + time,
		pipeline_input_record => $recorded_source || undef,
		pipeline_recorded_inputs => $recorded_paths,
		mapping_support_inputs => \@support_paths,
		provenance_rule => 'Manifest donor paths identify source files; each candidate QNAME and native-orientation SEQ digest must also match the assembly CRAM before methylation tags can transfer',
		shared_input_processing_this_run => {
			source_cram_scan_seconds => {
				map { ($_ => $source_scan_seconds{$_}) }
					grep { exists $source_scan_seconds{$_} } sort keys %sample_scopes
			},
			original_modbam_stream_seconds => {
				map { ($_ => $donor_stream_seconds{$_}) }
					grep { exists $donor_stream_seconds{$_} } sort keys %declared_donors
			},
			note => 'Shared scans may serve several MGS/mode units; do not sum these as per-unit exclusive CPU time',
		},
		filters => {
			historical_source_alignment => $HISTORICAL_SOURCE_FILTER,
			source_min_mapq => $source_min_mapq,
			source_min_coverage => $source_min_coverage,
			ont_target => \@filter_ont, pb_target => \@filter_pb,
		},
		mapper_presets => { ONT => $preset_ont, PB => $preset_pb },
		output_format => $output_format,
		policies => {
			supplementary_alignments => $supplementary_alignments,
			allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
		},
		modes => \%modes_for_log,
		warnings => \@warnings,
	});
	push @checkpoint_outputs, $log_file;
}

my @input_records = map {
	my @stat = stat($_);
	+{ path => $_, size => 0 + $stat[7], mtime => 0 + $stat[9] };
} @{$plan->{inputs}};
my $provenance = {
	format => 'matafiler-meth2rep-v2', component_version => $VERSION,
	created_epoch => 0 + time, modes => \@modes,
	selected_mgs => [sort keys %{$selected}],
	semantics => {
		candidate_source => 'primary nonduplicate QC-passing records retained in the already-filtered assembly CRAM and aligned to source-MAG contigs',
		competition_universe => 'the retained MATAFILER co-assembly represented by the source CRAM; no uniqueness claim is made against sequences absent from that assembly',
		donor_source => 'original manifest-declared modBAM primary records',
		identity => 'candidate QNAME plus exact native-orientation sequence must agree between the assembly CRAM and declared original modBAM',
		alignment => "fresh original-sequence minimap2 alignment to the representative MAG with -Y and --secondary=no; supplementary output policy is $supplementary_alignments",
		modification_projection => 'first-party exact full-native-sequence MM/ML transfer; MM/ML remain in original read coordinates, MN is regenerated, hard clipping and sequence mismatch fail',
	},
	filters => {
		historical_source_alignment => $HISTORICAL_SOURCE_FILTER,
		source_min_mapq => $source_min_mapq,
		source_min_coverage => $source_min_coverage,
		ont => join(' ', @filter_ont), pb => join(' ', @filter_pb),
	},
	mapper_presets => { ONT => $preset_ont, PB => $preset_pb },
	output => { format => $output_format, cram_reference => $output_format eq 'cram' ? 'embedded' : 'not_applicable' },
	policies => {
		supplementary_alignments => $supplementary_alignments,
		allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
	},
	resources => { threads => $threads, sort_memory_budget_gb => $memory_gb, samtools_sort_memory_per_thread => $sort_memory },
	tools => {
		samtools => { command => $samtools, version => $tool_versions{samtools}, identity => $tool_identities{samtools} },
		minimap2 => { command => $minimap2, version => $tool_versions{minimap2}, identity => $tool_identities{minimap2} },
		bam_filter => { command => $bam_filter, identity => $tool_identities{bam_filter} },
	},
	inputs => \@input_records,
	input_fingerprint => $checkpoint_parameters{input_fingerprint},
};
my $provenance_partial = "$provenance_file.part.$$";
push @publication_partials, $provenance_partial;
open my $provenance_fh, '>', $provenance_partial or die "Cannot create $provenance_partial: $!\n";
print {$provenance_fh} JSON::PP->new->ascii->canonical->pretty->encode($provenance)
	or die "Cannot write $provenance_partial: $!\n";
close $provenance_fh or die "Cannot close $provenance_partial: $!\n";
rename $provenance_partial, $provenance_file or die "Cannot publish $provenance_file: $!\n";
push @checkpoint_outputs, $summary_file, $provenance_file, $out_manifest;
write_checkpoint(
	$checkpoint_file,
	parameters => \%checkpoint_parameters,
	outputs => \@checkpoint_outputs,
);
print "meth2rep v$VERSION: completed; summary $summary_file\n";
#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../..";
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use File::Spec;
use File::Basename qw(dirname);
use File::Find qw(find);
use Cwd qw(abs_path);
use File::Temp qw(tempdir);
use Digest::SHA qw(sha1_hex sha256_hex);
use DB_File;
use Fcntl qw(O_CREAT O_RDWR LOCK_EX LOCK_NB);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP;
use Text::ParseWords qw(shellwords);
use Time::HiRes ();

use Mods::Checkpoint qw(write_checkpoint checkpoint_valid read_checkpoint);
use Mods::GenoMetaAss qw(getDirsPerAssmblGrp parseSupportReads);
use Mods::Meth2Rep qw(
	read_target_mgs read_modbam_manifest build_meth2rep_plan
);

my $VERSION = '0.6';
my @publication_partials;
END {
	for my $file (@publication_partials) {
		unlink $file if defined($file) && -f $file;
	}
}

sub usage {
	return <<'USAGE';
Usage:
  meth2rep.pl --mgs-dir DIR --map FILE[,FILE...]
    --modbam-manifest FILE (--mgs MGS.1,MGS.2 | --mgs-file FILE)
    (--mgs2rep | --rep2rep | both) [-o DIR]

Modes (either or both):
  --mgs2rep  Candidate reads supporting any MAG in an MGS are freshly aligned
             to that MGS representative MAG.
  --rep2rep  Only candidate reads supporting the representative MAG are freshly
             aligned to that same representative MAG.

Required manifest header (tab separated):
  sample  scope  technology  modbam

scope must be primary or support; technology must be ONT or PB. Repeated
sample/scope rows declare multiple donor BAMs. A QNAME must be unique across
those BAMs for its sample/scope; ambiguous identities stop the run.

Controls:
  --source-min-mapq INT       Minimum MAPQ in assembly CRAMs (default 10)
  --source-min-coverage FLOAT Minimum source aligned-query fraction (default 0.5)
  --target-min-mapq INT       Minimum representative MAPQ for both technologies
  --target-min-coverage FLOAT Minimum aligned query fraction for both technologies
  --target-max-edit-rate FLOAT Maximum representative edit rate for both technologies
  --target-min-end-clip INT   Minimum two-ended clipping rejected by bamFilter
  --mapper-filter-ont STRING  bamFilter arguments (default "0.15 0.5 10 0")
  --mapper-filter-pb STRING   bamFilter arguments (default "0.05 0.5 30 0")
  --minimap2-preset-ont NAME  ONT preset (default map-ont; modern option: lr:hq)
  --minimap2-preset-pb NAME   PacBio preset (default map-pb; HiFi option: map-hifi)
  --supplementary-alignments POLICY
                              drop (default) or keep non-overlapping records
  --allow-missing-mn          Permit legacy donors without MN (warning is recorded)
  --output-format FORMAT      bam (default) or self-contained cram
  --threads INT               Mapping/sorting threads (default 4)
  --memory-gb INT             Sort-memory sizing budget, not a process cap (default 32)
  --tmp DIR                   Parent for disposable intermediates
  --keep-read-ids             Retain gzip-compressed candidate-name audit files
  --out-manifest FILE         Completed-unit/donor TSV (default OUT/manifest.tsv)
  --override                  Rebuild selected units even when their checkpoints validate
  --redo                      Alias of --override
  --plan-only                 Run input/tool/output preflight and write a preview plan
  --samtools COMMAND          Override samtools (otherwise resolved from PATH)
  --minimap2 COMMAND          Override minimap2 (otherwise resolved from PATH)
  --bam-filter COMMAND        Override bundled first-party bamFilter.pl
USAGE
}

sub shell_quote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub tool_command {
	my ($tool, @arguments) = @_;
	my @prefix = shellwords($tool // '');
	die "Empty external-tool command\n" unless @prefix;
	return join(' ', map { shell_quote($_) } (@prefix, @arguments));
}

sub find_on_path {
	my ($program) = @_;
	return unless defined($program) && $program ne '';
	for my $directory (split /:/, ($ENV{PATH} // '')) {
		$directory = '.' if $directory eq '';
		my $candidate = File::Spec->catfile($directory, $program);
		return abs_path($candidate) || File::Spec->rel2abs($candidate)
			if -f $candidate && -x $candidate;
	}
	return;
}

sub file_sha256 {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot hash $file: $!\n";
	binmode $fh;
	my $digest = Digest::SHA->new(256);
	$digest->addfile($fh);
	close $fh or die "Cannot close $file after hashing: $!\n";
	return $digest->hexdigest;
}

sub canonical_command {
	my ($command, $label) = @_;
	my @tokens = shellwords($command // '');
	die "$label command is empty\n" unless @tokens;
	my $executable;
	if (File::Spec->file_name_is_absolute($tokens[0]) || $tokens[0] =~ m{/}) {
		my $candidate = File::Spec->rel2abs($tokens[0]);
		$executable = abs_path($candidate) || $candidate;
	} else {
		$executable = find_on_path($tokens[0]);
	}
	die "$label executable '$tokens[0]' is missing or not executable\n"
		unless defined($executable) && -f $executable && -x $executable;
	$tokens[0] = $executable;
	for my $index (1 .. $#tokens) {
		next unless -f $tokens[$index];
		$tokens[$index] = abs_path($tokens[$index]) || File::Spec->rel2abs($tokens[$index]);
	}
	return join(' ', map { shell_quote($_) } @tokens);
}

sub command_identity {
	my ($command) = @_;
	my @tokens = shellwords($command);
	my @parts = ('argv=' . join("\x1f", @tokens));
	for my $token (@tokens) {
		next unless -f $token;
		my @stat = stat($token);
		push @parts, join(':', $token, $stat[7], $stat[9], file_sha256($token));
	}
	return sha256_hex(join("\n", @parts));
}

sub run_shell {
	my ($description, $command) = @_;
	my $status = system('bash', '-o', 'pipefail', '-c', $command);
	return if $status == 0;
	my $detail = $status == -1 ? "could not execute: $!"
		: ($status & 127) ? 'terminated by signal ' . ($status & 127)
		: 'exit code ' . ($status >> 8);
	die "$description failed ($detail)\n";
}

sub open_command {
	my ($description, $command) = @_;
	open my $fh, '-|', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub open_sink_command {
	my ($description, $command) = @_;
	open my $fh, '|-', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub close_command {
	my ($fh, $description) = @_;
	return if close $fh;
	my $status = $?;
	my $detail = $status == -1 ? "could not execute: $!"
		: ($status & 127) ? 'terminated by signal ' . ($status & 127)
		: 'exit code ' . ($status >> 8);
	die "$description failed ($detail)\n";
}

sub capture_command {
	my ($description, $command) = @_;
	my $fh = open_command($description, $command);
	local $/;
	my $output = <$fh> // '';
	close_command($fh, $description);
	$output =~ s/[\r\n]+$//;
	return $output;
}

sub sort_unique {
	my ($input, $output, $tmp_parent) = @_;
	local $ENV{LC_ALL} = 'C';
	my @command = ('sort', '-u');
	push @command, ('-T', $tmp_parent) if defined($tmp_parent) && -d $tmp_parent;
	push @command, ('-o', $output, $input);
	my $status = system @command;
	die "Sorting read names failed (exit " . ($status == -1 ? -1 : $status >> 8) . ")\n"
		if $status != 0;
}

sub count_lines {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot count $file: $!\n";
	my $count = 0;
	$count++ while <$fh>;
	close $fh or die "Cannot close $file: $!\n";
	return $count;
}

sub concatenate_files {
	my ($output, @inputs) = @_;
	open my $out, '>', $output or die "Cannot create $output: $!\n";
	for my $input (@inputs) {
		open my $in, '<', $input or die "Cannot read $input: $!\n";
		while (read($in, my $buffer, 1024 * 1024)) {
			print {$out} $buffer or die "Cannot write $output: $!\n";
		}
		close $in or die "Cannot close $input: $!\n";
	}
	close $out or die "Cannot close $output: $!\n";
}

sub remove_temp_files {
	for my $file (@_) {
		next unless defined($file) && $file ne '' && -e $file;
		unlink $file or die "Cannot remove temporary file $file: $!\n";
	}
}

sub filter_bam_by_names {
	my (%options) = @_;
	my $database = "$options{workdir}/names." . sha1_hex(
		join("\0", $options{names}, $options{input}, $options{output}),
	) . '.db';
	my %wanted;
	tie %wanted, 'DB_File', $database, O_RDWR | O_CREAT, 0600, $DB_HASH
		or die "Cannot create disk-backed read-name index $database: $!\n";
	open my $names, '<', $options{names} or die "Cannot read $options{names}: $!\n";
	while (my $name = <$names>) {
		$name =~ s/[\r\n]+$//;
		$wanted{$name} = 1 if $name ne '';
	}
	close $names or die "Cannot close $options{names}: $!\n";

	my $source_description = "reading $options{input} for name filtering";
	my $source = open_command(
		$source_description,
		tool_command($options{samtools}, 'view', '-h', $options{input}),
	);
	my $sink_description = "writing filtered BAM $options{output}";
	my $sink = open_sink_command(
		$sink_description,
		tool_command($options{samtools}, 'view', '-b', '-o', $options{output}, '-'),
	);
	while (my $line = <$source>) {
		if ($line =~ /^\@/) {
			print {$sink} $line or die "Cannot stream header to $options{output}: $!\n";
			next;
		}
		my ($name, $flag) = split /\t/, $line, 3;
		next unless defined($name) && exists $wanted{$name};
		next if ($options{exclude_flags} || 0) && ($flag & $options{exclude_flags});
		print {$sink} $line or die "Cannot stream record to $options{output}: $!\n";
	}
	close_command($source, $source_description);
	close_command($sink, $sink_description);
	untie %wanted or die "Cannot close disk-backed read-name index $database: $!\n";
	unlink $database or die "Cannot remove temporary read-name index $database: $!\n" if -e $database;
}

sub assert_name_subset {
	my ($wanted_file, $found_file, $context) = @_;
	open my $wanted, '<', $wanted_file or die "Cannot read $wanted_file: $!\n";
	open my $found, '<', $found_file or die "Cannot read $found_file: $!\n";
	my $wanted_name = <$wanted>;
	my $found_name = <$found>;
	chomp $wanted_name if defined $wanted_name;
	chomp $found_name if defined $found_name;
	my (@missing, $missing_count);
	while (defined $wanted_name) {
		while (defined($found_name) && $found_name lt $wanted_name) {
			$found_name = <$found>;
			chomp $found_name if defined $found_name;
		}
		if (!defined($found_name) || $found_name ne $wanted_name) {
			$missing_count++;
			push @missing, $wanted_name if @missing < 5;
		}
		$wanted_name = <$wanted>;
		chomp $wanted_name if defined $wanted_name;
	}
	close $wanted or die "Cannot close $wanted_file: $!\n";
	close $found or die "Cannot close $found_file: $!\n";
	if ($missing_count) {
		die "$context: $missing_count candidate read name(s) are absent from the declared original modBAM set"
			. (@missing ? ' (examples: ' . join(', ', @missing) . ')' : '')
			. ". Refusing a silently incomplete methylation result.\n";
	}
}

sub file_component {
	my ($value) = @_;
	$value =~ s/([^A-Za-z0-9_.-])/sprintf('_%02X', ord($1))/ge;
	return $value;
}

sub native_sequence {
	my ($fields, $origin) = @_;
	my ($name, $flag, $cigar, $sequence) = @{$fields}[0, 1, 5, 9];
	die "$origin record '$name' has no complete SEQ\n"
		if !defined($sequence) || $sequence eq '' || $sequence eq '*';
	die "$origin record '$name' is hard-clipped ($cigar); its full native sequence is unavailable\n"
		if defined($cigar) && $cigar =~ /H/;
	if ($flag & 0x10) {
		$sequence = reverse $sequence;
		$sequence =~ tr/ACGTRYKMSWBDHVNacgtrykmswbdhvn/TGCAYRMKSWVHDBNtgcayrmkswvhdbn/;
	}
	return uc $sequence;
}

sub sort_identity_records {
	my ($input, $output, $tmp_parent) = @_;
	local $ENV{LC_ALL} = 'C';
	my @command = ('sort', '-k1,1');
	push @command, ('-T', $tmp_parent) if defined($tmp_parent) && -d $tmp_parent;
	push @command, ('-o', $output, $input);
	my $status = system @command;
	die "Sorting read-identity records failed (exit " . ($status == -1 ? -1 : $status >> 8) . ")\n"
		if $status != 0;
}

sub identity_record {
	my ($fh, $file) = @_;
	my $line = <$fh>;
	return unless defined $line;
	$line =~ s/[\r\n]+\z//;
	my ($name, $digest, @extra) = split /\t/, $line, -1;
	die "Malformed read-identity record in $file\n"
		unless defined($name) && $name ne '' && defined($digest)
		&& $digest =~ /\A[0-9a-f]{64}\z/ && !@extra;
	return [$name, $digest];
}

sub assert_source_donor_identity {
	my ($source_file, $donor_file, $context) = @_;
	open my $source_fh, '<', $source_file or die "Cannot read $source_file: $!\n";
	open my $donor_fh, '<', $donor_file or die "Cannot read $donor_file: $!\n";
	my $source = identity_record($source_fh, $source_file);
	my $donor = identity_record($donor_fh, $donor_file);
	my ($last_source, $last_donor);
	while (defined $source) {
		die "$context: source CRAM contains more than one primary record named '$source->[0]'\n"
			if defined($last_source) && $last_source eq $source->[0];
		while (defined($donor) && $donor->[0] lt $source->[0]) {
			die "$context: declared donor set contains more than one record named '$donor->[0]'\n"
				if defined($last_donor) && $last_donor eq $donor->[0];
			$last_donor = $donor->[0];
			$donor = identity_record($donor_fh, $donor_file);
		}
		die "$context: candidate read '$source->[0]' is absent from the declared original modBAM set\n"
			unless defined($donor) && $donor->[0] eq $source->[0];
		die "$context: native SEQ for '$source->[0]' differs between the assembly CRAM and declared original modBAM; QNAME alone is not sufficient provenance\n"
			unless $donor->[1] eq $source->[1];
		$last_source = $source->[0];
		$last_donor = $donor->[0];
		$source = identity_record($source_fh, $source_file);
		$donor = identity_record($donor_fh, $donor_file);
	}
	close $source_fh or die "Cannot close $source_file: $!\n";
	close $donor_fh or die "Cannot close $donor_file: $!\n";
}

sub bam_count {
	my ($samtools, $bam) = @_;
	my $count = capture_command(
		"counting records in $bam",
		tool_command($samtools, 'view', '-c', $bam),
	);
	die "samtools returned a non-integer record count for $bam: '$count'\n"
		unless $count =~ /^\d+$/;
	return 0 + $count;
}

sub bam_tag_counts {
	my ($samtools, $bam) = @_;
	my $fh = open_command("reading tags from $bam", tool_command($samtools, 'view', $bam));
	my ($records, $mm, $ml) = (0, 0, 0);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		my @fields = split /\t/, $line, -1;
		die "Malformed transferred SAM record from $bam\n" unless @fields >= 11;
		my (@mm_tags, @ml_tags, @mn_tags);
		for my $field (@fields[11 .. $#fields]) {
			push @mm_tags, $field if $field =~ /\AM[Mm]:/;
			push @ml_tags, $field if $field =~ /\AM[Ll]:/;
			push @mn_tags, $field if $field =~ /\AMN:/;
		}
		die "Transferred record '$fields[0]' in $bam does not contain exactly one MM:Z and one ML:B:C tag\n"
			unless @mm_tags == 1 && @ml_tags == 1
			&& $mm_tags[0] =~ /\AM[Mm]:Z:/ && $ml_tags[0] =~ /\AM[Ll]:B:C(?:,\d+)*\z/;
		die "Transferred record '$fields[0]' in $bam has missing or stale MN for SEQ length "
			. length($fields[9]) . "\n"
			unless @mn_tags == 1 && $mn_tags[0] =~ /\AMN:i:(\d+)\z/
			&& $1 == length($fields[9]);
		$records++;
		$mm++;
		$ml++;
	}
	close_command($fh, "reading tags from $bam");
	return ($records, $mm, $ml);
}

sub validate_donor_subset {
	my ($samtools, $bam, $found_raw, $identity_raw, $allow_missing_mn) = @_;
	my $fh = open_command("validating donor records in $bam", tool_command($samtools, 'view', $bam));
	open my $names, '>', $found_raw or die "Cannot create $found_raw: $!\n";
	my $identities;
	if (defined($identity_raw) && $identity_raw ne '') {
		open $identities, '>', $identity_raw or die "Cannot create $identity_raw: $!\n";
	}
	my ($previous, $records, $mm, $ml, $missing_mn);
	$records = $mm = $ml = $missing_mn = 0;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		my @fields = split /\t/, $line, -1;
		die "Malformed donor SAM record from $bam\n" unless @fields >= 11;
		my ($name, $flag) = @fields[0, 1];
		die "Donor modBAM $bam contains a paired record '$name'; meth2rep only supports singleton ONT/PB records\n"
			if $flag & 0x1;
		die "Donor modBAM $bam contains more than one primary record named '$name'; CRAM-to-donor identity is ambiguous\n"
			if defined($previous) && $previous eq $name;
		my $sequence = native_sequence(\@fields, 'Donor');
		my (@mm_fields, @ml_fields, @mn_fields);
		for my $field (@fields[11 .. $#fields]) {
			push @mm_fields, $field if $field =~ /\AM[Mm]:/;
			push @ml_fields, $field if $field =~ /\AM[Ll]:/;
			push @mn_fields, $field if $field =~ /\AMN:/;
		}
		die "Donor record '$name' must contain exactly one paired MM:Z/ML:B:C tag set\n"
			unless @mm_fields == 1 && @ml_fields == 1
			&& $mm_fields[0] =~ /\AM[Mm]:Z:(?:[ACGTUN][+-](?:[A-Za-z]+|\d+)[.?]?(?:,\d+)*;)*\z/
			&& $ml_fields[0] =~ /\AM[Ll]:B:C(?:,\d+)*\z/;
		die "Donor record '$name' contains duplicate MN tags\n" if @mn_fields > 1;
		if (!@mn_fields) {
			die "Donor record '$name' lacks MN; use --allow-missing-mn only for explicitly reviewed legacy input\n"
				unless $allow_missing_mn;
			$missing_mn++;
		} else {
			die "Donor record '$name' has malformed or stale MN\n"
				unless $mn_fields[0] =~ /\AMN:i:(\d+)\z/ && $1 == length($fields[9]);
		}
		$previous = $name;
		print {$names} "$name\n" or die "Cannot write $found_raw: $!\n";
		if ($identities) {
			print {$identities} "$name\t", sha256_hex($sequence), "\n"
				or die "Cannot write $identity_raw: $!\n";
		}
		$records++;
		$mm++;
		$ml++;
	}
	close $names or die "Cannot close $found_raw: $!\n";
	close $identities or die "Cannot close $identity_raw: $!\n" if $identities;
	close_command($fh, "validating donor records in $bam");
	return ($records, $mm, $ml, $missing_mn);
}

sub validate_filter {
	my ($value, $name) = @_;
	my @parts = split /\s+/, $value;
	die "$name must contain four bamFilter values: max_edit_rate min_query_coverage min_mapq min_end_clip\n"
		unless @parts == 4
		&& $parts[0] =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $parts[0] >= 0 && $parts[0] <= 1
		&& $parts[1] =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $parts[1] >= 0 && $parts[1] <= 1
		&& $parts[2] =~ /^\d+$/ && $parts[2] <= 255
		&& $parts[3] =~ /^\d+$/;
	return @parts;
}

sub cigar_query_coverage {
	my ($cigar, $context) = @_;
	my ($aligned, $length) = (0, 0);
	my $reconstructed = '';
	while ($cigar =~ /([1-9]\d*)([MIDNSHP=X])/g) {
		my ($span, $op) = ($1, $2);
		$reconstructed .= "$span$op";
		if ($op =~ /[MIS=XH]/) { $length += $span; }
		if ($op =~ /[MI=X]/) { $aligned += $span; }
	}
	die "Malformed source CIGAR '$cigar' for $context\n"
		unless $cigar ne '*' && $reconstructed eq $cigar && $length > 0;
	return $aligned / $length;
}

sub fingerprint_inputs {
	my ($inputs, $parameters) = @_;
	my @records;
	for my $file (@{$inputs}) {
		my @stat = stat($file);
		die "Cannot fingerprint meth2rep input $file\n" unless @stat;
		push @records, join("\t", $file, map { 0 + $stat[$_] } (0, 1, 7, 9, 10));
	}
	push @records, map { "option\t$_\t$parameters->{$_}" } sort keys %{$parameters};
	return sha256_hex(join("\n", @records));
}

sub require_mgs_readiness {
	my ($directory) = @_;
	die "MGS directory is missing: $directory\n" unless -d $directory;
	my %stages = (
		Stage1 => 'stage-1',
		BinExtr => 'extract-bin-contigs',
	);
	my %manifests;
	for my $name (sort keys %stages) {
		my $stone = File::Spec->catfile($directory, 'LOGandSUB', 'checkpoints', "$name.stone");
		die "Required MGS progress checkpoint is missing or legacy-empty: $stone\n" unless -s $stone;
		my $manifest = read_checkpoint($stone)
			or die "MGS progress checkpoint is malformed: $stone\n";
		die "MGS progress checkpoint has wrong stage: $stone\n"
			unless ($manifest->{parameters}{stage} // '') eq $stages{$name};
		die "MGS progress checkpoint has stale/missing outputs: $stone\n"
			unless checkpoint_valid($stone, parameters => { stage => $stages{$name} });
		$manifests{$name} = $manifest;
	}
	return \%manifests;
}

sub read_json_file {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read $file: $!\n";
	local $/;
	my $data = eval { JSON::PP->new->decode(<$fh> // '') };
	close $fh or die "Cannot close $file: $!\n";
	return $data if !$@ && ref($data) eq 'HASH';
	die "Invalid meth2rep unit JSON: $file\n";
}

sub write_json_file {
	my ($file, $data) = @_;
	my $partial = "$file.part.$$";
	push @publication_partials, $partial;
	open my $fh, '>', $partial or die "Cannot create $partial: $!\n";
	print {$fh} JSON::PP->new->ascii->canonical->pretty->encode($data)
		or die "Cannot write $partial: $!\n";
	close $fh or die "Cannot close $partial: $!\n";
	rename $partial, $file or die "Cannot publish $file: $!\n";
}

sub intersection_count {
	my ($left_file, $right_file) = @_;
	open my $left, '<', $left_file or die "Cannot read $left_file: $!\n";
	open my $right, '<', $right_file or die "Cannot read $right_file: $!\n";
	my ($a, $b, $count) = (scalar(<$left>), scalar(<$right>), 0);
	while (defined($a) && defined($b)) {
		chomp($a, $b);
		if ($a eq $b) { $count++; $a = <$left>; $b = <$right>; }
		elsif ($a lt $b) { $a = <$left>; }
		else { $b = <$right>; }
	}
	close $left or die "Cannot close $left_file: $!\n";
	close $right or die "Cannot close $right_file: $!\n";
	return $count;
}

sub write_intersection_origin_rows {
	my ($candidate_file, $donor_file, $sink, $sample, $scope, $modbam) = @_;
	open my $candidate, '<', $candidate_file or die "Cannot read $candidate_file: $!\n";
	open my $donor, '<', $donor_file or die "Cannot read $donor_file: $!\n";
	my ($a, $b) = (scalar(<$candidate>), scalar(<$donor>));
	while (defined($a) && defined($b)) {
		chomp($a, $b);
		if ($a eq $b) {
			print {$sink} join("\t", $sample, $scope, $a, $modbam), "\n"
				or die "Cannot write compressed read-origin rows: $GzipError\n";
			$a = <$candidate>; $b = <$donor>;
		} elsif ($a lt $b) { $a = <$candidate>; }
		else { $b = <$donor>; }
	}
	close $candidate or die "Cannot close $candidate_file: $!\n";
	close $donor or die "Cannot close $donor_file: $!\n";
}

sub unit_paths {
	my ($state_dir, $target, $sample) = @_;
	my $base = File::Spec->catfile($state_dir, 'units', $target->{mode}, $target->{mgs}, file_component($sample));
	return ("$base.json", "$base.stone");
}

sub output_alignment_path {
	my ($out_dir, $target, $sample, $format) = @_;
	my $name = join('__', file_component($sample), $target->{mode},
		file_component($target->{representative_mag})) . ".mod.$format";
	return File::Spec->catfile($out_dir, $target->{mgs}, $name);
}

sub alignment_index_path {
	my ($alignment, $format) = @_;
	return $alignment . ($format eq 'cram' ? '.crai' : '.bai');
}

sub reference_lengths {
	my ($fasta) = @_;
	my $fh = IO::Uncompress::Gunzip->new($fasta, Transparent => 1, MultiStream => 1)
		or die "Cannot read representative FASTA $fasta: $GunzipError\n";
	my (%lengths, $name);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		if ($line =~ /^>(\S+)/) {
			$name = $1;
			die "Representative FASTA $fasta repeats contig '$name'\n" if exists $lengths{$name};
			$lengths{$name} = 0;
		} elsif ($line ne '') {
			die "Representative FASTA $fasta has sequence before a header\n" unless defined $name;
			$line =~ s/\s+//g;
			$lengths{$name} += length($line);
		}
	}
	close $fh or die "Cannot close representative FASTA $fasta: $!\n";
	die "Representative FASTA $fasta contains no nonempty contigs\n"
		unless keys(%lengths) && !grep { $_ <= 0 } values %lengths;
	return \%lengths;
}

sub bam_coverage {
	my ($samtools, $bam, $lengths) = @_;
	my $total_bases = 0;
	$total_bases += $_ for values %{$lengths};
	my ($covered_bases, $depth_sum) = (0, 0);
	my $command = tool_command($samtools, 'depth', '-d', 0, '-q', 0, '-Q', 0, $bam);
	my $fh = open_command("computing reference coverage for $bam", $command);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		my ($contig, $position, $depth) = split /\t/, $line, -1;
		die "Malformed samtools depth result for $bam: $line\n"
			unless defined($depth) && exists($lengths->{$contig})
			&& $position =~ /^\d+$/ && $position >= 1 && $position <= $lengths->{$contig}
			&& $depth =~ /^\d+$/;
		$covered_bases++ if $depth > 0;
		$depth_sum += $depth;
	}
	close_command($fh, "computing reference coverage for $bam");
	die "Coverage calculation for $bam reported $covered_bases positions across only $total_bases reference bases\n"
		if $covered_bases > $total_bases;
	return {
		reference_bases => $total_bases, covered_bases => $covered_bases,
		breadth_fraction => $total_bases ? $covered_bases / $total_bases : 0,
		mean_depth => $total_bases ? $depth_sum / $total_bases : 0,
		depth_sum => $depth_sum,
		method => 'samtools depth -d 0 -q 0 -Q 0; reference denominator includes every representative contig',
	};
}

sub read_filter_report {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read alignment-filter report $file: $!\n";
	my %patterns = (
		input_records => qr/^Input records:\s*(\d+)/,
		retained_mapped => qr/^Retained mapped records:\s*(\d+)/,
		newly_filtered => qr/^Newly filtered records:\s*(\d+)/,
		malformed => qr/^Malformed SAM records skipped:\s*(\d+)/,
		mapq_rejected => qr/^\s*Mapping quality \([^)]*\):\s*(\d+)/,
		coverage_rejected => qr/^\s*Query coverage \([^)]*\):\s*(\d+)/,
		edit_rate_rejected => qr/^\s*Edit rate \([^)]*\):\s*(\d+)/,
		end_clip_rejected => qr/^\s*Both ends clipped \([^)]*\):\s*(\d+)/,
	);
	my (%values, @warnings);
	while (my $line = <$fh>) {
		for my $key (keys %patterns) {
			$values{$key} = 0 + $1 if $line =~ $patterns{$key};
		}
		push @warnings, $line if $line =~ /\b(?:warning|error|malformed)\b/i
			&& $line !~ /^Malformed SAM records skipped:/;
	}
	close $fh or die "Cannot close $file: $!\n";
	for my $required (qw(input_records retained_mapped newly_filtered malformed)) {
		die "Alignment-filter report $file lacks '$required' statistics\n"
			unless exists $values{$required};
	}
	die "Alignment filter skipped $values{malformed} malformed SAM record(s); refusing an incomplete output. "
		. join(' | ', @warnings[0 .. ($#warnings < 4 ? $#warnings : 4)]) . "\n"
		if $values{malformed};
	return (\%values, \@warnings);
}

sub read_transfer_stats {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read tag-transfer statistics $file: $!\n";
	my %values;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+\z//;
		my ($key, $value, @extra) = split /\t/, $line, -1;
		die "Malformed tag-transfer statistics in $file\n"
			unless defined($key) && $key =~ /\A[a-z_]+\z/
			&& defined($value) && $value =~ /\A\d+\z/ && !@extra
			&& !exists($values{$key});
		$values{$key} = 0 + $value;
	}
	close $fh or die "Cannot close $file: $!\n";
	for my $required (qw(input_alignments output_alignments output_read_names
		dropped_no_primary_groups dropped_no_primary_alignments
		dropped_supplementary_alignments stripped_sa_tags legacy_missing_mn_reads)) {
		die "Tag-transfer statistics $file lacks '$required'\n"
			unless exists $values{$required};
	}
	return \%values;
}

sub materialize_cram_reference {
	my ($source, $workdir, $samtools, $cache) = @_;
	return $cache->{$source} if exists $cache->{$source};
	my $destination = File::Spec->catfile($workdir,
		'reference.' . sha1_hex($source) . '.fa');
	my $input = IO::Uncompress::Gunzip->new($source, Transparent => 1, MultiStream => 1)
		or die "Cannot read representative FASTA $source for CRAM encoding: $GunzipError\n";
	open my $output, '>', $destination or die "Cannot create $destination: $!\n";
	while (my $line = <$input>) {
		print {$output} $line or die "Cannot write $destination: $!\n";
	}
	close $input or die "Cannot finish reading representative FASTA $source: $GunzipError\n";
	close $output or die "Cannot close $destination: $!\n";
	run_shell("indexing temporary CRAM reference", tool_command($samtools, 'faidx', $destination));
	$cache->{$source} = $destination;
	return $destination;
}

sub recorded_source_paths {
	my ($sample, $record) = @_;
	my $per_sample = File::Spec->catfile($record->{wrdir}, 'input_raw.txt');
	my @paths;
	my $source = '';
	if (-s $per_sample) {
		open my $fh, '<', $per_sample or die "Cannot read $per_sample: $!\n";
		local $/;
		@paths = split /;/, (<$fh> // '');
		close $fh or die "Cannot close $per_sample: $!\n";
		$source = $per_sample;
	} else {
		my $cohort = File::Spec->catfile(dirname($record->{wrdir}), 'Input_raw.txt');
		if (-s $cohort) {
			open my $fh, '<', $cohort or die "Cannot read $cohort: $!\n";
			while (my $line = <$fh>) {
				my ($row_sample, $value) = split /\t/, $line, 2;
				next unless defined($value) && $row_sample eq $sample;
				@paths = split /;/, $value;
				$source = $cohort;
				last;
			}
			close $fh or die "Cannot close $cohort: $!\n";
		}
	}
	for (@paths) { s/^\s+|\s+$//g; }
	@paths = grep { $_ ne '' } @paths;
	return ($source, \@paths);
}

sub write_plan {
	my ($file, $plan) = @_;
	my $partial = "$file.part.$$";
	push @publication_partials, $partial;
	open my $fh, '>', $partial or die "Cannot create $partial: $!\n";
	print {$fh} join("\t", qw(mgs mode representative_mag representative_fasta source_mag assembly_group assembly bin_assignment source_contigs eligible_scopes eligible_sample_scopes status)), "\n";
	for my $key (sort keys %{$plan->{targets}}) {
		my $target = $plan->{targets}{$key};
		next unless $target->{available};
		for my $source (@{$target->{sources}}) {
			my @source_scope_keys = sort grep {
				$plan->{scopes}{$_}{group} eq $source->{group}
			} keys %{$target->{scope_keys}};
			my @source_scopes = map {
				"$plan->{scopes}{$_}{sample}:$plan->{scopes}{$_}{scope}"
			} @source_scope_keys;
			print {$fh} join("\t",
				$target->{mgs}, $target->{mode}, $target->{representative_mag}, $target->{reference},
				$source->{mag}, $source->{group}, $source->{assembly}, $source->{assignment},
				scalar(@{$source->{contigs}}), scalar(@source_scopes), join(',', @source_scopes), 'ready',
			), "\n";
		}
	}
	for my $row (@{$plan->{unavailable}}) {
		print {$fh} join("\t", $row->{mgs}, $row->{mode}, '-', '-', '-', '-', '-', '-', 0, 0, '-', $row->{reason}), "\n";
	}
	close $fh or die "Cannot close $partial: $!\n";
	rename $partial, $file or die "Cannot publish $file: $!\n";
}

sub validate_cram_reference {
	my ($scope) = @_;
	(my $stat_file = $scope->{cram}) =~ s/\.cram\z/.reference.stat/;
	return unless -s $stat_file;
	open my $fh, '<', $stat_file or die "Cannot read $stat_file: $!\n";
	my $line = <$fh> // '';
	close $fh or die "Cannot close $stat_file: $!\n";
	$line =~ s/[\r\n]+$//;
	my ($recorded_size, $recorded_mtime) = split /\s+/, $line;
	my @stat = stat($scope->{assembly_reference});
	die "Assembly reference fingerprint in $stat_file does not match $scope->{assembly_reference}; the CRAM may target another assembly generation\n"
		unless defined($recorded_size) && defined($recorded_mtime)
		&& $recorded_size =~ /^\d+$/ && $recorded_mtime =~ /^\d+$/
		&& $recorded_size == $stat[7] && $recorded_mtime == $stat[9];
}

sub scan_scope_candidates {
	my (%options) = @_;
	my $scope = $options{scope};
	my $plan = $options{plan};
	my $assembly_targets = $plan->{contig_targets}{$scope->{assembly}} || {};
	my %target_files;
	for my $target_key (sort keys %{$plan->{targets}}) {
		my $target = $plan->{targets}{$target_key};
		next unless $target->{available} && $target->{scope_keys}{$scope->{key}}
			&& $options{active_target_scopes}{$target_key}{$scope->{key}};
		my $hash = sha1_hex($scope->{key} . "\0" . $target_key);
		$target_files{$target_key} = "$options{workdir}/candidate.$hash.raw";
	}
	validate_cram_reference($scope);
	run_shell("checking CRAM $scope->{cram}", tool_command($options{samtools}, 'quickcheck', $scope->{cram}));
	my $command = tool_command(
		$options{samtools}, 'view', '-@', $options{threads},
		'-T', $scope->{assembly_reference}, '-F', 3844,
		'-q', $options{source_min_mapq}, $scope->{cram},
	);
	my $fh = open_command("scanning candidate CRAM $scope->{cram}", $command);
	open my $identity, '>', $options{source_identity_raw}
		or die "Cannot create $options{source_identity_raw}: $!\n";
	my (%handles, %last_used, $clock);
	while (my $line = <$fh>) {
		my @fields = split /\t/, $line, 12;
		die "Malformed candidate CRAM alignment in $scope->{cram}\n" unless @fields >= 11;
		my ($name, $contig) = @fields[0, 2];
		next if $fields[4] == 255;
		next if cigar_query_coverage($fields[5], "$name in $scope->{cram}") < $options{source_min_coverage};
		next unless exists $assembly_targets->{$contig};
		my @matching_targets = grep { exists $target_files{$_} }
			keys %{$assembly_targets->{$contig}};
		next unless @matching_targets;
		my $source_native = native_sequence(\@fields, "Assembly CRAM $scope->{cram}");
		print {$identity} "$name\t", sha256_hex($source_native), "\n"
			or die "Cannot write $options{source_identity_raw}: $!\n";
		for my $target_key (@matching_targets) {
			$clock++;
			if (!exists $handles{$target_key}) {
				if (keys(%handles) >= 64) {
					my ($oldest) = sort { $last_used{$a} <=> $last_used{$b} } keys %handles;
					close $handles{$oldest} or die "Cannot close candidate-name partition: $!\n";
					delete $handles{$oldest};
					delete $last_used{$oldest};
				}
				open my $out, '>>', $target_files{$target_key}
					or die "Cannot append $target_files{$target_key}: $!\n";
				$handles{$target_key} = $out;
			}
			$last_used{$target_key} = $clock;
			print {$handles{$target_key}} "$name\n"
				or die "Cannot write $target_files{$target_key}: $!\n";
		}
	}
	close_command($fh, "scanning candidate CRAM $scope->{cram}");
	close $identity or die "Cannot close $options{source_identity_raw}: $!\n";
	close $handles{$_} or die "Cannot close candidate-name partition: $!\n" for keys %handles;
	my %result;
	for my $target_key (sort keys %target_files) {
		my $raw = $target_files{$target_key};
		if (!-e $raw) {
			open my $empty, '>', $raw or die "Cannot create $raw: $!\n";
			close $empty or die "Cannot close $raw: $!\n";
		}
		my $sorted = "$raw.names";
		sort_unique($raw, $sorted, $options{workdir});
		$result{$target_key} = { file => $sorted, count => count_lines($sorted) };
	}
	return \%result;
}

sub align_and_transfer {
	my (%options) = @_;
	my $started = Time::HiRes::time();
	my $hash = sha1_hex(join("\0", $options{scope}{key}, $options{target}{key}));
	my $prefix = "$options{workdir}/align.$hash";
	my $donor_target = "$prefix.donor.name.bam";
	filter_bam_by_names(
		samtools => $options{samtools}, names => $options{names},
		input => $options{donor_union}, output => $donor_target,
		workdir => $options{workdir}, exclude_flags => 0,
	);
	my $donor_count = bam_count($options{samtools}, $donor_target);
	die "Target donor extraction returned $donor_count records for $options{candidate_count} candidate names\n"
		unless $donor_count == $options{candidate_count};

	my $acceptor_unsorted = "$prefix.acceptor.unsorted.bam";
	my $preset = $options{preset};
	my @filter = $options{scope}{technology} eq 'ONT'
		? @{$options{filter_ont}} : @{$options{filter_pb}};
	my $rg_id = file_component("$options{scope}{sample}.$options{scope}{scope}");
	my $rg = "\@RG\\tID:$rg_id\\tSM:$options{scope}{sample}\\tPL:"
		. ($options{scope}{technology} eq 'ONT' ? 'ONT' : 'PACBIO')
		. "\\tLB:meth2rep.$options{scope}{scope}";
	my $filter_report = "$prefix.bamFilter.log";
	my $mapper_report = "$prefix.minimap2.log";
	my $mapping_command = tool_command(
		$options{samtools}, 'fastq', '-n', $donor_target,
	) . ' | ' . tool_command(
		$options{minimap2}, '-2', '-a', '-Y', '-t', $options{threads}, '--secondary=no',
		'-x', $preset, '-R', $rg, $options{reference_index}, '-',
	) . ' 2> ' . shell_quote($mapper_report)
		. ' | ' . tool_command($options{bam_filter}, @filter)
		. ' 2> ' . shell_quote($filter_report)
		. ' | ' . tool_command(
			$options{samtools}, 'view', '-b', '-F', 4,
			'-o', $acceptor_unsorted, '-',
		);
	my $mapping_ok = eval {
		run_shell("fresh representative mapping for $options{target}{mgs}", $mapping_command);
		1;
	};
	if (!$mapping_ok) {
		my $reason = $@ || "representative mapping failed\n";
		for my $log ($mapper_report, $filter_report) {
			next unless -s $log;
			open my $fh, '<', $log or next;
			my @lines = <$fh>;
			close $fh;
			@lines = @lines[-10 .. -1] if @lines > 10;
			$reason .= "\n$log:\n" . join('', @lines);
		}
		die $reason;
	}
	my ($filter_stats, $filter_warnings) = read_filter_report($filter_report);
	open my $mapper_fh, '<', $mapper_report or die "Cannot read $mapper_report: $!\n";
	my @mapper_warnings = grep { /\b(?:warning|error)\b/i } <$mapper_fh>;
	close $mapper_fh or die "Cannot close $mapper_report: $!\n";
	my @warnings = (@{$filter_warnings}, @mapper_warnings);
	s/\s+\z// for @warnings;
	my $accepted_mapping_count = bam_count($options{samtools}, $acceptor_unsorted);
	if (!$accepted_mapping_count) {
		remove_temp_files($donor_target, $acceptor_unsorted, $filter_report, $mapper_report);
		return {
			aligned => 0, donor => $donor_count, filter_stats => $filter_stats,
			warnings => \@warnings,
			elapsed_seconds => Time::HiRes::time() - $started,
		};
	}

	my $acceptor_name = "$prefix.acceptor.name.bam";
	run_shell(
		"name-sorting fresh representative mappings",
		tool_command($options{samtools}, 'sort', '-n', '-@', $options{threads},
			'-m', $options{sort_memory}, '-o', $acceptor_name, $acceptor_unsorted),
	);
	my $transferred_name = "$prefix.transferred.name.bam";
	my $transfer_stats_file = "$prefix.transfer.tsv";
	my @transfer_options = (
		'--samtools', $options{samtools}, '--donor', $donor_target,
		'--acceptor', $acceptor_name, '--output', $transferred_name,
		'--stats', $transfer_stats_file,
		'--supplementary', $options{supplementary_alignments},
	);
	push @transfer_options, '--allow-missing-mn' if $options{allow_missing_mn};
	run_shell(
		"transferring MM/ML tags for full-sequence representative alignments",
		tool_command($^X, "$Bin/transfer_mod_tags.pl", @transfer_options),
	);
	my $transfer_stats = read_transfer_stats($transfer_stats_file);
	die "Tag-transfer input count changed from $accepted_mapping_count to $transfer_stats->{input_alignments} for $options{scope}{sample} $options{target}{mgs}\n"
		unless $transfer_stats->{input_alignments} == $accepted_mapping_count;
	my $transferred_count = bam_count($options{samtools}, $transferred_name);
	die "Tag-transfer BAM contains $transferred_count records but reports $transfer_stats->{output_alignments} for $options{scope}{sample} $options{target}{mgs}\n"
		unless $transferred_count == $transfer_stats->{output_alignments};
	push @warnings, "$transfer_stats->{dropped_no_primary_alignments} orphan supplementary alignment(s) were dropped after their primary failed filtering"
		if $transfer_stats->{dropped_no_primary_alignments};
	push @warnings, "$transfer_stats->{dropped_supplementary_alignments} supplementary alignment(s) were removed by the '$options{supplementary_alignments}' output policy"
		if $transfer_stats->{dropped_supplementary_alignments};
	push @warnings, "$transfer_stats->{legacy_missing_mn_reads} donor read(s) lacked MN and were accepted under --allow-missing-mn"
		if $transfer_stats->{legacy_missing_mn_reads};
	remove_temp_files($donor_target, $acceptor_unsorted, $acceptor_name,
		$filter_report, $mapper_report, $transfer_stats_file);
	if (!$transferred_count) {
		remove_temp_files($transferred_name);
		return {
			donor => $donor_count, aligned => 0,
			filter_stats => $filter_stats, transfer_stats => $transfer_stats,
			warnings => \@warnings,
			elapsed_seconds => Time::HiRes::time() - $started,
		};
	}
	return {
		donor => $donor_count, aligned => $transferred_count,
		transferred_name => $transferred_name,
		filter_stats => $filter_stats, transfer_stats => $transfer_stats,
		warnings => \@warnings,
		elapsed_seconds => Time::HiRes::time() - $started,
	};
}

my ($mgs_dir, $map_file, $mgs_report, $representatives_dir, $binner, $manifest_file, $out_dir, $out_manifest);
my ($mgs_values, $mgs_file) = ('', '');
my ($mgs2rep, $rep2rep) = (0, 0);
my ($source_min_mapq, $source_min_coverage, $threads, $memory_gb) = (10, 0.5, 4, 32);
my ($target_min_mapq, $target_min_coverage, $target_max_edit_rate, $target_min_end_clip);
my ($filter_ont, $filter_pb) = ('0.15 0.5 10 0', '0.05 0.5 30 0');
my ($preset_ont, $preset_pb) = ('map-ont', 'map-pb');
my ($supplementary_alignments, $allow_missing_mn, $output_format) = ('drop', 0, 'bam');
my ($tmp_parent, $keep_read_ids, $redo, $plan_only, $help) = ('', 0, 0, 0, 0);
my ($samtools, $minimap2, $bam_filter) = ('', '', '');

GetOptions(
	'mgs-dir=s' => \$mgs_dir,
	'map=s' => \$map_file,
	'mgs-report=s' => \$mgs_report,
	'representatives-dir=s' => \$representatives_dir,
	'binner=s' => \$binner,
	'modbam-manifest=s' => \$manifest_file,
	'out|o=s' => \$out_dir,
	'out-manifest=s' => \$out_manifest,
	'mgs=s' => \$mgs_values,
	'mgs-file=s' => \$mgs_file,
	'mgs2rep!' => \$mgs2rep,
	'rep2rep!' => \$rep2rep,
	'source-min-mapq=i' => \$source_min_mapq,
	'source-min-coverage=f' => \$source_min_coverage,
	'target-min-mapq=i' => \$target_min_mapq,
	'target-min-coverage=f' => \$target_min_coverage,
	'target-max-edit-rate=f' => \$target_max_edit_rate,
	'target-min-end-clip=i' => \$target_min_end_clip,
	'mapper-filter-ont=s' => \$filter_ont,
	'mapper-filter-pb=s' => \$filter_pb,
	'minimap2-preset-ont|mapper-preset-ont=s' => \$preset_ont,
	'minimap2-preset-pb|mapper-preset-pb=s' => \$preset_pb,
	'supplementary-alignments=s' => \$supplementary_alignments,
	'allow-missing-mn!' => \$allow_missing_mn,
	'output-format=s' => \$output_format,
	'threads=i' => \$threads,
	'memory-gb=i' => \$memory_gb,
	'tmp=s' => \$tmp_parent,
	'keep-read-ids!' => \$keep_read_ids,
	'redo!' => \$redo,
	'override!' => \$redo,
	'plan-only!' => \$plan_only,
	'samtools=s' => \$samtools,
	'minimap2=s' => \$minimap2,
	'bam-filter=s' => \$bam_filter,
	'help|h' => \$help,
) or die usage();
if ($help) {
	print usage();
	exit 0;
}
die usage() if @ARGV;
die "--mgs-dir is required\n" unless defined($mgs_dir) && $mgs_dir ne '';
$mgs_dir = File::Spec->rel2abs($mgs_dir);
$mgs_dir =~ s{/\z}{} unless $mgs_dir eq '/';
$mgs_dir = abs_path($mgs_dir) || $mgs_dir;
my $mgs_progress = require_mgs_readiness($mgs_dir);
$mgs_report ||= File::Spec->catfile($mgs_dir, 'MAGvsGC.txt.gz');
$representatives_dir ||= File::Spec->catdir($mgs_dir, 'Genomes', 'MGS_ctg');
if (-e $mgs_report) {
	my $stage1_stone = File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone');
	die "MGS membership report was modified after the Stage1 completion checkpoint: $mgs_report\n"
		if (stat($mgs_report))[9] > (stat($stage1_stone))[9];
}
if (!defined($binner) || $binner eq '') {
	($binner) = $mgs_dir =~ m{(?:^|/)Bin_([A-Za-z0-9_.-]+)/?\z};
}
$out_dir ||= File::Spec->catdir($mgs_dir, 'Meth2Rep');
for my $required (
	['--map', $map_file], ['--mgs-report', $mgs_report],
	['--representatives-dir', $representatives_dir], ['--binner', $binner],
	['--modbam-manifest', $manifest_file], ['--out', $out_dir],
) {
	die "$required->[0] is required\n" unless defined($required->[1]) && $required->[1] ne '';
}
die "Enable --mgs2rep, --rep2rep, or both\n" unless $mgs2rep || $rep2rep;
die "--threads must be a positive integer\n" unless $threads > 0;
die "--memory-gb must be a positive integer\n" unless $memory_gb > 0;
die "--source-min-mapq must be between 0 and 254 (255 means unknown)\n"
	unless $source_min_mapq >= 0 && $source_min_mapq <= 254;
die "--source-min-coverage must be between 0 and 1\n"
	unless $source_min_coverage >= 0 && $source_min_coverage <= 1;
die "--output-format must be 'bam' or 'cram'\n"
	unless $output_format eq 'bam' || $output_format eq 'cram';
die "--supplementary-alignments must be 'drop' or 'keep'\n"
	unless $supplementary_alignments eq 'drop' || $supplementary_alignments eq 'keep';
for my $preset_spec (['--minimap2-preset-ont', $preset_ont], ['--minimap2-preset-pb', $preset_pb]) {
	die "$preset_spec->[0] contains unsafe or unsupported characters\n"
		unless $preset_spec->[1] =~ /\A[A-Za-z0-9_.:+-]+\z/;
}
my @filter_ont = validate_filter($filter_ont, '--mapper-filter-ont');
my @filter_pb = validate_filter($filter_pb, '--mapper-filter-pb');
if (defined $target_min_mapq) {
	die "--target-min-mapq must be between 0 and 254\n"
		unless $target_min_mapq >= 0 && $target_min_mapq <= 254;
	$filter_ont[2] = $filter_pb[2] = $target_min_mapq;
}
if (defined $target_min_coverage) {
	die "--target-min-coverage must be between 0 and 1\n"
		unless $target_min_coverage >= 0 && $target_min_coverage <= 1;
	$filter_ont[1] = $filter_pb[1] = $target_min_coverage;
}
if (defined $target_max_edit_rate) {
	die "--target-max-edit-rate must be between 0 and 1\n"
		unless $target_max_edit_rate >= 0 && $target_max_edit_rate <= 1;
	$filter_ont[0] = $filter_pb[0] = $target_max_edit_rate;
}
if (defined $target_min_end_clip) {
	die "--target-min-end-clip must be a nonnegative integer\n"
		unless $target_min_end_clip >= 0;
	$filter_ont[3] = $filter_pb[3] = $target_min_end_clip;
}
my $sort_memory_mb = int(($memory_gb * 1024 * 0.5) / $threads);
$sort_memory_mb = 1 if $sort_memory_mb < 1;
$sort_memory_mb = 1024 if $sort_memory_mb > 1024;
my $sort_memory = $sort_memory_mb . 'M';

my $selected = read_target_mgs(mgs => $mgs_values, mgs_file => $mgs_file);
my $manifest = read_modbam_manifest($manifest_file);
my ($assembly_groups, $map) = getDirsPerAssmblGrp($map_file);
my @modes = grep { $_->[1] } (['mgs2rep', $mgs2rep], ['rep2rep', $rep2rep]);
@modes = map { $_->[0] } @modes;
my $plan = build_meth2rep_plan(
	map => $map, assembly_groups => $assembly_groups, manifest => $manifest,
	manifest_file => File::Spec->rel2abs($manifest_file), selected => $selected,
	mgs_report => File::Spec->rel2abs($mgs_report),
	representatives_dir => File::Spec->rel2abs($representatives_dir),
	binner => $binner, modes => \@modes,
);
my %checkpoint_representatives = map {
	my $path = $_->{path} // '';
	my $canonical = $path ne '' ? abs_path($path) : undef;
	defined($canonical) ? ($canonical => 1) : ()
}
	@{$mgs_progress->{BinExtr}{outputs}};
for my $target (values %{$plan->{targets}}) {
	next unless $target->{available};
	die "Selected representative $target->{reference} is not recorded by the completed MGS bin-extraction checkpoint\n"
		unless $checkpoint_representatives{$target->{reference}};
}
my %encoded_samples;
for my $scope (values %{$plan->{scopes}}) {
	my $encoded = file_component($scope->{sample});
	die "Sample names '$encoded_samples{$encoded}' and '$scope->{sample}' collide in meth2rep output filenames\n"
		if exists($encoded_samples{$encoded}) && $encoded_samples{$encoded} ne $scope->{sample};
	$encoded_samples{$encoded} = $scope->{sample};
}
my %plan_inputs = map { $_ => 1 } @{$plan->{inputs}};
$plan_inputs{File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone')} = 1;
$plan_inputs{File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone')} = 1;
for my $map_path (split /,/, $map_file) {
	$map_path = File::Spec->rel2abs($map_path);
	die "Resolved mapping input is missing: $map_path\n" unless -e $map_path;
	$plan_inputs{$map_path} = 1;
}
if ($mgs_file ne '') {
	my $selection_path = File::Spec->rel2abs($mgs_file);
	$plan_inputs{$selection_path} = 1;
}
$plan->{inputs} = [sort keys %plan_inputs];

$out_dir = File::Spec->rel2abs($out_dir);
$out_dir =~ s{/\z}{} unless $out_dir eq '/';
$out_dir = abs_path($out_dir) if -d $out_dir;
die "--out cannot be the MGS directory or one of its ancestors\n"
	if index("$mgs_dir/", "$out_dir/") == 0;
for my $input (@{$plan->{inputs}}) {
	die "meth2rep output directory contains a required input ($input); choose an isolated --out\n"
		if index($input, "$out_dir/") == 0;
}
make_path($out_dir);
$out_dir = abs_path($out_dir) || $out_dir;
die "--out resolves to the MGS directory or one of its ancestors\n"
	if index("$mgs_dir/", "$out_dir/") == 0;
for my $input (@{$plan->{inputs}}) {
	my $canonical_input = abs_path($input) || $input;
	die "meth2rep output directory contains a required input ($canonical_input); choose an isolated --out\n"
		if index($canonical_input, "$out_dir/") == 0;
}
my $state_dir = File::Spec->catdir($out_dir, '.meth2rep');
die "Legacy meth2rep mode-first directories are present in $out_dir; choose a fresh -o directory to avoid mixing output layouts\n"
	if -d File::Spec->catdir($out_dir, 'mgs2rep')
	|| -d File::Spec->catdir($out_dir, 'rep2rep');
make_path($state_dir);
$out_manifest ||= File::Spec->catfile($out_dir, 'manifest.tsv');
$out_manifest = File::Spec->rel2abs($out_manifest);
die "--out-manifest collides with a required input file: $out_manifest\n"
	if grep { $_ eq $out_manifest } @{$plan->{inputs}};
die "--out-manifest collides with a meth2rep control file: $out_manifest\n"
	if grep { $out_manifest eq File::Spec->catfile($state_dir, $_) }
		qw(plan.tsv plan.preview.tsv summary.tsv provenance.json complete.stone .meth2rep.lock);
die "--out-manifest cannot be placed inside a meth2rep data subdirectory: $out_manifest\n"
	if grep { index($out_manifest, File::Spec->catdir($out_dir, $_) . '/') == 0 }
		qw(.meth2rep read_ids);
die "--out-manifest cannot be placed inside a per-MGS output directory or another nested meth2rep data directory: $out_manifest\n"
	if index($out_manifest, "$out_dir/") == 0
	&& dirname($out_manifest) ne $out_dir;
die "--out-manifest names an existing directory rather than a TSV file: $out_manifest\n"
	if -d $out_manifest;
make_path(dirname($out_manifest)) unless -d dirname($out_manifest);
make_path(File::Spec->catdir($state_dir, 'units'));
my $lock_path = File::Spec->catfile($state_dir, '.meth2rep.lock');
open my $lock_fh, '>>', $lock_path or die "Cannot open meth2rep lock $lock_path: $!\n";
flock($lock_fh, LOCK_EX | LOCK_NB) or die "Another meth2rep run is using $out_dir\n";
my $plan_file = File::Spec->catfile($state_dir, $plan_only ? 'plan.preview.tsv' : 'plan.tsv');
my $summary_file = File::Spec->catfile($state_dir, 'summary.tsv');
my $provenance_file = File::Spec->catfile($state_dir, 'provenance.json');
my $checkpoint_file = File::Spec->catfile($state_dir, 'complete.stone');
my (%tool_versions, %tool_identities);
my $has_runnable_target = grep { $_->{available} } values %{$plan->{targets}};
if ($has_runnable_target) {
	$samtools ||= find_on_path('samtools')
		or die "samtools was not found on PATH; provide --samtools COMMAND\n";
	$minimap2 ||= find_on_path('minimap2')
		or die "minimap2 was not found on PATH; provide --minimap2 COMMAND\n";
	$bam_filter ||= join(' ', shell_quote($^X),
		shell_quote(File::Spec->catfile($Bin, '..', 'assemblies', 'bamFilter.pl')));
	$samtools = canonical_command($samtools, 'samtools');
	$minimap2 = canonical_command($minimap2, 'minimap2');
	$bam_filter = canonical_command($bam_filter, 'bamFilter');
	for my $tool_spec (
		['samtools', $samtools, '--version'], ['minimap2', $minimap2, '--version'],
	) {
		my $version = capture_command(
			"checking $tool_spec->[0]", tool_command($tool_spec->[1], $tool_spec->[2]),
		);
		($version) = split /\n/, $version, 2;
		$tool_versions{$tool_spec->[0]} = $version;
	}
	my %required_presets = map {
		my $technology = $_->{technology};
		($technology => ($technology eq 'ONT' ? $preset_ont : $preset_pb));
	} values %{$plan->{scopes}};
	for my $technology (sort keys %required_presets) {
		capture_command(
			"validating minimap2 preset $required_presets{$technology} for $technology",
			tool_command($minimap2, '-x', $required_presets{$technology}, '--version'),
		);
	}
	$tool_identities{samtools} = command_identity($samtools);
	$tool_identities{minimap2} = command_identity($minimap2);
	$tool_identities{bam_filter} = command_identity($bam_filter);
} else {
	($samtools, $minimap2, $bam_filter) = ('not_used_this_invocation') x 3;
	%tool_versions = map { $_ => 'not_used_this_invocation' } qw(samtools minimap2);
	%tool_identities = map { $_ => 'not_used_this_invocation' } qw(samtools minimap2 bam_filter);
}
my %checkpoint_parameters = (
	component_version => $VERSION,
	modes => join(',', @modes),
	mgs => join(',', sort keys %{$selected}),
	binner => $binner,
	source_min_mapq => $source_min_mapq,
	source_min_coverage => $source_min_coverage,
	mapper_filter_ont => join(' ', @filter_ont),
	mapper_filter_pb => join(' ', @filter_pb),
	minimap2_preset_ont => $preset_ont,
	minimap2_preset_pb => $preset_pb,
	supplementary_alignments => $supplementary_alignments,
	allow_missing_mn => $allow_missing_mn ? 1 : 0,
	output_format => $output_format,
	samtools_identity => $tool_identities{samtools},
	minimap2_identity => $tool_identities{minimap2},
	bam_filter_identity => $tool_identities{bam_filter},
	samtools_version => $tool_versions{samtools},
	minimap2_version => $tool_versions{minimap2},
	keep_read_ids => $keep_read_ids ? 1 : 0,
);
$checkpoint_parameters{input_fingerprint} = fingerprint_inputs($plan->{inputs}, \%checkpoint_parameters);

my (%active_target_scopes, %cached_units, %unit_context);
for my $target_key (sort keys %{$plan->{targets}}) {
	my $target = $plan->{targets}{$target_key};
	next unless $target->{available};
	my %sample_scopes;
	for my $scope_key (keys %{$target->{scope_keys}}) {
		my $sample = $plan->{scopes}{$scope_key}{sample};
		push @{$sample_scopes{$sample}}, $scope_key;
	}
	for my $sample (sort keys %sample_scopes) {
		my @scope_keys = sort @{$sample_scopes{$sample}};
		my %unit_groups = map { $plan->{scopes}{$_}{group} => 1 } @scope_keys;
		my %files = map { $_ => 1 } (
			File::Spec->rel2abs(__FILE__),
			File::Spec->catfile($Bin, 'transfer_mod_tags.pl'),
			File::Spec->catfile($Bin, '..', 'assemblies', 'bamFilter.pl'),
			File::Spec->catfile($Bin, '..', '..', 'Mods', 'Meth2Rep.pm'),
			File::Spec->catfile($Bin, '..', '..', 'Mods', 'Checkpoint.pm'),
			File::Spec->rel2abs($mgs_report),
			File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone'),
			File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone'),
			$target->{reference},
			(map { $_->{assignment} }
				grep { $unit_groups{$_->{group}} } @{$target->{sources}}),
			(map { File::Spec->rel2abs($_) } split /,/, $map_file),
		);
		for my $scope_key (@scope_keys) {
			my $scope = $plan->{scopes}{$scope_key};
			$files{$_} = 1 for ($scope->{cram}, "$scope->{cram}.sto",
				$scope->{assembly_reference}, @{$scope->{modbams}});
			(my $reference_stat = $scope->{cram}) =~ s/\.cram\z/.reference.stat/;
			$files{$reference_stat} = 1 if -e $reference_stat;
		}
		my %unit_parameters = (
			component_version => $VERSION, mode => $target->{mode},
			mgs => $target->{mgs}, sample => $sample,
			binner => $binner,
			source_min_mapq => $source_min_mapq,
			source_min_coverage => $source_min_coverage,
			mapper_filter_ont => join(' ', @filter_ont),
			mapper_filter_pb => join(' ', @filter_pb),
			minimap2_preset_ont => $preset_ont,
			minimap2_preset_pb => $preset_pb,
			supplementary_alignments => $supplementary_alignments,
			allow_missing_mn => $allow_missing_mn ? 1 : 0,
			output_format => $output_format,
			samtools_identity => $tool_identities{samtools},
			minimap2_identity => $tool_identities{minimap2},
			bam_filter_identity => $tool_identities{bam_filter},
			samtools_version => $tool_versions{samtools},
			minimap2_version => $tool_versions{minimap2},
			keep_read_ids => $keep_read_ids ? 1 : 0,
		);
		$unit_parameters{input_fingerprint} = fingerprint_inputs([sort keys %files], \%unit_parameters);
		my ($json, $stone) = unit_paths($state_dir, $target, $sample);
		my $unit_key = "$target_key\t$sample";
		my $expected_alignment = output_alignment_path($out_dir, $target, $sample, $output_format);
		my $expected_index = alignment_index_path($expected_alignment, $output_format);
		my $alternate_format = $output_format eq 'bam' ? 'cram' : 'bam';
		my $alternate_alignment = output_alignment_path($out_dir, $target, $sample, $alternate_format);
		my $alternate_index = alignment_index_path($alternate_alignment, $alternate_format);
		die "Untracked output alignment or index already exists: $expected_alignment; inspect it and use --override only if replacement is intended\n"
			if (-e $expected_alignment || -e $expected_index) && !-e $stone && !$redo;
		die "An alternate-format meth2rep output already exists ($alternate_alignment); use --override to replace it without leaving ambiguous BAM/CRAM siblings\n"
			if (-e $alternate_alignment || -e $alternate_index) && !$redo;
		$unit_context{$unit_key} = {
			json => $json, stone => $stone, parameters => \%unit_parameters,
			scope_keys => \@scope_keys,
			expected_alignment => $expected_alignment,
			expected_index => $expected_index,
			alternate_alignment => $alternate_alignment,
			alternate_index => $alternate_index,
		};
		if (!$redo && -s $stone && checkpoint_valid($stone, parameters => \%unit_parameters) && -s $json) {
			my $record = read_json_file($json);
			die "Unit checkpoint metadata mismatch: $json\n"
				unless ($record->{mgs} // '') eq $target->{mgs}
				&& ($record->{mode} // '') eq $target->{mode}
				&& ($record->{sample} // '') eq $sample;
			$cached_units{$unit_key} = $record;
			next;
		}
		$active_target_scopes{$target_key}{$_} = 1 for @scope_keys;
	}
}

if ($plan_only) {
	write_plan($plan_file, $plan);
	print "meth2rep v$VERSION: input, tool, cache, and output preflight passed; preview written to $plan_file\n";
	exit 0;
}
unlink $checkpoint_file or die "Cannot invalidate stale meth2rep checkpoint $checkpoint_file: $!\n"
	if -e $checkpoint_file;
write_plan($plan_file, $plan);

my $workdir = '';
if (keys %active_target_scopes) {
	$tmp_parent = $ENV{SLURM_TMPDIR}
		if $tmp_parent eq '' && defined($ENV{SLURM_TMPDIR}) && -d $ENV{SLURM_TMPDIR};
	$tmp_parent = File::Spec->rel2abs($tmp_parent) if $tmp_parent ne '';
	make_path($tmp_parent) if $tmp_parent ne '' && !-d $tmp_parent;
	my %temp_options = (CLEANUP => 1);
	$temp_options{DIR} = $tmp_parent if $tmp_parent ne '';
	$workdir = tempdir('meth2rep.XXXXXX', %temp_options);
}

print "meth2rep v$VERSION: " . scalar(keys %{$selected}) . " selected MGS; modes "
	. join(', ', @modes) . "; " . scalar(keys %{$plan->{scopes}}) . " eligible sample/scope CRAM(s)\n";
warn "Manifest row is not required by the selected MGS plan: $_\n"
	for @{$plan->{manifest_unused}};

my (%candidate, %scope_results, %reference_indexes, %nonempty_targets, %donor_scopes);
my (%source_identity, %cram_reference_cache);
my (%source_scan_seconds, %donor_stream_seconds);
for my $scope_key (sort keys %{$plan->{scopes}}) {
	my $scope = $plan->{scopes}{$scope_key};
	next unless grep { $active_target_scopes{$_}{$scope_key} } keys %active_target_scopes;
	my $scan_started = Time::HiRes::time();
	my $scope_hash = sha1_hex($scope_key);
	my $source_identity_raw = "$workdir/source.$scope_hash.identity.raw";
	$candidate{$scope_key} = scan_scope_candidates(
		scope => $scope, plan => $plan, workdir => $workdir,
		samtools => $samtools, threads => $threads, source_min_mapq => $source_min_mapq,
		source_min_coverage => $source_min_coverage,
		active_target_scopes => \%active_target_scopes,
		source_identity_raw => $source_identity_raw,
	);
	my $source_identity_sorted = "$source_identity_raw.sorted";
	sort_identity_records($source_identity_raw, $source_identity_sorted, $workdir);
	$source_identity{$scope_key} = $source_identity_sorted;
	remove_temp_files($source_identity_raw);
	$source_scan_seconds{$scope_key} = Time::HiRes::time() - $scan_started;
	my @nonempty_targets = grep { $candidate{$scope_key}{$_}{count} > 0 }
		sort keys %{$candidate{$scope_key}};
	next unless @nonempty_targets;
	$nonempty_targets{$scope_key} = \@nonempty_targets;
	push @{$donor_scopes{$_}}, $scope_key for @{$scope->{modbams}};
}

# Keep candidate membership scope-specific, but defer donor access until every
# CRAM has been scanned. This guarantees one record-streaming pass through each
# distinct physical modBAM even when a manifest explicitly reuses it.
my (%donor_name_bam, %scope_donor_names);
for my $modbam (sort keys %donor_scopes) {
	my $donor_started = Time::HiRes::time();
	my @scope_keys = sort @{$donor_scopes{$modbam}};
	my %technologies = map { $plan->{scopes}{$_}{technology} => 1 } @scope_keys;
	die "Original modBAM '$modbam' is declared for more than one sequencing technology ("
		. join(', ', sort keys %technologies) . "); split it into unambiguous manifest donors\n"
		if keys(%technologies) > 1;
	my $donor_hash = sha1_hex($modbam);
	my $union_raw = "$workdir/donor.$donor_hash.union.raw";
	concatenate_files(
		$union_raw,
		map {
			my $scope_key = $_;
			map { $candidate{$scope_key}{$_}{file} } @{$nonempty_targets{$scope_key}}
		} @scope_keys,
	);
	my $union_names = "$workdir/donor.$donor_hash.union.names";
	sort_unique($union_raw, $union_names, $workdir);
	# Raw, unmapped modBAMs legitimately have no @SQ records, which makes older
	# samtools quickcheck versions reject them.  Parsing the full selected stream
	# below catches truncation; this header read supplies the early format check.
	run_shell("checking donor modBAM header $modbam",
		tool_command($samtools, 'view', '-H', $modbam) . ' >/dev/null');
	my $donor_unsorted = "$workdir/donor.$donor_hash.primary.bam";
	filter_bam_by_names(
		samtools => $samtools, names => $union_names,
		input => $modbam, output => $donor_unsorted,
		workdir => $workdir, exclude_flags => 3840,
	);
	my $donor_name = "$workdir/donor.$donor_hash.name.bam";
	run_shell(
		"name-sorting donor subset $modbam",
		tool_command($samtools, 'sort', '-n', '-@', $threads,
			'-m', $sort_memory, '-o', $donor_name, $donor_unsorted),
	);
	my $found_raw = "$workdir/donor.$donor_hash.found.raw";
	my ($donor_count) = validate_donor_subset(
		$samtools, $donor_name, $found_raw, undef, $allow_missing_mn,
	);
	my $found_names = "$workdir/donor.$donor_hash.found.names";
	sort_unique($found_raw, $found_names, $workdir);
	$donor_name_bam{$modbam} = $donor_name;
	$donor_stream_seconds{$modbam} = Time::HiRes::time() - $donor_started;
	remove_temp_files($union_raw, $union_names, $donor_unsorted, $found_raw, $found_names);
}

for my $scope_key (sort keys %nonempty_targets) {
		my $scope = $plan->{scopes}{$scope_key};
		my $scope_hash = sha1_hex($scope_key);
		my $scope_union_raw = "$workdir/scope.$scope_hash.union.raw";
		concatenate_files($scope_union_raw,
			map { $candidate{$scope_key}{$_}{file} } @{$nonempty_targets{$scope_key}});
		my $scope_union_names = "$workdir/scope.$scope_hash.union.names";
		sort_unique($scope_union_raw, $scope_union_names, $workdir);
		my @scope_donor_bams;
		for my $modbam (@{$scope->{modbams}}) {
			my $donor_hash = sha1_hex("$scope_key\0$modbam");
			my $scope_donor_unsorted = "$workdir/scope.donor.$donor_hash.unsorted.bam";
			filter_bam_by_names(
				samtools => $samtools, names => $scope_union_names,
				input => $donor_name_bam{$modbam}, output => $scope_donor_unsorted,
				workdir => $workdir, exclude_flags => 0,
			);
			my $scope_donor_bam = "$workdir/scope.donor.$donor_hash.name.bam";
			run_shell("name-sorting scope donor subset",
				tool_command($samtools, 'sort', '-n', '-@', $threads,
					'-m', $sort_memory, '-o', $scope_donor_bam, $scope_donor_unsorted));
			my $raw_names = "$workdir/scope.donor.$donor_hash.names.raw";
			validate_donor_subset(
				$samtools, $scope_donor_bam, $raw_names, undef, $allow_missing_mn,
			);
			my $sorted_names = "$raw_names.sorted";
			sort_unique($raw_names, $sorted_names, $workdir);
			$scope_donor_names{$scope_key}{$modbam} = $sorted_names;
			push @scope_donor_bams, $scope_donor_bam;
			remove_temp_files($scope_donor_unsorted, $raw_names);
		}
		my $scope_donor_union = "$workdir/scope.$scope_hash.donors.name.bam";
		if (@scope_donor_bams == 1) {
			$scope_donor_union = $scope_donor_bams[0];
		} else {
			run_shell("merging selected original modBAM donors for $scope->{sample}:$scope->{scope}",
				tool_command($samtools, 'merge', '-n', '-@', $threads, '-f',
					$scope_donor_union, @scope_donor_bams));
		}
		my $combined_names_raw = "$workdir/scope.$scope_hash.found.raw";
		my $combined_identity_raw = "$workdir/scope.$scope_hash.identity.raw";
		my ($combined_count) = validate_donor_subset(
			$samtools, $scope_donor_union, $combined_names_raw,
			$combined_identity_raw, $allow_missing_mn,
		);
		my $combined_names = "$combined_names_raw.sorted";
		sort_unique($combined_names_raw, $combined_names, $workdir);
		my $combined_identity = "$combined_identity_raw.sorted";
		sort_identity_records($combined_identity_raw, $combined_identity, $workdir);
		assert_name_subset($scope_union_names, $combined_names,
			"sample '$scope->{sample}' scope '$scope->{scope}' across all declared modBAM donors");
		die "Donor identity collision in sample '$scope->{sample}' scope '$scope->{scope}'\n"
			unless $combined_count == count_lines($scope_union_names);
		assert_source_donor_identity(
			$source_identity{$scope_key}, $combined_identity,
			"sample '$scope->{sample}' scope '$scope->{scope}'",
		);
		remove_temp_files($scope_union_raw, $scope_union_names,
			$combined_names_raw, $combined_names,
			$combined_identity_raw, $combined_identity,
			$source_identity{$scope_key});
		my %mgs_for_scope;
		for my $target_key (@{$nonempty_targets{$scope_key}}) {
			my $target = $plan->{targets}{$target_key};
			$mgs_for_scope{$target->{mgs}}{$target->{mode}} = $target_key;
		}
		for my $mgs (sort keys %mgs_for_scope) {
			my $by_mode = $mgs_for_scope{$mgs};
			my $primary_key = exists($by_mode->{mgs2rep}) ? $by_mode->{mgs2rep} : $by_mode->{rep2rep};
			my $primary_target = $plan->{targets}{$primary_key};
			my $reference = $primary_target->{reference};
			my $preset = $scope->{technology} eq 'ONT' ? $preset_ont : $preset_pb;
			my $index_key = join("\0", $reference, $preset);
			if (!exists $reference_indexes{$index_key}) {
				my $index = "$workdir/reference." . sha1_hex($index_key) . '.mmi';
				run_shell(
					"indexing representative $mgs with preset $preset",
					tool_command($minimap2, '-x', $preset, '-t', $threads,
						'-d', $index, $reference),
				);
				$reference_indexes{$index_key} = $index;
			}
			my $primary_result = align_and_transfer(
				scope => $scope, target => $primary_target,
				names => $candidate{$scope_key}{$primary_key}{file},
				candidate_count => $candidate{$scope_key}{$primary_key}{count},
				donor_union => $scope_donor_union, reference_index => $reference_indexes{$index_key},
				samtools => $samtools, minimap2 => $minimap2,
				bam_filter => $bam_filter, filter_ont => \@filter_ont, filter_pb => \@filter_pb,
				preset => $preset, supplementary_alignments => $supplementary_alignments,
				allow_missing_mn => $allow_missing_mn,
				threads => $threads, sort_memory => $sort_memory, workdir => $workdir,
			);
			$scope_results{$primary_key}{$scope_key} = $primary_result;

			if (exists($by_mode->{rep2rep}) && $primary_key ne $by_mode->{rep2rep}) {
				my $rep_key = $by_mode->{rep2rep};
				my $rep_candidate = $candidate{$scope_key}{$rep_key};
				if (!$primary_result->{aligned}) {
					$scope_results{$rep_key}{$scope_key} = {
						donor => $rep_candidate->{count}, aligned => 0,
						reused_from => $primary_key, elapsed_seconds => 0,
					};
				} else {
					my $subset_started = Time::HiRes::time();
					my $subset = "$workdir/transferred." . sha1_hex("$scope_key\0$rep_key") . '.name.bam';
					filter_bam_by_names(
						samtools => $samtools, names => $rep_candidate->{file},
						input => $primary_result->{transferred_name}, output => $subset,
						workdir => $workdir, exclude_flags => 0,
					);
					my $subset_count = bam_count($samtools, $subset);
					$scope_results{$rep_key}{$scope_key} = {
						donor => $rep_candidate->{count}, aligned => $subset_count,
						($subset_count ? (transferred_name => $subset) : ()),
						reused_from => $primary_key,
						elapsed_seconds => Time::HiRes::time() - $subset_started,
					};
					remove_temp_files($subset) unless $subset_count;
				}
			}
			for my $completed_key (values %{$by_mode}) {
				my $result = $scope_results{$completed_key}{$scope_key};
				next unless $result && $result->{aligned};
				my $coordinate = "$workdir/coordinate."
					. sha1_hex("$scope_key\0$completed_key") . '.bam';
				run_shell(
					"coordinate-sorting transferred alignment",
					tool_command($samtools, 'sort', '-@', $threads,
						'-m', $sort_memory, '-o', $coordinate, $result->{transferred_name}),
				);
				$result->{coordinate_bam} = $coordinate;
				remove_temp_files($result->{transferred_name});
				delete $result->{transferred_name};
			}
		}
		my %scope_bams = map { $_ => 1 } (@scope_donor_bams, $scope_donor_union);
		remove_temp_files(sort keys %scope_bams);
}
remove_temp_files(values %donor_name_bam);

my @summary_rows;
push @summary_rows, {
	mgs => $_->{mgs}, mode => $_->{mode}, sample => '-', scopes => '-',
	status => $_->{reason}, candidates => 0, donor => 0, aligned => 0,
	mm => 0, ml => 0, alignment => '-', alignment_format => '-', index => '-',
} for @{$plan->{unavailable}};
my @checkpoint_outputs = ($plan_file);
my %reference_length_cache;
for my $target_key (sort keys %{$plan->{targets}}) {
	my $target = $plan->{targets}{$target_key};
	next unless $target->{available};
	my %sample_scopes;
	for my $scope_key (sort keys %{$target->{scope_keys}}) {
		my $scope = $plan->{scopes}{$scope_key};
		$sample_scopes{$scope->{sample}}{$scope->{scope}} = $scope_key;
	}
	for my $sample (sort keys %sample_scopes) {
		my $unit_key = "$target_key\t$sample";
		if (exists $cached_units{$unit_key}) {
			push @summary_rows, $cached_units{$unit_key}{summary};
			push @checkpoint_outputs, $unit_context{$unit_key}{json};
			my $cached_alignment = $cached_units{$unit_key}{summary}{alignment};
			my $cached_index = $cached_units{$unit_key}{summary}{index};
			push @checkpoint_outputs, $cached_alignment, $cached_index
				if defined($cached_alignment) && $cached_alignment ne '-';
			next;
		}
		my (@coordinate_bams, @scope_names, @name_files, @scope_reports, @unit_warnings);
		my ($candidate_count, $donor_count, $aligned_count, $processing_seconds) = (0, 0, 0, 0);
		my @unit_donors;
		for my $scope_name (sort keys %{$sample_scopes{$sample}}) {
			my $scope_key = $sample_scopes{$sample}{$scope_name};
			push @scope_names, $scope_name;
			my $candidate_record = $candidate{$scope_key}{$target_key};
			my $count = $candidate_record ? $candidate_record->{count} : 0;
			$candidate_count += $count;
			push @name_files, [$scope_name, $candidate_record->{file}] if $candidate_record && $count;
			my $scope_report = {
				scope => $scope_name,
				technology => $plan->{scopes}{$scope_key}{technology},
				assembly_cram => $plan->{scopes}{$scope_key}{cram},
				candidate_reads => $count,
				donor_reads => 0, accepted_alignment_records => 0,
				minimap2_preset => $plan->{scopes}{$scope_key}{technology} eq 'ONT'
					? $preset_ont : $preset_pb,
			};
			push @scope_reports, $scope_report;
			for my $modbam (@{$plan->{scopes}{$scope_key}{modbams}}) {
				my $donor_names = $scope_donor_names{$scope_key}{$modbam};
				my $donor_candidates = $count && $donor_names
					? intersection_count($candidate_record->{file}, $donor_names) : 0;
				push @unit_donors, {
					scope => $scope_name,
					technology => $plan->{scopes}{$scope_key}{technology},
					modbam => $modbam,
					candidate_reads => $donor_candidates,
				};
			}
			next unless $count;
			my $result = $scope_results{$target_key}{$scope_key}
				or die "Internal error: missing alignment result for $scope_key / $target_key\n";
			$scope_report->{donor_reads} = $result->{donor};
			$scope_report->{accepted_alignment_records} = $result->{aligned};
			$scope_report->{filter_stats} = $result->{filter_stats} if $result->{filter_stats};
			$scope_report->{transfer_stats} = $result->{transfer_stats} if $result->{transfer_stats};
			$scope_report->{reused_mapping_from} = 'mgs2rep' if $result->{reused_from};
			push @unit_warnings, @{$result->{warnings} || []};
			$processing_seconds += $result->{elapsed_seconds} || 0;
			$donor_count += $result->{donor};
			$aligned_count += $result->{aligned};
			next unless $result->{aligned};
			die "Internal error: missing coordinate-sorted transfer for $scope_key / $target_key\n"
				unless defined($result->{coordinate_bam}) && -s $result->{coordinate_bam};
			push @coordinate_bams, $result->{coordinate_bam};
		}
		my $attributed = 0;
		$attributed += $_->{candidate_reads} for @unit_donors;
		die "Donor attribution for $target->{mgs} $target->{mode} $sample accounts for $attributed of $candidate_count candidates\n"
			unless $attributed == $candidate_count;
		for my $left_index (0 .. $#name_files) {
			for my $right_index ($left_index + 1 .. $#name_files) {
				my $shared = intersection_count($name_files[$left_index][1], $name_files[$right_index][1]);
				die "Read identity collision for $target->{mgs} $target->{mode} $sample: $shared QNAME(s) appear in both $name_files[$left_index][0] and $name_files[$right_index][0] scopes\n"
					if $shared;
			}
		}

		my $target_dir = File::Spec->catdir($out_dir, $target->{mgs});
		my $unit = $unit_context{$unit_key};
		my $expected_alignment = $unit->{expected_alignment};
		my $expected_index = $unit->{expected_index};
		my ($status, $final_alignment, $final_index, $mm_count, $ml_count) =
			('no_candidates', '-', '-', 0, 0);
		my $lengths = $reference_length_cache{$target->{reference}}
			||= reference_lengths($target->{reference});
		my $reference_bases = 0;
		$reference_bases += $_ for values %{$lengths};
		my $coverage = {
			reference_bases => $reference_bases, covered_bases => 0,
			breadth_fraction => 0, mean_depth => 0, depth_sum => 0,
			method => 'no accepted representative alignments',
		};
		if ($candidate_count && !$aligned_count) {
			$status = 'no_target_alignment';
			push @unit_warnings, 'Candidate reads passed source filtering, but none passed representative-alignment filtering';
		} elsif ($aligned_count) {
			$status = 'complete';
			make_path($target_dir);
			$final_alignment = $expected_alignment;
			$final_index = $expected_index;
			my $partial_alignment = "$final_alignment.part.$$";
			my $partial_index = alignment_index_path($partial_alignment, $output_format);
			push @publication_partials, $partial_alignment, $partial_index;
			my $merged_bam = $output_format eq 'bam'
				? $partial_alignment
				: "$workdir/final." . sha1_hex("$target_key\0$sample") . '.bam';
			run_shell(
				"merging transferred sample scopes",
				tool_command($samtools, 'merge', '-@', $threads, '-f', $merged_bam, @coordinate_bams),
			);
			if ($output_format eq 'cram') {
				my $cram_reference = materialize_cram_reference(
					$target->{reference}, $workdir, $samtools, \%cram_reference_cache,
				);
				run_shell(
					"encoding self-contained final modCRAM",
					tool_command($samtools, 'view', '-@', $threads, '-C',
						'-T', $cram_reference, '--output-fmt-option', 'embed_ref=1',
						'-o', $partial_alignment, $merged_bam),
				);
			}
			run_shell("checking final mod$output_format",
				tool_command($samtools, 'quickcheck', $partial_alignment));
			my ($published_count, $published_mm, $published_ml) =
				bam_tag_counts($samtools, $partial_alignment);
			die "Published mod$output_format record count changed from $aligned_count to $published_count for $target->{mgs} $sample\n"
				unless $published_count == $aligned_count;
			($mm_count, $ml_count) = ($published_mm, $published_ml);
			run_shell("indexing final mod$output_format",
				tool_command($samtools, 'index', $partial_alignment, $partial_index));
			run_shell("validating final alignment index",
				tool_command($samtools, 'idxstats', $partial_alignment) . ' >/dev/null');
			$coverage = bam_coverage($samtools, $partial_alignment, $lengths);
			rename $partial_alignment, $final_alignment
				or die "Cannot publish $final_alignment: $!\n";
			rename $partial_index, $final_index
				or die "Cannot publish $final_index: $!\n";
			if ($redo) {
				unlink $unit->{alternate_alignment}
					or die "Cannot remove replaced alternate alignment $unit->{alternate_alignment}: $!\n"
					if -e $unit->{alternate_alignment};
				unlink $unit->{alternate_index}
					or die "Cannot remove replaced alternate index $unit->{alternate_index}: $!\n"
					if -e $unit->{alternate_index};
			}
			push @checkpoint_outputs, $final_alignment, $final_index;
		} else {
			push @unit_warnings, 'No source-MAG candidate reads passed the source filters';
			unlink $expected_alignment
				or die "Cannot remove stale derived alignment $expected_alignment: $!\n"
				if -e $expected_alignment;
			unlink $expected_index
				or die "Cannot remove stale derived index $expected_index: $!\n"
				if -e $expected_index;
			if ($redo) {
				unlink $unit->{alternate_alignment}
					or die "Cannot remove replaced alternate alignment $unit->{alternate_alignment}: $!\n"
					if -e $unit->{alternate_alignment};
				unlink $unit->{alternate_index}
					or die "Cannot remove replaced alternate index $unit->{alternate_index}: $!\n"
					if -e $unit->{alternate_index};
			}
		}
		remove_temp_files(@coordinate_bams);

		my $published_origin_file = '';
		if ($keep_read_ids && @name_files) {
			my $id_dir = File::Spec->catdir($state_dir, 'read_ids', $target->{mode}, $target->{mgs});
			make_path($id_dir);
			my $id_file = File::Spec->catfile($id_dir, file_component($sample) . '.read_origins.tsv.gz');
			my $id_partial = "$id_file.part.$$";
			push @publication_partials, $id_partial;
			my $gzip_fh = IO::Compress::Gzip->new($id_partial)
				or die "Cannot create $id_partial: $GzipError\n";
			print {$gzip_fh} "sample\tscope\tqname\toriginal_modbam\n"
				or die "Cannot write $id_partial: $GzipError\n";
			for my $spec (@name_files) {
				my $scope_key = $sample . "\t" . $spec->[0];
				for my $modbam (@{$plan->{scopes}{$scope_key}{modbams}}) {
					my $donor_names = $scope_donor_names{$scope_key}{$modbam};
					next unless $donor_names;
					write_intersection_origin_rows($spec->[1], $donor_names,
						$gzip_fh, $sample, $spec->[0], $modbam);
				}
			}
			close $gzip_fh or die "Cannot close $id_partial: $GzipError\n";
			rename $id_partial, $id_file or die "Cannot publish $id_file: $!\n";
			push @checkpoint_outputs, $id_file;
			$published_origin_file = $id_file;
		}
		my $row = {
			mgs => $target->{mgs}, mode => $target->{mode}, sample => $sample,
			scopes => join(',', @scope_names), status => $status,
			candidates => $candidate_count, donor => $donor_count,
			aligned => $aligned_count, mm => $mm_count, ml => $ml_count,
			alignment => $final_alignment,
			alignment_format => $final_alignment eq '-' ? '-' : $output_format,
			index => $final_index,
		};
		push @summary_rows, $row;
		make_path(dirname($unit->{json}));
		write_json_file($unit->{json}, {
			format => 'matafiler-meth2rep-unit-v2',
			mgs => $target->{mgs}, mode => $target->{mode}, sample => $sample,
			created_epoch => 0 + time,
			input_fingerprint => $unit->{parameters}{input_fingerprint},
			representative_mag => $target->{representative_mag},
			reference_fasta => $target->{reference},
			coverage => $coverage,
			processing => { alignment_transfer_wall_seconds => $processing_seconds },
			warnings => \@unit_warnings,
			scope_reports => \@scope_reports,
			tools => {
				samtools => { command => $samtools, version => $tool_versions{samtools}, identity => $tool_identities{samtools} },
				minimap2 => { command => $minimap2, version => $tool_versions{minimap2}, identity => $tool_identities{minimap2} },
				bam_filter => { command => $bam_filter, identity => $tool_identities{bam_filter} },
			},
			output => { alignment => $final_alignment, format => $row->{alignment_format}, index => $final_index },
			policies => {
				supplementary_alignments => $supplementary_alignments,
				allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
			},
			summary => $row, donors => \@unit_donors,
		});
		my @unit_outputs = ($unit->{json});
		push @unit_outputs, $final_alignment, $final_index if $final_alignment ne '-';
		push @unit_outputs, $published_origin_file if $published_origin_file ne '';
		write_checkpoint($unit->{stone},
			parameters => $unit->{parameters}, outputs => \@unit_outputs);
		push @checkpoint_outputs, $unit->{json};
	}
}

my $summary_partial = "$summary_file.part.$$";
push @publication_partials, $summary_partial;
open my $summary, '>', $summary_partial or die "Cannot create $summary_partial: $!\n";
print {$summary} join("\t", qw(mgs mode sample scopes status candidate_reads donor_reads aligned_records MM_records ML_records alignment_format alignment index)), "\n";
for my $row (sort {
	$a->{mgs} cmp $b->{mgs} || $a->{mode} cmp $b->{mode} || $a->{sample} cmp $b->{sample}
} @summary_rows) {
	print {$summary} join("\t", @{$row}{qw(mgs mode sample scopes status candidates donor aligned mm ml alignment_format alignment index)}), "\n";
}
close $summary or die "Cannot close $summary_partial: $!\n";
rename $summary_partial, $summary_file or die "Cannot publish $summary_file: $!\n";

my @manifest_units;
find({
	no_chdir => 1,
	wanted => sub {
		return unless -f $File::Find::name && $File::Find::name =~ /\.json\z/;
		my $json = $File::Find::name;
		(my $stone = $json) =~ s/\.json\z/.stone/;
		return unless -s $stone && checkpoint_valid($stone);
		my $unit = read_json_file($json);
		return unless ($unit->{format} // '') eq 'matafiler-meth2rep-unit-v2'
			&& ref($unit->{donors}) eq 'ARRAY' && ref($unit->{summary}) eq 'HASH';
		push @manifest_units, $unit;
	},
}, File::Spec->catdir($state_dir, 'units'));
my $manifest_partial = "$out_manifest.part.$$";
push @publication_partials, $manifest_partial;
open my $manifest_out, '>', $manifest_partial or die "Cannot create $manifest_partial: $!\n";
print {$manifest_out} join("\t", qw(mgs mode sample scope technology original_modbam donor_candidate_reads total_candidate_reads transferred_records MM_records ML_records status input_validation alignment_format alignment index input_fingerprint completed_epoch)), "\n";
for my $unit (sort {
	$a->{mgs} cmp $b->{mgs} || $a->{mode} cmp $b->{mode} || $a->{sample} cmp $b->{sample}
} @manifest_units) {
	for my $donor (sort {
		$a->{scope} cmp $b->{scope} || $a->{modbam} cmp $b->{modbam}
	} @{$unit->{donors}}) {
		my $row = $unit->{summary};
		my $unit_key = "$unit->{mode}\t$unit->{mgs}\t$unit->{sample}";
		my $input_validation = exists $unit_context{$unit_key}
			? 'current' : 'recorded_only';
		print {$manifest_out} join("\t",
			$unit->{mgs}, $unit->{mode}, $unit->{sample},
			$donor->{scope}, $donor->{technology}, $donor->{modbam},
			$donor->{candidate_reads}, $row->{candidates}, $row->{aligned},
			$row->{mm}, $row->{ml}, $row->{status}, $input_validation,
			$row->{alignment_format}, $row->{alignment}, $row->{index},
			$unit->{input_fingerprint}, $unit->{created_epoch},
		), "\n";
	}
}
close $manifest_out or die "Cannot close $manifest_partial: $!\n";
rename $manifest_partial, $out_manifest or die "Cannot publish $out_manifest: $!\n";

my %units_by_key = map { (join("\t", @{$_}{qw(mode mgs sample)}) => $_) } @manifest_units;
my %mgs_samples;
for my $unit_key (keys %unit_context) {
	my ($mode, $mgs, $sample) = split /\t/, $unit_key, 3;
	$mgs_samples{"$mgs\t$sample"} = 1;
}
for my $mgs_sample (sort keys %mgs_samples) {
	my ($mgs, $sample) = split /\t/, $mgs_sample, 2;
	my ($current_target) = grep { $_->{available} && $_->{mgs} eq $mgs }
		values %{$plan->{targets}};
	next unless $current_target;
	my ($recorded_source, $recorded_paths) = recorded_source_paths($sample, $map->{$sample});
	my @support_paths;
	if (($map->{$sample}{SupportReads} // '') ne '') {
		(undef, my $paths) = parseSupportReads($map->{$sample}{SupportReads});
		@support_paths = @{$paths};
	}
	my %recorded = map {
		my $canonical = -e $_ ? abs_path($_) : File::Spec->rel2abs($_);
		$canonical => 1
	} (@{$recorded_paths}, @support_paths);
	my (%modes_for_log, %declared_donors, @warnings);
	for my $mode (qw(mgs2rep rep2rep)) {
		my $unit_key = "$mode\t$mgs\t$sample";
		my $unit = $units_by_key{$unit_key} or next;
		if (($unit->{representative_mag} // '') ne $current_target->{representative_mag}) {
			push @warnings, "Historical $mode output targets a different representative and is not included in this log";
			next;
		}
		my $row = $unit->{summary};
		my $validation = exists($unit_context{$unit_key}) ? 'current' : 'recorded_only';
		my $alignment = $row->{alignment};
		$modes_for_log{$mode} = {
			status => $row->{status}, input_validation => $validation,
			checkpoint_reused_this_run => exists($cached_units{$unit_key}) ? JSON::PP::true : JSON::PP::false,
			candidate_reads => $row->{candidates}, donor_reads => $row->{donor},
			accepted_alignment_records => $row->{aligned},
			MM_records => $row->{mm}, ML_records => $row->{ml},
			alignment => $alignment, alignment_format => $row->{alignment_format},
			index => $row->{index},
			coverage => $unit->{coverage},
			processing => $unit->{processing},
			tools => $unit->{tools},
			scopes => $unit->{scope_reports},
			donors => $unit->{donors},
			completed_epoch => $unit->{created_epoch},
		};
		$declared_donors{$_->{modbam}} = 1 for @{$unit->{donors}};
		push @warnings, @{$unit->{warnings} || []};
	}
	for my $donor (sort keys %declared_donors) {
		push @warnings, "Declared original modBAM is absent from MATAFILER's recorded input paths (possibly relocated or generated separately): $donor"
			if (%recorded && !$recorded{$donor});
	}
	push @warnings, 'Neither per-sample input_raw.txt nor cohort Input_raw.txt was found; original donor paths could not be independently corroborated'
		if $recorded_source eq '';
	my %unique_warnings;
	@warnings = grep { !$unique_warnings{$_}++ } @warnings;
	my %sample_scopes;
	for my $target (values %{$plan->{targets}}) {
		next unless $target->{available} && $target->{mgs} eq $mgs;
		for my $scope_key (keys %{$target->{scope_keys}}) {
			$sample_scopes{$scope_key} = 1
				if $plan->{scopes}{$scope_key}{sample} eq $sample;
		}
	}
	my $log_file = File::Spec->catfile($out_dir, $mgs, file_component($sample) . '.meth2rep.json');
	make_path(dirname($log_file));
	write_json_file($log_file, {
		format => 'matafiler-meth2rep-sample-v2', component_version => $VERSION,
		mgs => $mgs, sample => $sample,
		representative_mag => $current_target->{representative_mag},
		representative_fasta => $current_target->{reference},
		generated_epoch => 0 + time,
		pipeline_input_record => $recorded_source || undef,
		pipeline_recorded_inputs => $recorded_paths,
		mapping_support_inputs => \@support_paths,
		provenance_rule => 'Manifest donor paths identify source files; each candidate QNAME and native-orientation SEQ digest must also match the assembly CRAM before methylation tags can transfer',
		shared_input_processing_this_run => {
			source_cram_scan_seconds => {
				map { ($_ => $source_scan_seconds{$_}) }
					grep { exists $source_scan_seconds{$_} } sort keys %sample_scopes
			},
			original_modbam_stream_seconds => {
				map { ($_ => $donor_stream_seconds{$_}) }
					grep { exists $donor_stream_seconds{$_} } sort keys %declared_donors
			},
			note => 'Shared scans may serve several MGS/mode units; do not sum these as per-unit exclusive CPU time',
		},
		filters => {
			source_min_mapq => $source_min_mapq,
			source_min_coverage => $source_min_coverage,
			ont_target => \@filter_ont, pb_target => \@filter_pb,
		},
		mapper_presets => { ONT => $preset_ont, PB => $preset_pb },
		output_format => $output_format,
		policies => {
			supplementary_alignments => $supplementary_alignments,
			allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
		},
		modes => \%modes_for_log,
		warnings => \@warnings,
	});
	push @checkpoint_outputs, $log_file;
}

my @input_records = map {
	my @stat = stat($_);
	+{ path => $_, size => 0 + $stat[7], mtime => 0 + $stat[9] };
} @{$plan->{inputs}};
my $provenance = {
	format => 'matafiler-meth2rep-v2', component_version => $VERSION,
	created_epoch => 0 + time, modes => \@modes,
	selected_mgs => [sort keys %{$selected}],
	semantics => {
		candidate_source => 'primary nonduplicate QC-passing CRAM records on source-MAG contigs',
		donor_source => 'original manifest-declared modBAM primary records',
		identity => 'candidate QNAME plus exact native-orientation sequence must agree between the assembly CRAM and declared original modBAM',
		alignment => "fresh original-sequence minimap2 alignment to the representative MAG with -Y and --secondary=no; supplementary output policy is $supplementary_alignments",
		modification_projection => 'first-party exact full-native-sequence MM/ML transfer; MM/ML remain in original read coordinates, MN is regenerated, hard clipping and sequence mismatch fail',
	},
	filters => {
		source_min_mapq => $source_min_mapq,
		ont => join(' ', @filter_ont), pb => join(' ', @filter_pb),
	},
	mapper_presets => { ONT => $preset_ont, PB => $preset_pb },
	output => { format => $output_format, cram_reference => $output_format eq 'cram' ? 'embedded' : 'not_applicable' },
	policies => {
		supplementary_alignments => $supplementary_alignments,
		allow_missing_mn => $allow_missing_mn ? JSON::PP::true : JSON::PP::false,
	},
	resources => { threads => $threads, sort_memory_budget_gb => $memory_gb, samtools_sort_memory_per_thread => $sort_memory },
	tools => {
		samtools => { command => $samtools, version => $tool_versions{samtools}, identity => $tool_identities{samtools} },
		minimap2 => { command => $minimap2, version => $tool_versions{minimap2}, identity => $tool_identities{minimap2} },
		bam_filter => { command => $bam_filter, identity => $tool_identities{bam_filter} },
	},
	inputs => \@input_records,
	input_fingerprint => $checkpoint_parameters{input_fingerprint},
};
my $provenance_partial = "$provenance_file.part.$$";
push @publication_partials, $provenance_partial;
open my $provenance_fh, '>', $provenance_partial or die "Cannot create $provenance_partial: $!\n";
print {$provenance_fh} JSON::PP->new->ascii->canonical->pretty->encode($provenance)
	or die "Cannot write $provenance_partial: $!\n";
close $provenance_fh or die "Cannot close $provenance_partial: $!\n";
rename $provenance_partial, $provenance_file or die "Cannot publish $provenance_file: $!\n";
push @checkpoint_outputs, $summary_file, $provenance_file, $out_manifest;
write_checkpoint(
	$checkpoint_file,
	parameters => \%checkpoint_parameters,
	outputs => \@checkpoint_outputs,
);
print "meth2rep v$VERSION: completed; summary $summary_file\n";
