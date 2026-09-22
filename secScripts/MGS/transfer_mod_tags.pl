#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Both inputs are query-name-sorted BAMs. This intentionally implements only the
# full-original-sequence case; it is not a general clipping/tag-repair engine.
my ($samtools, $donor, $acceptor, $output, $stats) = ('samtools', '', '', '', '');
my $supplementary = 'drop';
my $allow_missing_mn = 0;
GetOptions(
	'samtools=s' => \$samtools,
	'donor=s' => \$donor,
	'acceptor=s' => \$acceptor,
	'output=s' => \$output,
	'stats=s' => \$stats,
	'supplementary=s' => \$supplementary,
	'allow-missing-mn!' => \$allow_missing_mn,
) or die "Usage: transfer_mod_tags.pl --donor NAME.bam --acceptor NAME.bam --output FILE.bam [--samtools PATH] [--stats FILE] [--supplementary drop|keep] [--allow-missing-mn]\n";
die "Donor, acceptor, and output BAM paths are required\n"
	unless $donor ne '' && $acceptor ne '' && $output ne '' && !@ARGV;
die "--supplementary must be 'keep' or 'drop'\n"
	unless $supplementary eq 'keep' || $supplementary eq 'drop';

sub shell_quote {
	my ($value) = @_;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub command {
	my ($tool, @args) = @_;
	return join(' ', $tool, map { shell_quote($_) } @args);
}

sub start_reader {
	my ($description, $command) = @_;
	open my $fh, '-|', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub start_writer {
	my ($description, $command) = @_;
	open my $fh, '|-', 'bash', '-o', 'pipefail', '-c', $command
		or die "Cannot start $description: $!\n";
	return $fh;
}

sub finish {
	my ($fh, $description) = @_;
	return if close $fh;
	die "$description failed (exit " . ($? >> 8) . ")\n";
}

sub fields {
	my ($line, $origin) = @_;
	$line =~ s/[\r\n]+\z//;
	my @fields = split /\t/, $line, -1;
	die "Malformed $origin SAM record\n" unless @fields >= 11;
	return @fields;
}

sub forward_sequence {
	my ($record, $origin) = @_;
	my ($name, $flag, $cigar, $sequence) = @{$record}[0, 1, 5, 9];
	die "$origin record '$name' has no complete SEQ\n"
		if $sequence eq '*' || $sequence eq '';
	die "$origin record '$name' is hard-clipped ($cigar); modification offsets cannot be copied safely\n"
		if $cigar =~ /H/;
	if ($flag & 0x10) {
		$sequence = reverse $sequence;
		$sequence =~ tr/ACGTRYKMSWBDHVNacgtrykmswbdhvn/TGCAYRMKSWVHDBNtgcayrmkswvhdbn/;
	}
	return uc $sequence;
}

sub donor_tags {
	my ($record, $forward_sequence) = @_;
	my ($name, $sequence) = @{$record}[0, 9];
	my ($mm, $ml, $mn);
	for my $field (@{$record}[11 .. $#{$record}]) {
		if ($field =~ /\AM[Mm]:/) {
			die "Donor '$name' has duplicate MM tags\n" if defined $mm;
			$mm = $field;
		} elsif ($field =~ /\AM[Ll]:/) {
			die "Donor '$name' has duplicate ML tags\n" if defined $ml;
			$ml = $field;
		} elsif ($field =~ /\AMN:/) {
			die "Donor '$name' has duplicate MN tags\n" if defined $mn;
			$mn = $field;
		}
	}
	die "Donor '$name' lacks a paired MM/ML tag set\n"
		unless defined($mm) && defined($ml);
	my ($groups) = $mm =~ /\AM[Mm]:Z:(.*)\z/s;
	die "Donor '$name' has malformed MM or ML tag type\n"
		unless defined($groups) && $ml =~ /\AM[Ll]:B:C(?:,\d+)*\z/;
	my @probabilities = split /,/, $ml;
	shift @probabilities; # ML:B:C
	my %base_positions = (map { $_ => [] } qw(A C G T U N));
	for my $position (0 .. length($forward_sequence) - 1) {
		my $base = substr($forward_sequence, $position, 1);
		push @{$base_positions{N}}, $position;
		push @{$base_positions{$base}}, $position if exists $base_positions{$base};
		# SAM stores RNA as T in SEQ in most producers; modkit likewise treats U
		# as the T fundamental base when resolving MM deltas.
		push @{$base_positions{U}}, $position if $base eq 'T';
	}
	my ($expected_probabilities, $probability_index) = (0, 0);
	my %position_probability;
	if ($groups ne '') {
		die "Donor '$name' has an unterminated MM group\n" unless $groups =~ /;\z/;
		for my $group (split /;/, $groups) {
			my ($base, $strand, $codes, $deltas) = $group =~ /\A([ACGTUN])([+-])([A-Za-z]+|\d+)[.?]?(?:,(\d+(?:,\d+)*))?\z/;
			die "Donor '$name' has malformed MM group '$group'\n" unless defined $codes;
			die "Donor '$name' has an MM ChEBI code outside the unsigned 32-bit range\n"
				if $codes =~ /\A\d+\z/ && $codes > 4_294_967_295;
			my @offsets = defined($deltas) ? split(/,/, $deltas) : ();
			my @codes = $codes =~ /\A\d+\z/ ? ($codes) : split(//, $codes);
			my $base_cursor = -1;
			for my $offset (@offsets) {
				die "Donor '$name' has an MM delta outside the unsigned 32-bit range\n"
					if $offset > 4_294_967_295;
				$base_cursor += $offset + 1;
				die "Donor '$name' has an MM call beyond its complete original SEQ\n"
					if $base_cursor > $#{$base_positions{$base}};
				my $position = $base_positions{$base}[$base_cursor];
				for my $code (@codes) {
					die "Donor '$name' has an ML array shorter than its MM calls\n"
						if $probability_index > $#probabilities;
					my $probability = $probabilities[$probability_index++];
					die "Donor '$name' has an ML probability outside 0..255\n"
						if $probability > 255;
					my $key = join("\t", $position, $strand);
					$position_probability{$key} += ($probability + 0.5) / 256;
					die "Donor '$name' has mutually exclusive modification probabilities above 1 at original read position $position\n"
						if $position_probability{$key} > 1.01;
				}
			}
			my $calls = scalar @offsets;
			$expected_probabilities += $calls * scalar(@codes);
		}
	}
	die "Donor '$name' has mismatched MM/ML call counts: expected $expected_probabilities probabilities, found "
		. scalar(@probabilities) . "\n"
		unless @probabilities == $expected_probabilities
		&& $probability_index == @probabilities;
	die "Donor '$name' lacks MN; use --allow-missing-mn only for explicitly reviewed legacy input\n"
		if !defined($mn) && !$allow_missing_mn;
	if (defined $mn) {
		die "Donor '$name' has malformed or stale MN\n"
			unless $mn =~ /\AMN:i:(\d+)\z/ && $1 == length($sequence);
	}
	return ($mm, $ml, defined($mn) ? 0 : 1);
}

sub natural_qname_cmp {
	my ($left, $right) = @_;
	my ($i, $j) = (0, 0);
	my ($left_length, $right_length) = (length($left), length($right));
	while ($i < $left_length && $j < $right_length) {
		my ($a, $b) = (substr($left, $i, 1), substr($right, $j, 1));
		if ($a !~ /[0-9]/ || $b !~ /[0-9]/) {
			return ord($a) - ord($b) if $a ne $b;
			$i++; $j++;
			next;
		}
		$i++ while $i < $left_length && substr($left, $i, 1) eq '0';
		$j++ while $j < $right_length && substr($right, $j, 1) eq '0';
		while ($i < $left_length && $j < $right_length
			&& substr($left, $i, 1) =~ /[0-9]/
			&& substr($right, $j, 1) =~ /[0-9]/
			&& substr($left, $i, 1) eq substr($right, $j, 1)) {
			$i++; $j++;
		}
		my $left_value = $i < $left_length ? ord(substr($left, $i, 1)) : 0;
		my $right_value = $j < $right_length ? ord(substr($right, $j, 1)) : 0;
		my $difference = $left_value - $right_value;
		while ($i < $left_length && $j < $right_length
			&& substr($left, $i, 1) =~ /[0-9]/
			&& substr($right, $j, 1) =~ /[0-9]/) {
			$i++; $j++;
		}
		return 1 if $i < $left_length && substr($left, $i, 1) =~ /[0-9]/;
		return -1 if $j < $right_length && substr($right, $j, 1) =~ /[0-9]/;
		return $difference if $difference;
	}
	return $i < $left_length ? 1 : $j < $right_length ? -1 : 0;
}

sub aligned_native_interval {
	my ($record) = @_;
	my ($name, $flag, $cigar, $sequence) = @{$record}[0, 1, 5, 9];
	my ($query_offset, $start, $end) = (0, undef, undef);
	my $rebuilt = '';
	while ($cigar =~ /([1-9]\d*)([MIDNSHP=X])/g) {
		my ($span, $operation) = ($1, $2);
		$rebuilt .= "$span$operation";
		if ($operation =~ /[MI=X]/) {
			$start = $query_offset unless defined $start;
			$query_offset += $span;
			$end = $query_offset;
		} elsif ($operation eq 'S') {
			$query_offset += $span;
		}
	}
	die "Accepted read '$name' has malformed CIGAR '$cigar'\n"
		unless $cigar ne '*' && $rebuilt eq $cigar && defined($start)
		&& defined($end) && $query_offset == length($sequence);
	if ($flag & 0x10) {
		($start, $end) = (length($sequence) - $end, length($sequence) - $start);
	}
	return ($start, $end);
}

sub take_natural_class {
	my ($first, $fh, $origin) = @_;
	return ([], undef) unless defined $first;
	my $class_name = $first->[0];
	my @records = ($first);
	while (my $line = <$fh>) {
		my @record = fields($line, $origin);
		my $comparison = natural_qname_cmp($record[0], $class_name);
		die ucfirst($origin) . " BAM is not in samtools natural query-name order ('$record[0]' follows '$class_name')\n"
			if $comparison < 0;
		return (\@records, \@record) if $comparison > 0;
		push @records, \@record;
	}
	return (\@records, undef);
}

my $donor_description = "reading name-sorted donor $donor";
my $acceptor_description = "reading name-sorted acceptor $acceptor";
my $output_description = "writing transferred modBAM $output";
my $donor_fh = start_reader($donor_description, command($samtools, 'view', $donor));
my $acceptor_fh = start_reader($acceptor_description, command($samtools, 'view', '-h', $acceptor));
my $output_fh = start_writer($output_description, command($samtools, 'view', '-b', '-o', $output, '-'));
my %counts = (
	input_alignments => 0, output_alignments => 0, output_read_names => 0,
	dropped_no_primary_groups => 0, dropped_no_primary_alignments => 0,
	dropped_supplementary_alignments => 0, stripped_sa_tags => 0,
	legacy_missing_mn_reads => 0,
);

my $transfer_group = sub {
	my ($name, $donor_record, $acceptor_group) = @_;
	my @primary = grep { !(($_->[1] + 0) & (0x100 | 0x800 | 0x4)) } @{$acceptor_group};
	die "Accepted read '$name' contains a secondary alignment despite --secondary=no\n"
		if grep { ($_->[1] + 0) & 0x100 } @{$acceptor_group};
	die "Accepted read '$name' contains an unmapped record inside the mapped acceptor\n"
		if grep { ($_->[1] + 0) & 0x4 } @{$acceptor_group};
	die "Accepted read '$name' has more than one primary representative alignment\n"
		if @primary > 1;
	if (!@primary) {
		$counts{dropped_no_primary_groups}++;
		$counts{dropped_no_primary_alignments} += scalar @{$acceptor_group};
		return;
	}
	my @selected = $supplementary eq 'keep' ? @{$acceptor_group} : @primary;
	$counts{dropped_supplementary_alignments} += @{$acceptor_group} - @selected;
	if ($supplementary eq 'keep' && @selected > 1) {
		my @intervals;
		for my $record (@selected) {
			my ($start, $end) = aligned_native_interval($record);
			for my $prior (@intervals) {
				die "Accepted read '$name' has overlapping primary/supplementary query spans; refusing a double-countable methylation projection\n"
					if $start < $prior->[1] && $prior->[0] < $end;
			}
			push @intervals, [$start, $end];
		}
	}
	my $donor_forward = forward_sequence($donor_record, 'Donor');
	my ($mm, $ml, $missing_mn) = donor_tags($donor_record, $donor_forward);
	$counts{legacy_missing_mn_reads} += $missing_mn;
	for my $acceptor_record (@selected) {
		my $acceptor_forward = forward_sequence($acceptor_record, 'Acceptor');
		die "Accepted read '$name' differs from the complete original modBAM SEQ; cannot transfer MM/ML safely\n"
			unless $acceptor_forward eq $donor_forward;
		my @optional;
		for my $field (@{$acceptor_record}[11 .. $#{$acceptor_record}]) {
			die "Accepted read '$name' unexpectedly already has modification tags\n"
				if $field =~ /\A(?:M[Mm]:|M[Ll]:|MN:)/;
			if ($field =~ /\ASA:Z:/) {
				$counts{stripped_sa_tags}++;
				next;
			}
			push @optional, $field;
		}
		splice @{$acceptor_record}, 11, @{$acceptor_record} - 11, @optional;
		push @{$acceptor_record}, $mm, $ml, 'MN:i:' . length($acceptor_record->[9]);
		print {$output_fh} join("\t", @{$acceptor_record}), "\n"
			or die "Cannot write transferred record '$name': $!\n";
		$counts{output_alignments}++;
	}
	$counts{output_read_names}++;
};

my $donor_first_line = <$donor_fh>;
my $donor_first = defined($donor_first_line) ? [fields($donor_first_line, 'donor')] : undef;
my $acceptor_first;
while (my $line = <$acceptor_fh>) {
	if (!defined($acceptor_first) && $line =~ /\A\@/) {
		print {$output_fh} $line or die "Cannot write acceptor header: $!\n";
		next;
	}
	$acceptor_first = [fields($line, 'acceptor')];
	last;
}
while (defined $acceptor_first) {
	my ($acceptor_class, $next_acceptor) = take_natural_class($acceptor_first, $acceptor_fh, 'acceptor');
	my $class_name = $acceptor_class->[0][0];
	$counts{input_alignments} += scalar @{$acceptor_class};
	while (defined($donor_first) && natural_qname_cmp($donor_first->[0], $class_name) < 0) {
		my (undef, $next_donor) = take_natural_class($donor_first, $donor_fh, 'donor');
		$donor_first = $next_donor;
	}
	die "No original modBAM donor natural-name class for accepted read '$class_name'\n"
		unless defined($donor_first) && natural_qname_cmp($donor_first->[0], $class_name) == 0;
	my ($donor_class, $next_donor) = take_natural_class($donor_first, $donor_fh, 'donor');
	my %donors;
	for my $record (@{$donor_class}) {
		die "Original donor BAM contains more than one record named '$record->[0]'\n"
			if exists $donors{$record->[0]};
		$donors{$record->[0]} = $record;
	}
	my (%acceptors, @name_order);
	for my $record (@{$acceptor_class}) {
		push @name_order, $record->[0] unless exists $acceptors{$record->[0]};
		push @{$acceptors{$record->[0]}}, $record;
	}
	for my $name (@name_order) {
		die "No exact original modBAM donor for accepted read '$name'\n"
			unless exists $donors{$name};
		$transfer_group->($name, $donors{$name}, $acceptors{$name});
	}
	$donor_first = $next_donor;
	$acceptor_first = $next_acceptor;
}
finish($donor_fh, $donor_description);
finish($acceptor_fh, $acceptor_description);
finish($output_fh, $output_description);
if ($stats ne '') {
	open my $stats_fh, '>', $stats or die "Cannot create transfer statistics $stats: $!\n";
	print {$stats_fh} "$_\t$counts{$_}\n" for sort keys %counts;
	close $stats_fh or die "Cannot close transfer statistics $stats: $!\n";
}
print STDERR "Transferred MM/ML for $counts{output_alignments} alignments from $counts{output_read_names} original reads"
	. ($counts{dropped_no_primary_groups}
		? "; dropped $counts{dropped_no_primary_alignments} orphan supplementary alignment(s) from $counts{dropped_no_primary_groups} read(s)"
		: '') . "\n";
