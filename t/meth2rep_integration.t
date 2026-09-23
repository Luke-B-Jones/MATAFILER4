use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3 qw(open3);
use JSON::PP qw(decode_json);
use Symbol qw(gensym);
use Test::More;
use lib File::Spec->catdir($Bin, '..');
use Mods::Checkpoint qw(write_checkpoint);

my $samtools = `command -v samtools 2>/dev/null`;
my $minimap2 = `command -v minimap2 2>/dev/null`;
chomp($samtools, $minimap2);
plan skip_all => 'samtools and minimap2 are required for the integration test'
	unless $samtools ne '' && $minimap2 ne '';

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub run_ok {
	my ($description, @command) = @_;
	my $status = system @command;
	is($status, 0, $description) or BAIL_OUT("command failed: @command");
}

sub capture {
	my (@command) = @_;
	my $stderr = gensym;
	my $pid = open3(undef, my $stdout, $stderr, @command);
	local $/;
	my $out = <$stdout> // '';
	my $err = <$stderr> // '';
	waitpid($pid, 0);
	return ($? >> 8, $out, $err);
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $contents = <$fh> // '';
	close $fh or die "Cannot close $path: $!";
	return $contents;
}

sub reverse_complement {
	my ($sequence) = @_;
	$sequence = reverse $sequence;
	$sequence =~ tr/ACGTNacgtn/TGCANtgcan/;
	return $sequence;
}

my $tmp = tempdir(CLEANUP => 1);
my $run_root = File::Spec->catdir($tmp, 'output', 'fixture');
my $sample_dir = File::Spec->catdir($run_root, 'S1');
my $mapping_dir = File::Spec->catdir($sample_dir, 'mapping');
my $assembly = File::Spec->catdir($tmp, 'assembly');
my $representatives = File::Spec->catdir($tmp, 'representatives');
my $mgs_dir = File::Spec->catdir($tmp, 'mgs', 'Bin_SB');
my $out = File::Spec->catdir($tmp, 'meth2rep');
my $out_manifest = File::Spec->catfile($tmp, 'tracking', 'meth2rep.tsv');
make_path(
	$mapping_dir,
	File::Spec->catdir($sample_dir, 'assemblies', 'metag'),
	File::Spec->catdir($assembly, 'Binning', 'SB'),
	$representatives,
	File::Spec->catdir($mgs_dir, 'LOGandSUB', 'checkpoints'),
	File::Spec->catdir($tmp, 'raw'),
);

my $seed = 17;
my @bases = qw(A C G T);
my $reference_sequence = '';
for (1 .. 1200) {
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$reference_sequence .= $bases[($seed >> 16) % 4];
}
my $read_sequence = substr($reference_sequence, 200, 600);
$read_sequence =~ s/^[^C]*//; # ensure MM delta zero addresses a canonical C
my $read_length = length($read_sequence);
my $quality = 'I' x $read_length;
my $reverse_donor_sequence = reverse_complement($read_sequence);
my $member_read_sequence = substr($reference_sequence, 400, 600);
$member_read_sequence =~ s/^[^C]*//;
my $member_read_length = length($member_read_sequence);
my $member_quality = 'I' x $member_read_length;
my $short_source_cigar = '50S10M' . ($read_length - 60) . 'S';
my $assembly_fasta = File::Spec->catfile($assembly, 'scaffolds.fasta.filt');
write_file($assembly_fasta,
	">S1-contig\n$reference_sequence\n>S1-member-contig\n$reference_sequence\n");
write_file(File::Spec->catfile($sample_dir, 'assemblies', 'metag', 'assembly.txt'), "$assembly\n");
write_file(File::Spec->catfile($assembly, 'Binning', 'SB', 'S1'),
	"Sequence ID\tBin\nS1-contig\t1.fa.gz\nS1-member-contig\t2.fa.gz\n");

my $representative = File::Spec->catfile($representatives, 'MGS.1.ctgs.S1__1.fa.gz');
gzip(\(">S1-contig\n$reference_sequence\n") => $representative)
	or die "Cannot write $representative: $GzipError";
my $report = File::Spec->catfile($tmp, 'MAGvsGC.txt.gz');
gzip(\("MAG\tMGS\tRepresentative4MGS\n"
		. "S1__1.fa.gz\tMGS.1\t*\nS1__2.fa.gz\tMGS.1\t\n"
		. "Cano__7\tMGS.2\t*\n") => $report)
	or die "Cannot write $report: $GzipError";

my $donor_sam = File::Spec->catfile($tmp, 'donor.sam');
my $donor_bam = File::Spec->catfile($tmp, 'donor.bam');
my $donor2_sam = File::Spec->catfile($tmp, 'donor2.sam');
my $donor2_bam = File::Spec->catfile($tmp, 'donor2.bam');
write_file($donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:donor-ref\tLN:1200\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$read_sequence\t$quality\tMM:Z:C+m.,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r3\t16\tdonor-ref\t201\t60\t${read_length}M\t*\t0\t0\t$reverse_donor_sequence\t$quality\tMM:Z:C+m.,0;\tML:B:C,210\tMN:i:$read_length\tNM:i:0\n");
run_ok('synthetic donor modBAM is created',
	$samtools, 'view', '-b', '-o', $donor_bam, $donor_sam);
write_file($donor2_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n"
	. "r2\t4\t*\t0\t0\t*\t*\t0\t0\t$member_read_sequence\t$member_quality\tMM:Z:C+m.,0;\tML:B:C,180\tMN:i:$member_read_length\n");
run_ok('second original modBAM is created',
	$samtools, 'view', '-b', '-o', $donor2_bam, $donor2_sam);

my $candidate_sam = File::Spec->catfile($tmp, 'candidate.sam');
my $candidate_cram = File::Spec->catfile($mapping_dir, 'S1-smd.cram');
write_file($candidate_sam,
	"\@HD\tVN:1.6\tSO:coordinate\n\@SQ\tSN:S1-contig\tLN:1200\n\@SQ\tSN:S1-member-contig\tLN:1200\n"
	. "r1\t0\tS1-contig\t201\t60\t${read_length}M\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n"
	. "r2\t0\tS1-member-contig\t401\t60\t${member_read_length}M\t*\t0\t0\t$member_read_sequence\t$member_quality\tNM:i:0\n"
	. "r3\t0\tS1-contig\t201\t60\t${read_length}M\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n"
	. "r4\t0\tS1-contig\t201\t60\t$short_source_cigar\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n");
run_ok('synthetic assembly backmapping CRAM is created',
	$samtools, 'view', '-C', '-T', $assembly_fasta, '-o', $candidate_cram, $candidate_sam);
write_file("$candidate_cram.sto", "done\n");
my @assembly_stat = stat($assembly_fasta);
write_file(File::Spec->catfile($mapping_dir, 'S1-smd.reference.stat'),
	"$assembly_stat[7] $assembly_stat[9]\n");

my $map = File::Spec->catfile($tmp, 'fixture.map');
write_file($map,
	"#SmplID\tPath\tAssmblGrps\tSupportReads\tINFO\tSeqTech\n"
	. "#RunID\tfixture\n#OutPath\t" . File::Spec->catdir($tmp, 'output') . "/\n"
	. "#DirPath\t" . File::Spec->catdir($tmp, 'raw') . "/\n"
	. "#WARNING\tOFF\nS1\tinput\tG1\t\t\tONT\n");
my $manifest = File::Spec->catfile($tmp, 'modbams.tsv');
write_file($manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$donor_bam\n"
	. "S1\tprimary\tONT\t$donor2_bam\n");
write_file(File::Spec->catfile($sample_dir, 'input_raw.txt'), "$donor_bam;$donor2_bam");
write_checkpoint(File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone'),
	parameters => { stage => 'stage-1' }, outputs => [$report]);
my $checkpoint_representative = $representative;
$checkpoint_representative =~ s{/representatives/}{/representatives//};
write_checkpoint(File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone'),
	parameters => { stage => 'extract-bin-contigs' }, outputs => [$report, $checkpoint_representative]);

my $script = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'meth2rep.pl');
my $bam_filter = 'perl ' . File::Spec->catfile($Bin, '..', 'secScripts', 'assemblies', 'bamFilter.pl');
sub output_bam {
	my ($out_dir, $mode) = @_;
	return File::Spec->catfile($out_dir, 'MGS.1',
		join('__', 'S1', $mode, 'S1__1.fa.gz') . '.mod.bam');
}
sub output_alignment {
	my ($out_dir, $mode, $format) = @_;
	return File::Spec->catfile($out_dir, 'MGS.1',
		join('__', 'S1', $mode, 'S1__1.fa.gz') . ".mod.$format");
}

sub alignment_signatures {
	my ($sam) = @_;
	my @signatures;
	for my $line (grep { $_ ne '' && $_ !~ /^\@/ } split /\n/, $sam) {
		my @fields = split /\t/, $line, -1;
		my %tags;
		for my $field (@fields[11 .. $#fields]) {
			$tags{MM} = $field if $field =~ /^M[Mm]:/;
			$tags{ML} = $field if $field =~ /^M[Ll]:/;
			$tags{MN} = $field if $field =~ /^MN:/;
		}
		push @signatures, join("\t", @fields[0 .. 5], $fields[9], @tags{qw(MM ML MN)});
	}
	return [sort @signatures];
}
my @command = (
	$^X, '-I' . File::Spec->catdir($Bin, '..'), $script,
	'--mgs-dir', $mgs_dir,
	'--map', $map, '--mgs-report', $report,
	'--representatives-dir', $representatives, '--binner', 'SB',
	'--modbam-manifest', $manifest, '--out', $out, '--mgs', 'MGS.1',
	'--out-manifest', $out_manifest,
	'--mgs2rep', '--rep2rep', '--threads', 2,
	'--keep-read-ids',
	'--samtools', $samtools, '--minimap2', $minimap2,
	'--bam-filter', $bam_filter,
);
my ($status, $stdout, $stderr) = capture(@command);
is($status, 0, 'meth2rep completes through real CRAM/minimap2 and internal tag transfer')
	or diag($stdout, $stderr);
ok(-s $out_manifest, 'a custom completed-unit donor manifest is published');
my $tracking = slurp($out_manifest);
like($tracking, qr/^MGS\.1\tmgs2rep\tS1\tprimary\tONT\t\Q$donor_bam\E\t2\t3\t3\t3\t3\tcomplete\t/m,
	'the first donor has two attributed MGS candidate reads');
like($tracking, qr/^MGS\.1\tmgs2rep\tS1\tprimary\tONT\t\Q$donor2_bam\E\t1\t3\t3\t3\t3\tcomplete\t/m,
	'the second donor has one attributed MGS candidate read');
like($tracking, qr/^MGS\.1\trep2rep\tS1\tprimary\tONT\t\Q$donor2_bam\E\t0\t2\t2\t2\t2\tcomplete\t/m,
	'the donor manifest retains zero-use donor provenance per requested mode');
my $sample_log_path = File::Spec->catfile($out, 'MGS.1', 'S1.meth2rep.json');
ok(-s $sample_log_path, 'one sample-level interrogation log is published inside the MGS directory');
my $sample_log = decode_json(slurp($sample_log_path));
is_deeply([sort keys %{$sample_log->{modes}}], [qw(mgs2rep rep2rep)],
	'the single sample log reports both requested modes');
is($sample_log->{modes}{mgs2rep}{coverage}{reference_bases}, 1200,
	'the sample log records the representative reference denominator');
ok($sample_log->{modes}{mgs2rep}{coverage}{covered_bases} > 0
	&& $sample_log->{modes}{mgs2rep}{coverage}{breadth_fraction} > 0
	&& $sample_log->{modes}{mgs2rep}{coverage}{mean_depth} > 0,
	'the sample log records covered bases, breadth, and mean depth');
is($sample_log->{modes}{mgs2rep}{scopes}[0]{filter_stats}{malformed}, 0,
	'the per-scope alignment-filter diagnostics are retained');
like($sample_log->{modes}{mgs2rep}{scopes}[0]{historical_source_filter_provenance},
	qr/^inherited_unknown:/,
	'the log does not pretend Meth2Rep source thresholds recover the historical CRAM filter');
is_deeply($sample_log->{pipeline_recorded_inputs}, [$donor_bam, $donor2_bam],
	'pipeline input_raw provenance is cross-checked and exposed beside the donor manifest');
ok(!(grep { /absent from MATAFILER/ } @{$sample_log->{warnings}}),
	'matching recorded inputs produce no donor-provenance warning');
opendir my $mgs_dir_handle, File::Spec->catdir($out, 'MGS.1') or die $!;
my @mgs_files = sort grep { $_ ne '.' && $_ ne '..' } readdir $mgs_dir_handle;
closedir $mgs_dir_handle;
is(scalar(@mgs_files), 5, 'the MGS folder contains two mode BAMs, their indexes, and one sample log');
my $origins_gzip = File::Spec->catfile($out, '.meth2rep', 'read_ids', 'mgs2rep', 'MGS.1', 'S1.read_origins.tsv.gz');
my $origins = '';
gunzip($origins_gzip => \$origins) or die "Cannot read $origins_gzip: $GunzipError";
like($origins, qr/^S1\tprimary\tr1\t\Q$donor_bam\E$/m,
	'compressed per-read origins identify the first donor');
like($origins, qr/^S1\tprimary\tr2\t\Q$donor2_bam\E$/m,
	'compressed per-read origins identify the second donor');
my $initial_inode = (stat(output_bam($out, 'mgs2rep')))[1];

for my $mode (qw(mgs2rep rep2rep)) {
	my $bam = output_bam($out, $mode);
	if (!-s $bam) {
		my $summary_file = File::Spec->catfile($out, '.meth2rep', 'summary.tsv');
		my $summary_text = '';
		if (-s $summary_file) {
			open my $summary_fh, '<', $summary_file or die "Cannot read $summary_file: $!";
			{ local $/; $summary_text = <$summary_fh> // ''; }
			close $summary_fh;
		}
		diag("meth2rep stdout:\n$stdout\nstderr:\n$stderr\nsummary:\n$summary_text");
	}
	ok(-s $bam && -s "$bam.bai", "$mode publishes an indexed compact modBAM");
	my ($view_status, $view, $view_error) = capture($samtools, 'view', $bam);
	is($view_status, 0, "$mode modBAM is readable") or diag($view_error);
	my @records = grep { $_ ne '' } split /\n/, $view;
	is(scalar(@records), $mode eq 'mgs2rep' ? 3 : 2,
		"$mode contains exactly its intended all-member or representative-only read set");
	like($view, qr/^r1\t/m, "$mode retains the representative-MAG read");
	like($view, qr/^r3\t/m,
		"$mode recovers the original orientation from a reverse-aligned donor");
	if ($mode eq 'mgs2rep') {
		like($view, qr/^r2\t/m, 'mgs2rep includes the non-representative member-MAG read');
	} else {
		unlike($view, qr/^r2\t/m, 'rep2rep excludes the non-representative member-MAG read');
	}
	like($view, qr/\tMM:Z:C\+m\.,0;/, "$mode retains the MM methylation tag");
	like($view, qr/\tML:B:C,220(?:\t|\n)/, "$mode retains the ML probability tag");
	like($view, qr/^r3\t.*\tML:B:C,210(?:\t|$)/m,
		"$mode preserves the reverse-aligned donor's modification probability");
}

my $cram_out = File::Spec->catdir($tmp, 'cram-output');
my @cram_command = grep { $_ ne '--rep2rep' } @command;
for my $i (0 .. $#cram_command - 1) {
	$cram_command[$i + 1] = $cram_out if $cram_command[$i] eq '--out';
	$cram_command[$i + 1] = File::Spec->catfile($cram_out, 'manifest.tsv')
		if $cram_command[$i] eq '--out-manifest';
}
push @cram_command, '--output-format', 'cram';
my ($cram_run_status, $cram_run_stdout, $cram_run_stderr) = capture(@cram_command);
is($cram_run_status, 0, 'self-contained CRAM output completes from an ordinary gzip representative')
	or diag($cram_run_stdout, $cram_run_stderr);
my $cram = output_alignment($cram_out, 'mgs2rep', 'cram');
ok(-s $cram && -s "$cram.crai", 'CRAM output and CRAI use explicit predictable suffixes');
ok(!-e output_alignment($cram_out, 'mgs2rep', 'bam'),
	'CRAM selection does not leave a duplicate BAM output');
ok(!-e "$representative.fai" && !-e "$representative.gzi",
	'CRAM encoding does not write indexes beside the immutable gzip representative');
my ($cram_view_status, $cram_view, $cram_view_error) = capture($samtools, 'view', $cram);
is($cram_view_status, 0, 'embedded-reference CRAM decodes without an external -T reference')
	or diag($cram_view_error);
my (undef, $bam_view_for_cram, undef) = capture($samtools, 'view', output_bam($out, 'mgs2rep'));
is_deeply(alignment_signatures($cram_view), alignment_signatures($bam_view_for_cram),
	'BAM and CRAM preserve identical alignment fields and MM/ML/MN payloads');
my ($idxstats_status, undef, $idxstats_error) = capture($samtools, 'idxstats', $cram);
is($idxstats_status, 0, 'published CRAI is usable by samtools idxstats') or diag($idxstats_error);
my $cram_log = decode_json(slurp(File::Spec->catfile($cram_out, 'MGS.1', 'S1.meth2rep.json')));
is($cram_log->{modes}{mgs2rep}{alignment_format}, 'cram',
	'sample log identifies the actual compact output format');
is($cram_log->{modes}{mgs2rep}{index}, "$cram.crai",
	'sample log records the exact CRAI path');
my $cram_inode = (stat($cram))[1];
my ($cram_resume_status, undef, $cram_resume_error) = capture(@cram_command);
is($cram_resume_status, 0, 'an identical CRAM invocation resumes') or diag($cram_resume_error);
is((stat($cram))[1], $cram_inode, 'CRAM resume does not rewrite a validated alignment');
my ($format_switch_status, undef, $format_switch_error) = capture(
	@cram_command, '--output-format', 'bam');
ok($format_switch_status != 0, 'changing output format cannot silently leave BAM and CRAM siblings');
like($format_switch_error, qr/alternate-format meth2rep output already exists/,
	'format-switch refusal explains the required override');
my ($format_override_status, $format_override_stdout, $format_override_error) = capture(
	@cram_command, '--output-format', 'bam', '--override');
is($format_override_status, 0, 'explicit override safely replaces CRAM with BAM')
	or diag($format_override_stdout, $format_override_error);
ok(-s output_alignment($cram_out, 'mgs2rep', 'bam')
	&& !-e $cram && !-e "$cram.crai",
	'format replacement removes the prior tracked CRAM and CRAI only after BAM publication');

my ($resume_status, $resume_stdout, $resume_stderr) = capture(@command);
is($resume_status, 0, 'a second identical invocation resumes successfully')
	or diag($resume_stderr);
like($resume_stdout, qr/completed/,
	'identical input/options retain completed per-unit results');
is((stat(output_bam($out, 'mgs2rep')))[1], $initial_inode,
	'a valid per-unit resume does not remap or replace an existing modBAM');

my $checkpoint = File::Spec->catfile($out, '.meth2rep', 'complete.stone');
my $checkpoint_before_preview = slurp($checkpoint);
my ($preview_status, $preview_stdout, $preview_stderr) = capture(@command, '--plan-only');
is($preview_status, 0, 'plan-only validates without executing the mapper')
	or diag($preview_stdout, $preview_stderr);
my $preview = File::Spec->catfile($out, '.meth2rep', 'plan.preview.tsv');
ok(-s $preview, 'plan-only writes a separate preview plan');
is(slurp($checkpoint), $checkpoint_before_preview,
	'plan-only does not invalidate or rewrite a completed checkpoint');
my @preview_lines = grep { $_ ne '' } split /\n/, slurp($preview);
ok(!(grep { scalar(split /\t/, $_, -1) != 12 } @preview_lines),
	'preview plan remains a regular twelve-column TSV when sample scopes are embedded');
my ($invalid_preset_status, undef, $invalid_preset_error) = capture(
	@command, '--minimap2-preset-ont', 'definitely-invalid-preset', '--plan-only');
ok($invalid_preset_status != 0,
	'plan-only rejects a minimap2 preset unsupported by the installed executable');
like($invalid_preset_error, qr/(?:unknown preset|validating minimap2 preset).*definitely-invalid-preset/is,
	'unsupported-preset preflight identifies the requested preset');
my @path_discovery_command;
for (my $i = 0; $i <= $#command; $i++) {
	if ($command[$i] eq '--samtools' || $command[$i] eq '--minimap2'
		|| $command[$i] eq '--bam-filter') {
		$i++;
		next;
	}
	push @path_discovery_command, $command[$i];
}
for my $i (0 .. $#path_discovery_command - 1) {
	$path_discovery_command[$i + 1] = File::Spec->catdir($tmp, 'path-preflight')
		if $path_discovery_command[$i] eq '--out';
	$path_discovery_command[$i + 1] = File::Spec->catfile($tmp, 'path-preflight', 'manifest.tsv')
		if $path_discovery_command[$i] eq '--out-manifest';
}
my ($path_status, undef, $path_error);
{
	local $ENV{MGTKDIR};
	delete $ENV{MGTKDIR};
	($path_status, undef, $path_error) = capture(@path_discovery_command, '--plan-only');
}
is($path_status, 0,
	'standalone preflight resolves samtools/minimap2 from PATH and the co-located bamFilter without MGTKDIR')
	or diag($path_error);
my ($manifest_collision_status, undef, $manifest_collision_error) = capture(
	@command, '--out-manifest', $sample_log_path, '--plan-only');
ok($manifest_collision_status != 0,
	'custom manifest cannot overwrite a per-MGS alignment, index, or sample log');
like($manifest_collision_error, qr/cannot be placed inside a per-MGS output directory/,
	'manifest collision fails during preflight with the unsafe path');
my $historical_dir = File::Spec->catdir($out, 'MGS.historical');
make_path($historical_dir);
my $historical_alignment = File::Spec->catfile($historical_dir, 'prior.mod.bam');
write_file($historical_alignment, "historical-alignment-data\n");
my $historical_before = slurp($historical_alignment);
my $historical_inode = (stat($historical_alignment))[1];
my ($historical_collision_status, undef, $historical_collision_error) = capture(
	@command, '--out-manifest', $historical_alignment, '--plan-only');
ok($historical_collision_status != 0,
	'a custom manifest cannot target an unselected historical MGS data directory');
like($historical_collision_error, qr/cannot be placed inside a per-MGS output directory/,
	'historical-data collision is rejected during preflight');
is(slurp($historical_alignment), $historical_before,
	'rejected manifest placement leaves historical MGS data unchanged');
is((stat($historical_alignment))[1], $historical_inode,
	'rejected manifest placement does not replace the historical file');
my @short_out_command = @command;
for my $i (0 .. $#short_out_command - 1) {
	$short_out_command[$i] = '-o' if $short_out_command[$i] eq '--out';
	$short_out_command[$i + 1] = File::Spec->catdir($tmp, 'short-out')
		if $short_out_command[$i] eq '-o';
	$short_out_command[$i + 1] = File::Spec->catfile($tmp, 'short-out', 'manifest.tsv')
		if $short_out_command[$i] eq '--out-manifest';
}
my ($short_out_status, undef, $short_out_error) = capture(@short_out_command, '--plan-only');
is($short_out_status, 0, 'the single short -o output option is accepted')
	or diag($short_out_error);

my @mgs_only_command = grep { $_ ne '--rep2rep' } @command;
my ($mgs_only_status, $mgs_only_stdout, $mgs_only_stderr) = capture(@mgs_only_command);
is($mgs_only_status, 0, 'a narrower mode rerun completes')
	or diag($mgs_only_stdout, $mgs_only_stderr);
ok(-s output_bam($out, 'mgs2rep'),
	'the still-requested mgs2rep output remains published');
ok(-s output_bam($out, 'rep2rep'),
	'a narrower invocation preserves previously completed unselected modes');
like(slurp($out_manifest), qr/^MGS\.1\trep2rep\tS1\t/m,
	'the tracking manifest retains previously completed unselected modes');

my $unready_dir = File::Spec->catdir($tmp, 'unready', 'Bin_SB');
make_path($unready_dir);
my @unready_command = @command;
for my $i (0 .. $#unready_command - 1) {
	$unready_command[$i + 1] = $unready_dir if $unready_command[$i] eq '--mgs-dir';
}
my ($unready_status, undef, $unready_error) = capture(@unready_command, '--plan-only');
ok($unready_status != 0, 'missing MGS progress checkpoints block even a plan preview');
like($unready_error, qr/progress checkpoint is missing/, 'readiness failure names the missing checkpoint');
my $stale_report = File::Spec->catfile($tmp, 'post-checkpoint-MAGvsGC.txt.gz');
write_file($stale_report, slurp($report));
utime(time + 30, time + 30, $stale_report) or die "Cannot age $stale_report: $!";
my @stale_report_command = @command;
for my $i (0 .. $#stale_report_command - 1) {
	$stale_report_command[$i + 1] = $stale_report if $stale_report_command[$i] eq '--mgs-report';
}
my ($stale_report_status, undef, $stale_report_error) = capture(@stale_report_command, '--plan-only');
ok($stale_report_status != 0, 'a membership report newer than Stage1 completion is rejected');
like($stale_report_error, qr/membership report was modified after/, 'the stale-report error explains the provenance mismatch');

my $collision_sam = File::Spec->catfile($tmp, 'collision.sam');
my $collision_bam = File::Spec->catfile($tmp, 'collision.bam');
write_file($collision_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$read_sequence\t$quality\tMM:Z:C+m.,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r2\t4\t*\t0\t0\t*\t*\t0\t0\t$member_read_sequence\t$member_quality\tMM:Z:C+m.,0;\tML:B:C,180\tMN:i:$member_read_length\n");
run_ok('colliding donor modBAM is created', $samtools, 'view', '-b', '-o', $collision_bam, $collision_sam);
my $collision_manifest = File::Spec->catfile($tmp, 'collision-modbams.tsv');
write_file($collision_manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$donor_bam\n"
	. "S1\tprimary\tONT\t$collision_bam\n");
my @collision_command = @command;
for my $i (0 .. $#collision_command - 1) {
	$collision_command[$i + 1] = $collision_manifest if $collision_command[$i] eq '--modbam-manifest';
	$collision_command[$i + 1] = File::Spec->catdir($tmp, 'collision-out') if $collision_command[$i] eq '--out';
}
my ($collision_status, undef, $collision_error) = capture(@collision_command);
ok($collision_status != 0, 'a QNAME shared across original modBAMs fails closed');
like($collision_error, qr/more than one primary record named 'r1'/,
	'the collision error identifies the ambiguous read');
my $wrong_sequence = $read_sequence;
substr($wrong_sequence, -1, 1) = substr($wrong_sequence, -1, 1) eq 'A' ? 'C' : 'A';
my $wrong_donor_sam = File::Spec->catfile($tmp, 'wrong-origin.sam');
my $wrong_donor_bam = File::Spec->catfile($tmp, 'wrong-origin.bam');
write_file($wrong_donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:donor-ref\tLN:1200\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$wrong_sequence\t$quality\tMM:Z:C+m.,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r3\t16\tdonor-ref\t201\t60\t${read_length}M\t*\t0\t0\t$reverse_donor_sequence\t$quality\tMM:Z:C+m.,0;\tML:B:C,210\tMN:i:$read_length\tNM:i:0\n");
run_ok('same-QNAME wrong-sequence donor is created', $samtools, 'view', '-b',
	'-o', $wrong_donor_bam, $wrong_donor_sam);
my $wrong_manifest = File::Spec->catfile($tmp, 'wrong-origin.tsv');
write_file($wrong_manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$wrong_donor_bam\n");
my @wrong_origin_command = grep { $_ ne '--mgs2rep' } @command;
for my $i (0 .. $#wrong_origin_command - 1) {
	$wrong_origin_command[$i + 1] = $wrong_manifest
		if $wrong_origin_command[$i] eq '--modbam-manifest';
	$wrong_origin_command[$i + 1] = File::Spec->catdir($tmp, 'wrong-origin-out')
		if $wrong_origin_command[$i] eq '--out';
	$wrong_origin_command[$i + 1] = File::Spec->catfile($tmp, 'wrong-origin-out', 'manifest.tsv')
		if $wrong_origin_command[$i] eq '--out-manifest';
}
my ($wrong_origin_status, undef, $wrong_origin_error) = capture(@wrong_origin_command);
ok($wrong_origin_status != 0,
	'a coincidentally matching QNAME in the wrong donor cannot acquire methylation tags');
like($wrong_origin_error, qr/native SEQ.*differs between the assembly CRAM and declared original modBAM/,
	'wrong-donor failure requires sequence provenance in addition to QNAME');
my $untracked_out = File::Spec->catdir($tmp, 'untracked-out');
make_path(File::Spec->catdir($untracked_out, 'MGS.1'));
my $untracked_bam = output_bam($untracked_out, 'mgs2rep');
write_file($untracked_bam, 'user-owned fixture');
my @untracked_command = @command;
for my $i (0 .. $#untracked_command - 1) {
	$untracked_command[$i + 1] = $untracked_out if $untracked_command[$i] eq '--out';
	$untracked_command[$i + 1] = File::Spec->catfile($untracked_out, 'manifest.tsv')
		if $untracked_command[$i] eq '--out-manifest';
}
my ($untracked_status, undef, $untracked_error) = capture(@untracked_command);
ok($untracked_status != 0, 'an existing uncheckpointed BAM is not overwritten');
like($untracked_error, qr/Untracked output alignment or index already exists/, 'the overwrite refusal identifies the exact condition');
is(slurp($untracked_bam), 'user-owned fixture', 'the pre-existing file remains unchanged');
my ($invalid_coverage_status, undef, $invalid_coverage_error) = capture(
	@command, '--source-min-coverage', 1.2, '--plan-only');
ok($invalid_coverage_status != 0, 'an invalid source-coverage threshold is rejected');
like($invalid_coverage_error, qr/source-min-coverage must be between 0 and 1/,
	'coverage validation names the offending option');
my @no_source_floor_command = @command;
for my $i (0 .. $#no_source_floor_command - 1) {
	$no_source_floor_command[$i + 1] = File::Spec->catdir($tmp, 'no-source-floor')
		if $no_source_floor_command[$i] eq '--out';
	$no_source_floor_command[$i + 1] = File::Spec->catfile($tmp, 'no-source-floor', 'manifest.tsv')
		if $no_source_floor_command[$i] eq '--out-manifest';
}
my ($no_floor_status, undef, $no_floor_error) = capture(
	@no_source_floor_command, '--source-min-coverage', 0);
ok($no_floor_status != 0, 'disabling source coverage exposes a low-overlap read absent from the donors');
like($no_floor_error, qr/candidate read name\(s\) are absent.*r4/s,
	'the source-coverage failure identifies the otherwise excluded read');
my @high_mapq_command = @command;
for my $i (0 .. $#high_mapq_command - 1) {
	$high_mapq_command[$i + 1] = File::Spec->catdir($tmp, 'high-mapq')
		if $high_mapq_command[$i] eq '--out';
	$high_mapq_command[$i + 1] = File::Spec->catfile($tmp, 'high-mapq', 'manifest.tsv')
		if $high_mapq_command[$i] eq '--out-manifest';
}
my ($high_mapq_status, $high_mapq_stdout, $high_mapq_error) = capture(
	@high_mapq_command, '--target-min-mapq', 200, '--target-min-coverage', 0.7,
	'--target-max-edit-rate', 0.2, '--target-min-end-clip', 5);
is($high_mapq_status, 0, 'named target MAPQ, coverage, edit-rate, and clipping controls execute end to end')
	or diag($high_mapq_stdout, $high_mapq_error);
like(slurp(File::Spec->catfile($tmp, 'high-mapq', '.meth2rep', 'summary.tsv')),
	qr/^MGS\.1\tmgs2rep\tS1\tprimary\tno_target_alignment\t3\t3\t0\t/m,
	'a high MAPQ floor rejects every representative alignment while preserving donor accounting');
my $high_mapq_log = decode_json(slurp(File::Spec->catfile(
	$tmp, 'high-mapq', 'MGS.1', 'S1.meth2rep.json')));
is_deeply($high_mapq_log->{filters}{ont_target}, [0.2, 0.7, 200, 5],
	'sample log records the effective named target-filter overrides in bamFilter order');
my ($override_status, $override_stdout, $override_error) = capture(@command, '--override');
is($override_status, 0, '--override rebuilds previously complete selected units')
	or diag($override_stdout, $override_error);
isnt((stat(output_bam($out, 'mgs2rep')))[1], $initial_inode,
	'override replaces the selected output rather than reusing its old inode');

my $canopy_out = File::Spec->catdir($tmp, 'canopy-only');
my @canopy_command = (
	$^X, '-I' . File::Spec->catdir($Bin, '..'), $script,
	'--mgs-dir', $mgs_dir,
	'--map', $map, '--mgs-report', $report,
	'--representatives-dir', $representatives, '--binner', 'SB',
	'--modbam-manifest', $manifest, '--out', $canopy_out,
	'--mgs', 'MGS.2', '--rep2rep',
	'--samtools', '/not/required/samtools', '--minimap2', '/not/required/minimap2',
	'--bam-filter', '/not/required/bamFilter',
);
my ($canopy_status, $canopy_stdout, $canopy_stderr) = capture(@canopy_command);
is($canopy_status, 0, 'a Canopy-only selection reports biological unavailability without mapping tools')
	or diag($canopy_stdout, $canopy_stderr);
like(slurp(File::Spec->catfile($canopy_out, '.meth2rep', 'summary.tsv')),
	qr/^MGS\.2\trep2rep\t-\t-\tcanopy_only_no_mag_reference\t/m,
	'Canopy-only output carries the explicit unavailable state');
ok(-s File::Spec->catfile($canopy_out, '.meth2rep', 'complete.stone'),
	'Canopy-only unavailability is a checkpointed successful result');

# Exercise our deliberately narrower transfer contract independently of minimap2.
# Soft clipping preserves the complete SEQ, while hard clipping must fail closed.
my $clip_sequence = 'ACCCGTTACGCCATGCTAACCGT';
my $clip_quality = 'I' x length($clip_sequence);
my $trimmed_sequence = substr($clip_sequence, 7, 12);
my $trimmed_quality = 'I' x length($trimmed_sequence);
my $clip_donor_sam = File::Spec->catfile($tmp, 'clip-donor.sam');
my $clip_donor_bam = File::Spec->catfile($tmp, 'clip-donor.bam');
write_file($clip_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m.,0,2,2;\tML:B:C,201,202,250\tMN:i:" . length($clip_sequence) . "\n"
	. "indel\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m.,0;\tML:B:C,155\tMN:i:" . length($clip_sequence) . "\n"
	. "reverseclip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m.,0;\tML:B:C,177\tMN:i:" . length($clip_sequence) . "\n");
run_ok('clipping-contract donor BAM is created',
	$samtools, 'view', '-b', '-o', $clip_donor_bam, $clip_donor_sam);
my $clip_acceptor_sam = File::Spec->catfile($tmp, 'clip-acceptor.sam');
my $clip_acceptor_bam = File::Spec->catfile($tmp, 'clip-acceptor.bam');
write_file($clip_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
		. "clip\t0\tref\t1\t60\t12M11S\t*\t0\t0\t$clip_sequence\t$clip_quality\tNM:i:0\tSA:Z:ref,30,+,12S11M,60,0;\n"
		. "clip\t2048\tref\t30\t60\t12S11M\t*\t0\t0\t$clip_sequence\t$clip_quality\tNM:i:0\tSA:Z:ref,1,+,12M11S,60,0;\n"
		. "indel\t0\tref\t40\t60\t6M1I6M1D10M\t*\t0\t0\t$clip_sequence\t$clip_quality\tNM:i:2\n"
		. "reverseclip\t16\tref\t60\t60\t23M\t*\t0\t0\t" . reverse_complement($clip_sequence)
	. "\t$clip_quality\tNM:i:0\n");
run_ok('clipping-contract acceptor BAM is created',
	$samtools, 'view', '-b', '-o', $clip_acceptor_bam, $clip_acceptor_sam);
my $transfer_script = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'transfer_mod_tags.pl');
my $clip_transferred_bam = File::Spec->catfile($tmp, 'clip-transferred.bam');
my $clip_transfer_stats = File::Spec->catfile($tmp, 'clip-transfer.tsv');
run_ok('internal transfer accepts full-sequence soft-clipped, supplementary, and reverse records',
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $clip_acceptor_bam, '--output', $clip_transferred_bam,
	'--supplementary', 'keep', '--stats', $clip_transfer_stats);
my ($clip_status, $clip_view, $clip_error) = capture($samtools, 'view', $clip_transferred_bam);
is($clip_status, 0, 'the full-sequence transfer output BAM is readable') or diag($clip_error);
like($clip_view,
	qr/^clip\t0\t.*\t12M11S\t.*\tMM:Z:C\+m\.,0,2,2;\tML:B:C,201,202,250\tMN:i:23$/m,
	'soft clipping keeps the complete sequence and every original modification call');
like($clip_view,
	qr/^clip\t2048\t.*\t12S11M\t.*\tMM:Z:C\+m\.,0,2,2;\tML:B:C,201,202,250\tMN:i:23$/m,
	'full-sequence supplementary alignment receives the same unmodified read-relative tags');
like($clip_view,
	qr/^reverseclip\t16\t.*\t23M\t.*\tMM:Z:C\+m\.,0;\tML:B:C,177\tMN:i:23$/m,
	'reverse-strand acceptor also receives unchanged original-orientation tags');
like($clip_view,
	qr/^indel\t0\t.*\t6M1I6M1D10M\t.*\tMM:Z:C\+m\.,0;\tML:B:C,155\tMN:i:23$/m,
	'insertion/deletion alignment retains read-coordinate modification tags');
unlike($clip_view, qr/\tSA:Z:/, 'stale SA links are removed after alignment-set filtering');
like(slurp($clip_transfer_stats), qr/^output_alignments\t4$/m,
	'tag transfer publishes machine-readable post-policy record counts');

# Release validation only: modkit remains optional and is never a runtime
# dependency.  This fixture exercises forward, reverse, soft-clipped, split and
# indel-bearing records through a real default pileup without the legacy
# --force-allow-implicit compatibility switch.
SKIP: {
	my $modkit = $ENV{METH2REP_MODKIT} // '';
	skip 'set METH2REP_MODKIT=/path/to/modkit for the optional pileup interoperability test', 5
		unless $modkit ne '' && -x $modkit;
	my $clip_coordinate_bam = File::Spec->catfile($tmp, 'clip-transferred.coordinate.bam');
	run_ok('modkit fixture is coordinate sorted',
		$samtools, 'sort', '-o', $clip_coordinate_bam, $clip_transferred_bam);
	run_ok('modkit fixture is indexed', $samtools, 'index', $clip_coordinate_bam);
	my $clip_bed = File::Spec->catfile($tmp, 'clip-transferred.bed');
	my ($pileup_status, $pileup_stdout, $pileup_stderr) = capture(
		$modkit, 'pileup', $clip_coordinate_bam, $clip_bed,
		'--no-filtering', '--suppress-progress', '--threads', '1');
	is($pileup_status, 0,
		'real modkit pileup accepts normalized MM/ML across clipping, strand and indel cases')
		or diag($pileup_stdout, $pileup_stderr);
	ok(-s $clip_bed, 'real modkit pileup emits nonempty bedMethyl output');
	my @pileup_starts = sort { $a <=> $b } map { $_->[1] }
		grep { $_->[11] > 0 }
		map { [split /\t/] }
		grep { $_ ne '' && $_ !~ /^#/ } split /\n/, slurp($clip_bed);
	is_deeply(\@pileup_starts, [1, 8, 40, 80],
		'modkit projects forward, reverse, indel and aligned soft-clipped calls to expected reference coordinates')
		or diag(slurp($clip_bed));
}
my $hardclip_sam = File::Spec->catfile($tmp, 'hardclip-acceptor.sam');
my $hardclip_bam = File::Spec->catfile($tmp, 'hardclip-acceptor.bam');
write_file($hardclip_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t30\t60\t7H12M4H\t*\t0\t0\t$trimmed_sequence\t$trimmed_quality\tNM:i:0\n");
run_ok('hard-clipped acceptor fixture is created',
	$samtools, 'view', '-b', '-o', $hardclip_bam, $hardclip_sam);
my ($hardclip_status, undef, $hardclip_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $hardclip_bam, '--output', File::Spec->catfile($tmp, 'hardclip-rejected.bam'),
);
ok($hardclip_status != 0, 'hard-clipped acceptor fails instead of receiving guessed modification offsets');
like($hardclip_error, qr/hard-clipped/, 'hard-clipping failure identifies the unsafe condition');
my $edited_sequence = $clip_sequence;
substr($edited_sequence, -1, 1) = substr($edited_sequence, -1, 1) eq 'A' ? 'C' : 'A';
my $edited_sam = File::Spec->catfile($tmp, 'edited-acceptor.sam');
my $edited_bam = File::Spec->catfile($tmp, 'edited-acceptor.bam');
write_file($edited_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t1\t60\t23M\t*\t0\t0\t$edited_sequence\t$clip_quality\tNM:i:1\n");
run_ok('edited acceptor fixture is created', $samtools, 'view', '-b', '-o', $edited_bam, $edited_sam);
my ($edited_status, undef, $edited_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $edited_bam, '--output', File::Spec->catfile($tmp, 'edited-rejected.bam'),
);
ok($edited_status != 0, 'edited acceptor cannot inherit modification offsets');
like($edited_error, qr/differs from the complete original/, 'sequence mismatch is reported explicitly');
my $invalid_donor_sam = File::Spec->catfile($tmp, 'invalid-donor.sam');
my $invalid_donor_bam = File::Spec->catfile($tmp, 'invalid-donor.bam');
write_file($invalid_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m.,999;\tML:B:C,101\tMN:i:23\n");
run_ok('out-of-range MM donor fixture is created',
	$samtools, 'view', '-b', '-o', $invalid_donor_bam, $invalid_donor_sam);
my ($invalid_status, undef, $invalid_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $invalid_donor_bam,
	'--acceptor', $clip_acceptor_bam, '--output', File::Spec->catfile($tmp, 'invalid-rejected.bam'),
);
ok($invalid_status != 0, 'out-of-range MM coordinates cannot be published');
like($invalid_error, qr/MM call beyond/, 'out-of-range MM failure identifies the invalid donor');

# samtools -n uses natural numeric ordering and treats zero-padded numeric names
# as one equivalence class.  Primary/supplementary flag ordering can therefore
# interleave distinct exact QNAMEs inside that class.
my $numeric_sequence = 'ACCCGTTACGCC';
my $numeric_quality = 'I' x length($numeric_sequence);
my $numeric_donor_sam = File::Spec->catfile($tmp, 'numeric-donor.sam');
my $numeric_donor_unsorted = File::Spec->catfile($tmp, 'numeric-donor.unsorted.bam');
my $numeric_donor_bam = File::Spec->catfile($tmp, 'numeric-donor.name.bam');
write_file($numeric_donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:ref\tLN:100\n"
	. join('', map {
		"$_\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
		. "\tMM:Z:C+m.,0;\tML:B:C,120\tMN:i:12\n"
	} qw(read10 read002 read2 read02)));
run_ok('natural-QNAME donor fixture is created', $samtools, 'view', '-b',
	'-o', $numeric_donor_unsorted, $numeric_donor_sam);
run_ok('natural-QNAME donor fixture is name sorted', $samtools, 'sort', '-n',
	'-o', $numeric_donor_bam, $numeric_donor_unsorted);
my $numeric_acceptor_sam = File::Spec->catfile($tmp, 'numeric-acceptor.sam');
my $numeric_acceptor_unsorted = File::Spec->catfile($tmp, 'numeric-acceptor.unsorted.bam');
my $numeric_acceptor_bam = File::Spec->catfile($tmp, 'numeric-acceptor.name.bam');
write_file($numeric_acceptor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:ref\tLN:100\n"
	. "read2\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\tSA:Z:ref,30,+,12M,60,0;\n"
	. "read2\t2048\tref\t30\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read02\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read002\t2048\tref\t30\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read10\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n");
run_ok('natural-QNAME acceptor fixture is created', $samtools, 'view', '-b',
	'-o', $numeric_acceptor_unsorted, $numeric_acceptor_sam);
run_ok('natural-QNAME acceptor fixture is name sorted', $samtools, 'sort', '-n',
	'-o', $numeric_acceptor_bam, $numeric_acceptor_unsorted);
my $numeric_output = File::Spec->catfile($tmp, 'numeric-transferred.bam');
my $numeric_stats = File::Spec->catfile($tmp, 'numeric-transfer.tsv');
run_ok('natural-equivalence QNAME classes transfer without lexical merge errors',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $numeric_donor_bam, '--acceptor', $numeric_acceptor_bam,
	'--output', $numeric_output, '--stats', $numeric_stats);
my (undef, $numeric_view, undef) = capture($samtools, 'view', $numeric_output);
is_deeply([sort map { (split /\t/)[0] } grep { $_ ne '' } split /\n/, $numeric_view],
	[qw(read02 read10 read2)],
	'exact QNAMEs are joined correctly inside a zero-padded natural-name class');
unlike($numeric_view, qr/\tSA:Z:/, 'primary output cannot retain an SA link to a dropped supplementary');
like(slurp($numeric_stats), qr/^dropped_no_primary_groups\t1$/m,
	'an orphan supplementary-only exact QNAME group is counted and dropped');
like(slurp($numeric_stats), qr/^dropped_supplementary_alignments\t1$/m,
	'the safe default drops supplementary records from otherwise valid groups');

# Regression for a donor that is much larger than the accepted alignment set.
# The transfer must drain samtools' stdout instead of closing a still-writing
# pipe and turning a valid result into a dataset-size-dependent SIGPIPE error.
my $tail_donor_sam = File::Spec->catfile($tmp, 'large-tail-donor.sam');
my $tail_donor_bam = File::Spec->catfile($tmp, 'large-tail-donor.bam');
my $tail_acceptor_sam = File::Spec->catfile($tmp, 'large-tail-acceptor.sam');
my $tail_acceptor_bam = File::Spec->catfile($tmp, 'large-tail-acceptor.bam');
my $tail_records = '';
for my $index (1 .. 20000) {
	$tail_records .= "trail$index\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
		. "\tMM:Z:C+m.,0;\tML:B:C,100\tMN:i:12\n";
}
ok(length($tail_records) > 1_000_000,
	'trailing-donor regression exceeds ordinary pipe-buffer capacity');
write_file($tail_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "accepted\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+m.,0;\tML:B:C,120\tMN:i:12\n"
	. $tail_records);
write_file($tail_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "accepted\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n");
run_ok('large trailing donor BAM is created', $samtools, 'view', '-b',
	'-o', $tail_donor_bam, $tail_donor_sam);
run_ok('single-record large-tail acceptor BAM is created', $samtools, 'view', '-b',
	'-o', $tail_acceptor_bam, $tail_acceptor_sam);
run_ok('tag transfer drains a large unused donor tail without SIGPIPE',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $tail_donor_bam, '--acceptor', $tail_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'large-tail-transferred.bam'));

my $variant_donor_sam = File::Spec->catfile($tmp, 'variant-donor.sam');
my $variant_donor_bam = File::Spec->catfile($tmp, 'variant-donor.bam');
my $variant_acceptor_sam = File::Spec->catfile($tmp, 'variant-acceptor.sam');
my $variant_acceptor_bam = File::Spec->catfile($tmp, 'variant-acceptor.bam');
write_file($variant_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+mh?,0;A+123.,0;G+Z.,0;\tML:B:C,100,100,200,50\tMN:i:12\n");
write_file($variant_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n");
run_ok('multi-code and numeric-ChEBI donor is created', $samtools, 'view', '-b',
	'-o', $variant_donor_bam, $variant_donor_sam);
run_ok('multi-code acceptor is created', $samtools, 'view', '-b',
	'-o', $variant_acceptor_bam, $variant_acceptor_sam);
my $variant_output = File::Spec->catfile($tmp, 'variant-transferred.bam');
run_ok('formal MM alphabetic, numeric, question, and dot forms transfer intact',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $variant_donor_bam, '--acceptor', $variant_acceptor_bam,
	'--output', $variant_output);
my (undef, $variant_view, undef) = capture($samtools, 'view', $variant_output);
like($variant_view, qr/MM:Z:C\+mh\?,0;A\+123\.,0;G\+Z\.,0;/,
	'valid modern MM group syntax is preserved byte-for-byte');

my $legacy_tag_sam = File::Spec->catfile($tmp, 'legacy-tags.sam');
my $legacy_tag_bam = File::Spec->catfile($tmp, 'legacy-tags.bam');
write_file($legacy_tag_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMm:Z:C+m,0;\tMl:B:C,100\tMN:i:12\n");
run_ok('legacy-name and absent-mode donor is created', $samtools, 'view', '-b',
	'-o', $legacy_tag_bam, $legacy_tag_sam);
my $legacy_tag_output = File::Spec->catfile($tmp, 'legacy-tags-transferred.bam');
my $legacy_tag_stats = File::Spec->catfile($tmp, 'legacy-tags-transfer.tsv');
run_ok('legacy tag spelling and SAM-equivalent absent mode are normalized',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $legacy_tag_bam, '--acceptor', $variant_acceptor_bam,
	'--output', $legacy_tag_output, '--stats', $legacy_tag_stats);
my (undef, $legacy_tag_view, undef) = capture($samtools, 'view', $legacy_tag_output);
like($legacy_tag_view, qr/\tMM:Z:C\+m\.,0;\tML:B:C,100\tMN:i:12(?:\t|$)/,
	'new output uses standard MM/ML names and explicit equivalent dot mode');
unlike($legacy_tag_view, qr/\tM[ml]:/,
	'new output does not re-emit locally reserved draft tag names');
like(slurp($legacy_tag_stats), qr/^legacy_tag_names_normalized\t2$/m,
	'normalization statistics count both legacy tag names');
like(slurp($legacy_tag_stats), qr/^implicit_mode_groups_normalized\t1$/m,
	'normalization statistics count the formerly absent mode');

my $probability_donor_sam = File::Spec->catfile($tmp, 'probability-donor.sam');
my $probability_donor_bam = File::Spec->catfile($tmp, 'probability-donor.bam');
write_file($probability_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+mh.,0;\tML:B:C,200,200\tMN:i:12\n");
run_ok('invalid multi-code probability donor is created', $samtools, 'view', '-b',
	'-o', $probability_donor_bam, $probability_donor_sam);
my ($probability_status, undef, $probability_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $probability_donor_bam, '--acceptor', $variant_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'probability-rejected.bam'));
ok($probability_status != 0, 'mutually exclusive multi-code probabilities above one fail closed');
like($probability_error, qr/probabilities above 1/,
	'probability-sum failure identifies the inconsistent original read position');

my $missing_mn_sam = File::Spec->catfile($tmp, 'missing-mn.sam');
my $missing_mn_bam = File::Spec->catfile($tmp, 'missing-mn.bam');
write_file($missing_mn_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+m.,0;\tML:B:C,100\n");
run_ok('legacy donor without MN is created', $samtools, 'view', '-b',
	'-o', $missing_mn_bam, $missing_mn_sam);
my ($missing_mn_status, undef, $missing_mn_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $missing_mn_bam, '--acceptor', $variant_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'missing-mn-rejected.bam'));
ok($missing_mn_status != 0, 'MN is required by the fail-closed default');
like($missing_mn_error, qr/lacks MN/, 'missing-MN failure points to the legacy override');
my $legacy_output = File::Spec->catfile($tmp, 'missing-mn-allowed.bam');
my $legacy_stats = File::Spec->catfile($tmp, 'missing-mn-allowed.tsv');
run_ok('explicit legacy override accepts a reviewed donor without MN',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $missing_mn_bam, '--acceptor', $variant_acceptor_bam,
	'--output', $legacy_output, '--stats', $legacy_stats, '--allow-missing-mn');
my (undef, $legacy_view, undef) = capture($samtools, 'view', $legacy_output);
like($legacy_view, qr/\tMN:i:12(?:\t|$)/,
	'legacy override regenerates a current MN on the exact-sequence acceptor');
like(slurp($legacy_stats), qr/^legacy_missing_mn_reads\t1$/m,
	'legacy use is explicit in machine-readable transfer statistics');

my $overlap_acceptor_sam = File::Spec->catfile($tmp, 'overlap-acceptor.sam');
my $overlap_acceptor_bam = File::Spec->catfile($tmp, 'overlap-acceptor.bam');
write_file($overlap_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t1\t60\t5S13M5S\t*\t0\t0\t$clip_sequence\t$clip_quality\n"
	. "clip\t2048\tref\t30\t60\t7S12M4S\t*\t0\t0\t$clip_sequence\t$clip_quality\n");
run_ok('overlapping split-alignment fixture is created', $samtools, 'view', '-b',
	'-o', $overlap_acceptor_bam, $overlap_acceptor_sam);
my ($overlap_status, undef, $overlap_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $clip_donor_bam, '--acceptor', $overlap_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'overlap-rejected.bam'),
	'--supplementary', 'keep');
ok($overlap_status != 0, 'supplementary retention rejects overlapping native query spans');
like($overlap_error, qr/overlapping primary\/supplementary query spans/,
	'overlap failure explains the double-counting risk');

my ($split_left, $split_right) = ('', '');
for (1 .. 1800) {
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$split_left .= $bases[($seed >> 16) % 4];
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$split_right .= $bases[($seed >> 16) % 4];
}
my $split_read = substr($split_left, 300, 700) . substr($split_right, 500, 700);
my $split_reference = File::Spec->catfile($tmp, 'split-reference.fa');
my $split_fastq = File::Spec->catfile($tmp, 'split-read.fastq');
write_file($split_reference, ">left\n$split_left\n>right\n$split_right\n");
write_file($split_fastq,
	"\@split\n$split_read\n+\n" . ('I' x length($split_read)) . "\n");
my ($split_status, $split_sam, $split_error) = capture(
	$minimap2, '-a', '-Y', '--secondary=no', '-x', 'map-ont',
	$split_reference, $split_fastq,
);
is($split_status, 0, 'minimap2 creates the split-read clipping fixture')
	or diag($split_error);
my @split_records = map { [split /\t/, $_, -1] }
	grep { $_ ne '' && $_ !~ /^\@/ } split /\n/, $split_sam;
my @supplementary = grep { $_->[1] & 0x800 } @split_records;
ok(@supplementary, 'the chimeric read produces a supplementary alignment');
ok(!(grep { $_->[5] =~ /H/ } @supplementary),
	'-Y prevents hard clipping on every supplementary alignment');
ok(!(grep { length($_->[9]) != length($split_read) } @supplementary),
	'every supplementary alignment retains the complete read sequence for exact transfer');

done_testing();
use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3 qw(open3);
use JSON::PP qw(decode_json);
use Symbol qw(gensym);
use Test::More;
use lib File::Spec->catdir($Bin, '..');
use Mods::Checkpoint qw(write_checkpoint);

my $samtools = `command -v samtools 2>/dev/null`;
my $minimap2 = `command -v minimap2 2>/dev/null`;
chomp($samtools, $minimap2);
plan skip_all => 'samtools and minimap2 are required for the integration test'
	unless $samtools ne '' && $minimap2 ne '';

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub run_ok {
	my ($description, @command) = @_;
	my $status = system @command;
	is($status, 0, $description) or BAIL_OUT("command failed: @command");
}

sub capture {
	my (@command) = @_;
	my $stderr = gensym;
	my $pid = open3(undef, my $stdout, $stderr, @command);
	local $/;
	my $out = <$stdout> // '';
	my $err = <$stderr> // '';
	waitpid($pid, 0);
	return ($? >> 8, $out, $err);
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $contents = <$fh> // '';
	close $fh or die "Cannot close $path: $!";
	return $contents;
}

sub reverse_complement {
	my ($sequence) = @_;
	$sequence = reverse $sequence;
	$sequence =~ tr/ACGTNacgtn/TGCANtgcan/;
	return $sequence;
}

my $tmp = tempdir(CLEANUP => 1);
my $run_root = File::Spec->catdir($tmp, 'output', 'fixture');
my $sample_dir = File::Spec->catdir($run_root, 'S1');
my $mapping_dir = File::Spec->catdir($sample_dir, 'mapping');
my $assembly = File::Spec->catdir($tmp, 'assembly');
my $representatives = File::Spec->catdir($tmp, 'representatives');
my $mgs_dir = File::Spec->catdir($tmp, 'mgs', 'Bin_SB');
my $out = File::Spec->catdir($tmp, 'meth2rep');
my $out_manifest = File::Spec->catfile($tmp, 'tracking', 'meth2rep.tsv');
make_path(
	$mapping_dir,
	File::Spec->catdir($sample_dir, 'assemblies', 'metag'),
	File::Spec->catdir($assembly, 'Binning', 'SB'),
	$representatives,
	File::Spec->catdir($mgs_dir, 'LOGandSUB', 'checkpoints'),
	File::Spec->catdir($tmp, 'raw'),
);

my $seed = 17;
my @bases = qw(A C G T);
my $reference_sequence = '';
for (1 .. 1200) {
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$reference_sequence .= $bases[($seed >> 16) % 4];
}
my $read_sequence = substr($reference_sequence, 200, 600);
$read_sequence =~ s/^[^C]*//; # ensure MM delta zero addresses a canonical C
my $read_length = length($read_sequence);
my $quality = 'I' x $read_length;
my $reverse_donor_sequence = reverse_complement($read_sequence);
my $member_read_sequence = substr($reference_sequence, 400, 600);
$member_read_sequence =~ s/^[^C]*//;
my $member_read_length = length($member_read_sequence);
my $member_quality = 'I' x $member_read_length;
my $short_source_cigar = '50S10M' . ($read_length - 60) . 'S';
my $assembly_fasta = File::Spec->catfile($assembly, 'scaffolds.fasta.filt');
write_file($assembly_fasta,
	">S1-contig\n$reference_sequence\n>S1-member-contig\n$reference_sequence\n");
write_file(File::Spec->catfile($sample_dir, 'assemblies', 'metag', 'assembly.txt'), "$assembly\n");
write_file(File::Spec->catfile($assembly, 'Binning', 'SB', 'S1'),
	"Sequence ID\tBin\nS1-contig\t1.fa.gz\nS1-member-contig\t2.fa.gz\n");

my $representative = File::Spec->catfile($representatives, 'MGS.1.ctgs.S1__1.fa.gz');
gzip(\(">S1-contig\n$reference_sequence\n") => $representative)
	or die "Cannot write $representative: $GzipError";
my $report = File::Spec->catfile($tmp, 'MAGvsGC.txt.gz');
gzip(\("MAG\tMGS\tRepresentative4MGS\n"
		. "S1__1.fa.gz\tMGS.1\t*\nS1__2.fa.gz\tMGS.1\t\n"
		. "Cano__7\tMGS.2\t*\n") => $report)
	or die "Cannot write $report: $GzipError";

my $donor_sam = File::Spec->catfile($tmp, 'donor.sam');
my $donor_bam = File::Spec->catfile($tmp, 'donor.bam');
my $donor2_sam = File::Spec->catfile($tmp, 'donor2.sam');
my $donor2_bam = File::Spec->catfile($tmp, 'donor2.bam');
write_file($donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:donor-ref\tLN:1200\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$read_sequence\t$quality\tMM:Z:C+m,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r3\t16\tdonor-ref\t201\t60\t${read_length}M\t*\t0\t0\t$reverse_donor_sequence\t$quality\tMM:Z:C+m,0;\tML:B:C,210\tMN:i:$read_length\tNM:i:0\n");
run_ok('synthetic donor modBAM is created',
	$samtools, 'view', '-b', '-o', $donor_bam, $donor_sam);
write_file($donor2_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n"
	. "r2\t4\t*\t0\t0\t*\t*\t0\t0\t$member_read_sequence\t$member_quality\tMM:Z:C+m,0;\tML:B:C,180\tMN:i:$member_read_length\n");
run_ok('second original modBAM is created',
	$samtools, 'view', '-b', '-o', $donor2_bam, $donor2_sam);

my $candidate_sam = File::Spec->catfile($tmp, 'candidate.sam');
my $candidate_cram = File::Spec->catfile($mapping_dir, 'S1-smd.cram');
write_file($candidate_sam,
	"\@HD\tVN:1.6\tSO:coordinate\n\@SQ\tSN:S1-contig\tLN:1200\n\@SQ\tSN:S1-member-contig\tLN:1200\n"
	. "r1\t0\tS1-contig\t201\t60\t${read_length}M\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n"
	. "r2\t0\tS1-member-contig\t401\t60\t${member_read_length}M\t*\t0\t0\t$member_read_sequence\t$member_quality\tNM:i:0\n"
	. "r3\t0\tS1-contig\t201\t60\t${read_length}M\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n"
	. "r4\t0\tS1-contig\t201\t60\t$short_source_cigar\t*\t0\t0\t$read_sequence\t$quality\tNM:i:0\n");
run_ok('synthetic assembly backmapping CRAM is created',
	$samtools, 'view', '-C', '-T', $assembly_fasta, '-o', $candidate_cram, $candidate_sam);
write_file("$candidate_cram.sto", "done\n");
my @assembly_stat = stat($assembly_fasta);
write_file(File::Spec->catfile($mapping_dir, 'S1-smd.reference.stat'),
	"$assembly_stat[7] $assembly_stat[9]\n");

my $map = File::Spec->catfile($tmp, 'fixture.map');
write_file($map,
	"#SmplID\tPath\tAssmblGrps\tSupportReads\tINFO\tSeqTech\n"
	. "#RunID\tfixture\n#OutPath\t" . File::Spec->catdir($tmp, 'output') . "/\n"
	. "#DirPath\t" . File::Spec->catdir($tmp, 'raw') . "/\n"
	. "#WARNING\tOFF\nS1\tinput\tG1\t\t\tONT\n");
my $manifest = File::Spec->catfile($tmp, 'modbams.tsv');
write_file($manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$donor_bam\n"
	. "S1\tprimary\tONT\t$donor2_bam\n");
write_file(File::Spec->catfile($sample_dir, 'input_raw.txt'), "$donor_bam;$donor2_bam");
write_checkpoint(File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'Stage1.stone'),
	parameters => { stage => 'stage-1' }, outputs => [$report]);
my $checkpoint_representative = $representative;
$checkpoint_representative =~ s{/representatives/}{/representatives//};
write_checkpoint(File::Spec->catfile($mgs_dir, 'LOGandSUB', 'checkpoints', 'BinExtr.stone'),
	parameters => { stage => 'extract-bin-contigs' }, outputs => [$report, $checkpoint_representative]);

my $script = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'meth2rep.pl');
my $bam_filter = 'perl ' . File::Spec->catfile($Bin, '..', 'secScripts', 'assemblies', 'bamFilter.pl');
sub output_bam {
	my ($out_dir, $mode) = @_;
	return File::Spec->catfile($out_dir, 'MGS.1',
		join('__', 'S1', $mode, 'S1__1.fa.gz') . '.mod.bam');
}
sub output_alignment {
	my ($out_dir, $mode, $format) = @_;
	return File::Spec->catfile($out_dir, 'MGS.1',
		join('__', 'S1', $mode, 'S1__1.fa.gz') . ".mod.$format");
}

sub alignment_signatures {
	my ($sam) = @_;
	my @signatures;
	for my $line (grep { $_ ne '' && $_ !~ /^\@/ } split /\n/, $sam) {
		my @fields = split /\t/, $line, -1;
		my %tags;
		for my $field (@fields[11 .. $#fields]) {
			$tags{MM} = $field if $field =~ /^M[Mm]:/;
			$tags{ML} = $field if $field =~ /^M[Ll]:/;
			$tags{MN} = $field if $field =~ /^MN:/;
		}
		push @signatures, join("\t", @fields[0 .. 5], $fields[9], @tags{qw(MM ML MN)});
	}
	return [sort @signatures];
}
my @command = (
	$^X, '-I' . File::Spec->catdir($Bin, '..'), $script,
	'--mgs-dir', $mgs_dir,
	'--map', $map, '--mgs-report', $report,
	'--representatives-dir', $representatives, '--binner', 'SB',
	'--modbam-manifest', $manifest, '--out', $out, '--mgs', 'MGS.1',
	'--out-manifest', $out_manifest,
	'--mgs2rep', '--rep2rep', '--threads', 2,
	'--keep-read-ids',
	'--samtools', $samtools, '--minimap2', $minimap2,
	'--bam-filter', $bam_filter,
);
my ($status, $stdout, $stderr) = capture(@command);
is($status, 0, 'meth2rep completes through real CRAM/minimap2 and internal tag transfer')
	or diag($stdout, $stderr);
ok(-s $out_manifest, 'a custom completed-unit donor manifest is published');
my $tracking = slurp($out_manifest);
like($tracking, qr/^MGS\.1\tmgs2rep\tS1\tprimary\tONT\t\Q$donor_bam\E\t2\t3\t3\t3\t3\tcomplete\t/m,
	'the first donor has two attributed MGS candidate reads');
like($tracking, qr/^MGS\.1\tmgs2rep\tS1\tprimary\tONT\t\Q$donor2_bam\E\t1\t3\t3\t3\t3\tcomplete\t/m,
	'the second donor has one attributed MGS candidate read');
like($tracking, qr/^MGS\.1\trep2rep\tS1\tprimary\tONT\t\Q$donor2_bam\E\t0\t2\t2\t2\t2\tcomplete\t/m,
	'the donor manifest retains zero-use donor provenance per requested mode');
my $sample_log_path = File::Spec->catfile($out, 'MGS.1', 'S1.meth2rep.json');
ok(-s $sample_log_path, 'one sample-level interrogation log is published inside the MGS directory');
my $sample_log = decode_json(slurp($sample_log_path));
is_deeply([sort keys %{$sample_log->{modes}}], [qw(mgs2rep rep2rep)],
	'the single sample log reports both requested modes');
is($sample_log->{modes}{mgs2rep}{coverage}{reference_bases}, 1200,
	'the sample log records the representative reference denominator');
ok($sample_log->{modes}{mgs2rep}{coverage}{covered_bases} > 0
	&& $sample_log->{modes}{mgs2rep}{coverage}{breadth_fraction} > 0
	&& $sample_log->{modes}{mgs2rep}{coverage}{mean_depth} > 0,
	'the sample log records covered bases, breadth, and mean depth');
is($sample_log->{modes}{mgs2rep}{scopes}[0]{filter_stats}{malformed}, 0,
	'the per-scope alignment-filter diagnostics are retained');
is_deeply($sample_log->{pipeline_recorded_inputs}, [$donor_bam, $donor2_bam],
	'pipeline input_raw provenance is cross-checked and exposed beside the donor manifest');
ok(!(grep { /absent from MATAFILER/ } @{$sample_log->{warnings}}),
	'matching recorded inputs produce no donor-provenance warning');
opendir my $mgs_dir_handle, File::Spec->catdir($out, 'MGS.1') or die $!;
my @mgs_files = sort grep { $_ ne '.' && $_ ne '..' } readdir $mgs_dir_handle;
closedir $mgs_dir_handle;
is(scalar(@mgs_files), 5, 'the MGS folder contains two mode BAMs, their indexes, and one sample log');
my $origins_gzip = File::Spec->catfile($out, '.meth2rep', 'read_ids', 'mgs2rep', 'MGS.1', 'S1.read_origins.tsv.gz');
my $origins = '';
gunzip($origins_gzip => \$origins) or die "Cannot read $origins_gzip: $GunzipError";
like($origins, qr/^S1\tprimary\tr1\t\Q$donor_bam\E$/m,
	'compressed per-read origins identify the first donor');
like($origins, qr/^S1\tprimary\tr2\t\Q$donor2_bam\E$/m,
	'compressed per-read origins identify the second donor');
my $initial_inode = (stat(output_bam($out, 'mgs2rep')))[1];

for my $mode (qw(mgs2rep rep2rep)) {
	my $bam = output_bam($out, $mode);
	if (!-s $bam) {
		my $summary_file = File::Spec->catfile($out, '.meth2rep', 'summary.tsv');
		my $summary_text = '';
		if (-s $summary_file) {
			open my $summary_fh, '<', $summary_file or die "Cannot read $summary_file: $!";
			{ local $/; $summary_text = <$summary_fh> // ''; }
			close $summary_fh;
		}
		diag("meth2rep stdout:\n$stdout\nstderr:\n$stderr\nsummary:\n$summary_text");
	}
	ok(-s $bam && -s "$bam.bai", "$mode publishes an indexed compact modBAM");
	my ($view_status, $view, $view_error) = capture($samtools, 'view', $bam);
	is($view_status, 0, "$mode modBAM is readable") or diag($view_error);
	my @records = grep { $_ ne '' } split /\n/, $view;
	is(scalar(@records), $mode eq 'mgs2rep' ? 3 : 2,
		"$mode contains exactly its intended all-member or representative-only read set");
	like($view, qr/^r1\t/m, "$mode retains the representative-MAG read");
	like($view, qr/^r3\t/m,
		"$mode recovers the original orientation from a reverse-aligned donor");
	if ($mode eq 'mgs2rep') {
		like($view, qr/^r2\t/m, 'mgs2rep includes the non-representative member-MAG read');
	} else {
		unlike($view, qr/^r2\t/m, 'rep2rep excludes the non-representative member-MAG read');
	}
	like($view, qr/\tMM:Z:C\+m,0;/, "$mode retains the MM methylation tag");
	like($view, qr/\tML:B:C,220(?:\t|\n)/, "$mode retains the ML probability tag");
	like($view, qr/^r3\t.*\tML:B:C,210(?:\t|$)/m,
		"$mode preserves the reverse-aligned donor's modification probability");
}

my $cram_out = File::Spec->catdir($tmp, 'cram-output');
my @cram_command = grep { $_ ne '--rep2rep' } @command;
for my $i (0 .. $#cram_command - 1) {
	$cram_command[$i + 1] = $cram_out if $cram_command[$i] eq '--out';
	$cram_command[$i + 1] = File::Spec->catfile($cram_out, 'manifest.tsv')
		if $cram_command[$i] eq '--out-manifest';
}
push @cram_command, '--output-format', 'cram';
my ($cram_run_status, $cram_run_stdout, $cram_run_stderr) = capture(@cram_command);
is($cram_run_status, 0, 'self-contained CRAM output completes from an ordinary gzip representative')
	or diag($cram_run_stdout, $cram_run_stderr);
my $cram = output_alignment($cram_out, 'mgs2rep', 'cram');
ok(-s $cram && -s "$cram.crai", 'CRAM output and CRAI use explicit predictable suffixes');
ok(!-e output_alignment($cram_out, 'mgs2rep', 'bam'),
	'CRAM selection does not leave a duplicate BAM output');
ok(!-e "$representative.fai" && !-e "$representative.gzi",
	'CRAM encoding does not write indexes beside the immutable gzip representative');
my ($cram_view_status, $cram_view, $cram_view_error) = capture($samtools, 'view', $cram);
is($cram_view_status, 0, 'embedded-reference CRAM decodes without an external -T reference')
	or diag($cram_view_error);
my (undef, $bam_view_for_cram, undef) = capture($samtools, 'view', output_bam($out, 'mgs2rep'));
is_deeply(alignment_signatures($cram_view), alignment_signatures($bam_view_for_cram),
	'BAM and CRAM preserve identical alignment fields and MM/ML/MN payloads');
my ($idxstats_status, undef, $idxstats_error) = capture($samtools, 'idxstats', $cram);
is($idxstats_status, 0, 'published CRAI is usable by samtools idxstats') or diag($idxstats_error);
my $cram_log = decode_json(slurp(File::Spec->catfile($cram_out, 'MGS.1', 'S1.meth2rep.json')));
is($cram_log->{modes}{mgs2rep}{alignment_format}, 'cram',
	'sample log identifies the actual compact output format');
is($cram_log->{modes}{mgs2rep}{index}, "$cram.crai",
	'sample log records the exact CRAI path');
my $cram_inode = (stat($cram))[1];
my ($cram_resume_status, undef, $cram_resume_error) = capture(@cram_command);
is($cram_resume_status, 0, 'an identical CRAM invocation resumes') or diag($cram_resume_error);
is((stat($cram))[1], $cram_inode, 'CRAM resume does not rewrite a validated alignment');
my ($format_switch_status, undef, $format_switch_error) = capture(
	@cram_command, '--output-format', 'bam');
ok($format_switch_status != 0, 'changing output format cannot silently leave BAM and CRAM siblings');
like($format_switch_error, qr/alternate-format meth2rep output already exists/,
	'format-switch refusal explains the required override');
my ($format_override_status, $format_override_stdout, $format_override_error) = capture(
	@cram_command, '--output-format', 'bam', '--override');
is($format_override_status, 0, 'explicit override safely replaces CRAM with BAM')
	or diag($format_override_stdout, $format_override_error);
ok(-s output_alignment($cram_out, 'mgs2rep', 'bam')
	&& !-e $cram && !-e "$cram.crai",
	'format replacement removes the prior tracked CRAM and CRAI only after BAM publication');

my ($resume_status, $resume_stdout, $resume_stderr) = capture(@command);
is($resume_status, 0, 'a second identical invocation resumes successfully')
	or diag($resume_stderr);
like($resume_stdout, qr/completed/,
	'identical input/options retain completed per-unit results');
is((stat(output_bam($out, 'mgs2rep')))[1], $initial_inode,
	'a valid per-unit resume does not remap or replace an existing modBAM');

my $checkpoint = File::Spec->catfile($out, '.meth2rep', 'complete.stone');
my $checkpoint_before_preview = slurp($checkpoint);
my ($preview_status, $preview_stdout, $preview_stderr) = capture(@command, '--plan-only');
is($preview_status, 0, 'plan-only validates without executing the mapper')
	or diag($preview_stdout, $preview_stderr);
my $preview = File::Spec->catfile($out, '.meth2rep', 'plan.preview.tsv');
ok(-s $preview, 'plan-only writes a separate preview plan');
is(slurp($checkpoint), $checkpoint_before_preview,
	'plan-only does not invalidate or rewrite a completed checkpoint');
my @preview_lines = grep { $_ ne '' } split /\n/, slurp($preview);
ok(!(grep { scalar(split /\t/, $_, -1) != 12 } @preview_lines),
	'preview plan remains a regular twelve-column TSV when sample scopes are embedded');
my ($invalid_preset_status, undef, $invalid_preset_error) = capture(
	@command, '--minimap2-preset-ont', 'definitely-invalid-preset', '--plan-only');
ok($invalid_preset_status != 0,
	'plan-only rejects a minimap2 preset unsupported by the installed executable');
like($invalid_preset_error, qr/(?:unknown preset|validating minimap2 preset).*definitely-invalid-preset/is,
	'unsupported-preset preflight identifies the requested preset');
my @path_discovery_command;
for (my $i = 0; $i <= $#command; $i++) {
	if ($command[$i] eq '--samtools' || $command[$i] eq '--minimap2'
		|| $command[$i] eq '--bam-filter') {
		$i++;
		next;
	}
	push @path_discovery_command, $command[$i];
}
for my $i (0 .. $#path_discovery_command - 1) {
	$path_discovery_command[$i + 1] = File::Spec->catdir($tmp, 'path-preflight')
		if $path_discovery_command[$i] eq '--out';
	$path_discovery_command[$i + 1] = File::Spec->catfile($tmp, 'path-preflight', 'manifest.tsv')
		if $path_discovery_command[$i] eq '--out-manifest';
}
my ($path_status, undef, $path_error);
{
	local $ENV{MGTKDIR};
	delete $ENV{MGTKDIR};
	($path_status, undef, $path_error) = capture(@path_discovery_command, '--plan-only');
}
is($path_status, 0,
	'standalone preflight resolves samtools/minimap2 from PATH and the co-located bamFilter without MGTKDIR')
	or diag($path_error);
my ($manifest_collision_status, undef, $manifest_collision_error) = capture(
	@command, '--out-manifest', $sample_log_path, '--plan-only');
ok($manifest_collision_status != 0,
	'custom manifest cannot overwrite a per-MGS alignment, index, or sample log');
like($manifest_collision_error, qr/cannot be placed inside a per-MGS output directory/,
	'manifest collision fails during preflight with the unsafe path');
my $historical_dir = File::Spec->catdir($out, 'MGS.historical');
make_path($historical_dir);
my $historical_alignment = File::Spec->catfile($historical_dir, 'prior.mod.bam');
write_file($historical_alignment, "historical-alignment-data\n");
my $historical_before = slurp($historical_alignment);
my $historical_inode = (stat($historical_alignment))[1];
my ($historical_collision_status, undef, $historical_collision_error) = capture(
	@command, '--out-manifest', $historical_alignment, '--plan-only');
ok($historical_collision_status != 0,
	'a custom manifest cannot target an unselected historical MGS data directory');
like($historical_collision_error, qr/cannot be placed inside a per-MGS output directory/,
	'historical-data collision is rejected during preflight');
is(slurp($historical_alignment), $historical_before,
	'rejected manifest placement leaves historical MGS data unchanged');
is((stat($historical_alignment))[1], $historical_inode,
	'rejected manifest placement does not replace the historical file');
my @short_out_command = @command;
for my $i (0 .. $#short_out_command - 1) {
	$short_out_command[$i] = '-o' if $short_out_command[$i] eq '--out';
	$short_out_command[$i + 1] = File::Spec->catdir($tmp, 'short-out')
		if $short_out_command[$i] eq '-o';
	$short_out_command[$i + 1] = File::Spec->catfile($tmp, 'short-out', 'manifest.tsv')
		if $short_out_command[$i] eq '--out-manifest';
}
my ($short_out_status, undef, $short_out_error) = capture(@short_out_command, '--plan-only');
is($short_out_status, 0, 'the single short -o output option is accepted')
	or diag($short_out_error);

my @mgs_only_command = grep { $_ ne '--rep2rep' } @command;
my ($mgs_only_status, $mgs_only_stdout, $mgs_only_stderr) = capture(@mgs_only_command);
is($mgs_only_status, 0, 'a narrower mode rerun completes')
	or diag($mgs_only_stdout, $mgs_only_stderr);
ok(-s output_bam($out, 'mgs2rep'),
	'the still-requested mgs2rep output remains published');
ok(-s output_bam($out, 'rep2rep'),
	'a narrower invocation preserves previously completed unselected modes');
like(slurp($out_manifest), qr/^MGS\.1\trep2rep\tS1\t/m,
	'the tracking manifest retains previously completed unselected modes');

my $unready_dir = File::Spec->catdir($tmp, 'unready', 'Bin_SB');
make_path($unready_dir);
my @unready_command = @command;
for my $i (0 .. $#unready_command - 1) {
	$unready_command[$i + 1] = $unready_dir if $unready_command[$i] eq '--mgs-dir';
}
my ($unready_status, undef, $unready_error) = capture(@unready_command, '--plan-only');
ok($unready_status != 0, 'missing MGS progress checkpoints block even a plan preview');
like($unready_error, qr/progress checkpoint is missing/, 'readiness failure names the missing checkpoint');
my $stale_report = File::Spec->catfile($tmp, 'post-checkpoint-MAGvsGC.txt.gz');
write_file($stale_report, slurp($report));
utime(time + 30, time + 30, $stale_report) or die "Cannot age $stale_report: $!";
my @stale_report_command = @command;
for my $i (0 .. $#stale_report_command - 1) {
	$stale_report_command[$i + 1] = $stale_report if $stale_report_command[$i] eq '--mgs-report';
}
my ($stale_report_status, undef, $stale_report_error) = capture(@stale_report_command, '--plan-only');
ok($stale_report_status != 0, 'a membership report newer than Stage1 completion is rejected');
like($stale_report_error, qr/membership report was modified after/, 'the stale-report error explains the provenance mismatch');

my $collision_sam = File::Spec->catfile($tmp, 'collision.sam');
my $collision_bam = File::Spec->catfile($tmp, 'collision.bam');
write_file($collision_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$read_sequence\t$quality\tMM:Z:C+m,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r2\t4\t*\t0\t0\t*\t*\t0\t0\t$member_read_sequence\t$member_quality\tMM:Z:C+m,0;\tML:B:C,180\tMN:i:$member_read_length\n");
run_ok('colliding donor modBAM is created', $samtools, 'view', '-b', '-o', $collision_bam, $collision_sam);
my $collision_manifest = File::Spec->catfile($tmp, 'collision-modbams.tsv');
write_file($collision_manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$donor_bam\n"
	. "S1\tprimary\tONT\t$collision_bam\n");
my @collision_command = @command;
for my $i (0 .. $#collision_command - 1) {
	$collision_command[$i + 1] = $collision_manifest if $collision_command[$i] eq '--modbam-manifest';
	$collision_command[$i + 1] = File::Spec->catdir($tmp, 'collision-out') if $collision_command[$i] eq '--out';
}
my ($collision_status, undef, $collision_error) = capture(@collision_command);
ok($collision_status != 0, 'a QNAME shared across original modBAMs fails closed');
like($collision_error, qr/more than one primary record named 'r1'/,
	'the collision error identifies the ambiguous read');
my $wrong_sequence = $read_sequence;
substr($wrong_sequence, -1, 1) = substr($wrong_sequence, -1, 1) eq 'A' ? 'C' : 'A';
my $wrong_donor_sam = File::Spec->catfile($tmp, 'wrong-origin.sam');
my $wrong_donor_bam = File::Spec->catfile($tmp, 'wrong-origin.bam');
write_file($wrong_donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:donor-ref\tLN:1200\n"
	. "r1\t4\t*\t0\t0\t*\t*\t0\t0\t$wrong_sequence\t$quality\tMM:Z:C+m,0;\tML:B:C,220\tMN:i:$read_length\n"
	. "r3\t16\tdonor-ref\t201\t60\t${read_length}M\t*\t0\t0\t$reverse_donor_sequence\t$quality\tMM:Z:C+m,0;\tML:B:C,210\tMN:i:$read_length\tNM:i:0\n");
run_ok('same-QNAME wrong-sequence donor is created', $samtools, 'view', '-b',
	'-o', $wrong_donor_bam, $wrong_donor_sam);
my $wrong_manifest = File::Spec->catfile($tmp, 'wrong-origin.tsv');
write_file($wrong_manifest,
	"sample\tscope\ttechnology\tmodbam\nS1\tprimary\tONT\t$wrong_donor_bam\n");
my @wrong_origin_command = grep { $_ ne '--mgs2rep' } @command;
for my $i (0 .. $#wrong_origin_command - 1) {
	$wrong_origin_command[$i + 1] = $wrong_manifest
		if $wrong_origin_command[$i] eq '--modbam-manifest';
	$wrong_origin_command[$i + 1] = File::Spec->catdir($tmp, 'wrong-origin-out')
		if $wrong_origin_command[$i] eq '--out';
	$wrong_origin_command[$i + 1] = File::Spec->catfile($tmp, 'wrong-origin-out', 'manifest.tsv')
		if $wrong_origin_command[$i] eq '--out-manifest';
}
my ($wrong_origin_status, undef, $wrong_origin_error) = capture(@wrong_origin_command);
ok($wrong_origin_status != 0,
	'a coincidentally matching QNAME in the wrong donor cannot acquire methylation tags');
like($wrong_origin_error, qr/native SEQ.*differs between the assembly CRAM and declared original modBAM/,
	'wrong-donor failure requires sequence provenance in addition to QNAME');
my $untracked_out = File::Spec->catdir($tmp, 'untracked-out');
make_path(File::Spec->catdir($untracked_out, 'MGS.1'));
my $untracked_bam = output_bam($untracked_out, 'mgs2rep');
write_file($untracked_bam, 'user-owned fixture');
my @untracked_command = @command;
for my $i (0 .. $#untracked_command - 1) {
	$untracked_command[$i + 1] = $untracked_out if $untracked_command[$i] eq '--out';
	$untracked_command[$i + 1] = File::Spec->catfile($untracked_out, 'manifest.tsv')
		if $untracked_command[$i] eq '--out-manifest';
}
my ($untracked_status, undef, $untracked_error) = capture(@untracked_command);
ok($untracked_status != 0, 'an existing uncheckpointed BAM is not overwritten');
like($untracked_error, qr/Untracked output alignment or index already exists/, 'the overwrite refusal identifies the exact condition');
is(slurp($untracked_bam), 'user-owned fixture', 'the pre-existing file remains unchanged');
my ($invalid_coverage_status, undef, $invalid_coverage_error) = capture(
	@command, '--source-min-coverage', 1.2, '--plan-only');
ok($invalid_coverage_status != 0, 'an invalid source-coverage threshold is rejected');
like($invalid_coverage_error, qr/source-min-coverage must be between 0 and 1/,
	'coverage validation names the offending option');
my @no_source_floor_command = @command;
for my $i (0 .. $#no_source_floor_command - 1) {
	$no_source_floor_command[$i + 1] = File::Spec->catdir($tmp, 'no-source-floor')
		if $no_source_floor_command[$i] eq '--out';
	$no_source_floor_command[$i + 1] = File::Spec->catfile($tmp, 'no-source-floor', 'manifest.tsv')
		if $no_source_floor_command[$i] eq '--out-manifest';
}
my ($no_floor_status, undef, $no_floor_error) = capture(
	@no_source_floor_command, '--source-min-coverage', 0);
ok($no_floor_status != 0, 'disabling source coverage exposes a low-overlap read absent from the donors');
like($no_floor_error, qr/candidate read name\(s\) are absent.*r4/s,
	'the source-coverage failure identifies the otherwise excluded read');
my @high_mapq_command = @command;
for my $i (0 .. $#high_mapq_command - 1) {
	$high_mapq_command[$i + 1] = File::Spec->catdir($tmp, 'high-mapq')
		if $high_mapq_command[$i] eq '--out';
	$high_mapq_command[$i + 1] = File::Spec->catfile($tmp, 'high-mapq', 'manifest.tsv')
		if $high_mapq_command[$i] eq '--out-manifest';
}
my ($high_mapq_status, $high_mapq_stdout, $high_mapq_error) = capture(
	@high_mapq_command, '--target-min-mapq', 200, '--target-min-coverage', 0.7,
	'--target-max-edit-rate', 0.2, '--target-min-end-clip', 5);
is($high_mapq_status, 0, 'named target MAPQ, coverage, edit-rate, and clipping controls execute end to end')
	or diag($high_mapq_stdout, $high_mapq_error);
like(slurp(File::Spec->catfile($tmp, 'high-mapq', '.meth2rep', 'summary.tsv')),
	qr/^MGS\.1\tmgs2rep\tS1\tprimary\tno_target_alignment\t3\t3\t0\t/m,
	'a high MAPQ floor rejects every representative alignment while preserving donor accounting');
my $high_mapq_log = decode_json(slurp(File::Spec->catfile(
	$tmp, 'high-mapq', 'MGS.1', 'S1.meth2rep.json')));
is_deeply($high_mapq_log->{filters}{ont_target}, [0.2, 0.7, 200, 5],
	'sample log records the effective named target-filter overrides in bamFilter order');
my ($override_status, $override_stdout, $override_error) = capture(@command, '--override');
is($override_status, 0, '--override rebuilds previously complete selected units')
	or diag($override_stdout, $override_error);
isnt((stat(output_bam($out, 'mgs2rep')))[1], $initial_inode,
	'override replaces the selected output rather than reusing its old inode');

my $canopy_out = File::Spec->catdir($tmp, 'canopy-only');
my @canopy_command = (
	$^X, '-I' . File::Spec->catdir($Bin, '..'), $script,
	'--mgs-dir', $mgs_dir,
	'--map', $map, '--mgs-report', $report,
	'--representatives-dir', $representatives, '--binner', 'SB',
	'--modbam-manifest', $manifest, '--out', $canopy_out,
	'--mgs', 'MGS.2', '--rep2rep',
	'--samtools', '/not/required/samtools', '--minimap2', '/not/required/minimap2',
	'--bam-filter', '/not/required/bamFilter',
);
my ($canopy_status, $canopy_stdout, $canopy_stderr) = capture(@canopy_command);
is($canopy_status, 0, 'a Canopy-only selection reports biological unavailability without mapping tools')
	or diag($canopy_stdout, $canopy_stderr);
like(slurp(File::Spec->catfile($canopy_out, '.meth2rep', 'summary.tsv')),
	qr/^MGS\.2\trep2rep\t-\t-\tcanopy_only_no_mag_reference\t/m,
	'Canopy-only output carries the explicit unavailable state');
ok(-s File::Spec->catfile($canopy_out, '.meth2rep', 'complete.stone'),
	'Canopy-only unavailability is a checkpointed successful result');

# Exercise our deliberately narrower transfer contract independently of minimap2.
# Soft clipping preserves the complete SEQ, while hard clipping must fail closed.
my $clip_sequence = 'ACCCGTTACGCCATGCTAACCGT';
my $clip_quality = 'I' x length($clip_sequence);
my $trimmed_sequence = substr($clip_sequence, 7, 12);
my $trimmed_quality = 'I' x length($trimmed_sequence);
my $clip_donor_sam = File::Spec->catfile($tmp, 'clip-donor.sam');
my $clip_donor_bam = File::Spec->catfile($tmp, 'clip-donor.bam');
write_file($clip_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m,0,2,2;\tML:B:C,101,202,250\tMN:i:" . length($clip_sequence) . "\n"
	. "reverseclip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m,0;\tML:B:C,77\tMN:i:" . length($clip_sequence) . "\n");
run_ok('clipping-contract donor BAM is created',
	$samtools, 'view', '-b', '-o', $clip_donor_bam, $clip_donor_sam);
my $clip_acceptor_sam = File::Spec->catfile($tmp, 'clip-acceptor.sam');
my $clip_acceptor_bam = File::Spec->catfile($tmp, 'clip-acceptor.bam');
write_file($clip_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t1\t60\t12M11S\t*\t0\t0\t$clip_sequence\t$clip_quality\tNM:i:0\tSA:Z:ref,30,+,12S11M,60,0;\n"
	. "clip\t2048\tref\t30\t60\t12S11M\t*\t0\t0\t$clip_sequence\t$clip_quality\tNM:i:0\tSA:Z:ref,1,+,12M11S,60,0;\n"
	. "reverseclip\t16\tref\t60\t60\t23M\t*\t0\t0\t" . reverse_complement($clip_sequence)
	. "\t$clip_quality\tNM:i:0\n");
run_ok('clipping-contract acceptor BAM is created',
	$samtools, 'view', '-b', '-o', $clip_acceptor_bam, $clip_acceptor_sam);
my $transfer_script = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'transfer_mod_tags.pl');
my $clip_transferred_bam = File::Spec->catfile($tmp, 'clip-transferred.bam');
my $clip_transfer_stats = File::Spec->catfile($tmp, 'clip-transfer.tsv');
run_ok('internal transfer accepts full-sequence soft-clipped, supplementary, and reverse records',
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $clip_acceptor_bam, '--output', $clip_transferred_bam,
	'--supplementary', 'keep', '--stats', $clip_transfer_stats);
my ($clip_status, $clip_view, $clip_error) = capture($samtools, 'view', $clip_transferred_bam);
is($clip_status, 0, 'the full-sequence transfer output BAM is readable') or diag($clip_error);
like($clip_view,
	qr/^clip\t0\t.*\t12M11S\t.*\tMM:Z:C\+m,0,2,2;\tML:B:C,101,202,250\tMN:i:23$/m,
	'soft clipping keeps the complete sequence and every original modification call');
like($clip_view,
	qr/^clip\t2048\t.*\t12S11M\t.*\tMM:Z:C\+m,0,2,2;\tML:B:C,101,202,250\tMN:i:23$/m,
	'full-sequence supplementary alignment receives the same unmodified read-relative tags');
like($clip_view,
	qr/^reverseclip\t16\t.*\t23M\t.*\tMM:Z:C\+m,0;\tML:B:C,77\tMN:i:23$/m,
	'reverse-strand acceptor also receives unchanged original-orientation tags');
unlike($clip_view, qr/\tSA:Z:/, 'stale SA links are removed after alignment-set filtering');
like(slurp($clip_transfer_stats), qr/^output_alignments\t3$/m,
	'tag transfer publishes machine-readable post-policy record counts');
my $hardclip_sam = File::Spec->catfile($tmp, 'hardclip-acceptor.sam');
my $hardclip_bam = File::Spec->catfile($tmp, 'hardclip-acceptor.bam');
write_file($hardclip_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t30\t60\t7H12M4H\t*\t0\t0\t$trimmed_sequence\t$trimmed_quality\tNM:i:0\n");
run_ok('hard-clipped acceptor fixture is created',
	$samtools, 'view', '-b', '-o', $hardclip_bam, $hardclip_sam);
my ($hardclip_status, undef, $hardclip_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $hardclip_bam, '--output', File::Spec->catfile($tmp, 'hardclip-rejected.bam'),
);
ok($hardclip_status != 0, 'hard-clipped acceptor fails instead of receiving guessed modification offsets');
like($hardclip_error, qr/hard-clipped/, 'hard-clipping failure identifies the unsafe condition');
my $edited_sequence = $clip_sequence;
substr($edited_sequence, -1, 1) = substr($edited_sequence, -1, 1) eq 'A' ? 'C' : 'A';
my $edited_sam = File::Spec->catfile($tmp, 'edited-acceptor.sam');
my $edited_bam = File::Spec->catfile($tmp, 'edited-acceptor.bam');
write_file($edited_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t1\t60\t23M\t*\t0\t0\t$edited_sequence\t$clip_quality\tNM:i:1\n");
run_ok('edited acceptor fixture is created', $samtools, 'view', '-b', '-o', $edited_bam, $edited_sam);
my ($edited_status, undef, $edited_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $clip_donor_bam,
	'--acceptor', $edited_bam, '--output', File::Spec->catfile($tmp, 'edited-rejected.bam'),
);
ok($edited_status != 0, 'edited acceptor cannot inherit modification offsets');
like($edited_error, qr/differs from the complete original/, 'sequence mismatch is reported explicitly');
my $invalid_donor_sam = File::Spec->catfile($tmp, 'invalid-donor.sam');
my $invalid_donor_bam = File::Spec->catfile($tmp, 'invalid-donor.bam');
write_file($invalid_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t4\t*\t0\t0\t*\t*\t0\t0\t$clip_sequence\t$clip_quality"
	. "\tMM:Z:C+m,999;\tML:B:C,101\tMN:i:23\n");
run_ok('out-of-range MM donor fixture is created',
	$samtools, 'view', '-b', '-o', $invalid_donor_bam, $invalid_donor_sam);
my ($invalid_status, undef, $invalid_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools, '--donor', $invalid_donor_bam,
	'--acceptor', $clip_acceptor_bam, '--output', File::Spec->catfile($tmp, 'invalid-rejected.bam'),
);
ok($invalid_status != 0, 'out-of-range MM coordinates cannot be published');
like($invalid_error, qr/MM call beyond/, 'out-of-range MM failure identifies the invalid donor');

# samtools -n uses natural numeric ordering and treats zero-padded numeric names
# as one equivalence class.  Primary/supplementary flag ordering can therefore
# interleave distinct exact QNAMEs inside that class.
my $numeric_sequence = 'ACCCGTTACGCC';
my $numeric_quality = 'I' x length($numeric_sequence);
my $numeric_donor_sam = File::Spec->catfile($tmp, 'numeric-donor.sam');
my $numeric_donor_unsorted = File::Spec->catfile($tmp, 'numeric-donor.unsorted.bam');
my $numeric_donor_bam = File::Spec->catfile($tmp, 'numeric-donor.name.bam');
write_file($numeric_donor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:ref\tLN:100\n"
	. join('', map {
		"$_\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
		. "\tMM:Z:C+m,0;\tML:B:C,120\tMN:i:12\n"
	} qw(read10 read002 read2 read02)));
run_ok('natural-QNAME donor fixture is created', $samtools, 'view', '-b',
	'-o', $numeric_donor_unsorted, $numeric_donor_sam);
run_ok('natural-QNAME donor fixture is name sorted', $samtools, 'sort', '-n',
	'-o', $numeric_donor_bam, $numeric_donor_unsorted);
my $numeric_acceptor_sam = File::Spec->catfile($tmp, 'numeric-acceptor.sam');
my $numeric_acceptor_unsorted = File::Spec->catfile($tmp, 'numeric-acceptor.unsorted.bam');
my $numeric_acceptor_bam = File::Spec->catfile($tmp, 'numeric-acceptor.name.bam');
write_file($numeric_acceptor_sam,
	"\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:ref\tLN:100\n"
	. "read2\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\tSA:Z:ref,30,+,12M,60,0;\n"
	. "read2\t2048\tref\t30\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read02\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read002\t2048\tref\t30\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n"
	. "read10\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n");
run_ok('natural-QNAME acceptor fixture is created', $samtools, 'view', '-b',
	'-o', $numeric_acceptor_unsorted, $numeric_acceptor_sam);
run_ok('natural-QNAME acceptor fixture is name sorted', $samtools, 'sort', '-n',
	'-o', $numeric_acceptor_bam, $numeric_acceptor_unsorted);
my $numeric_output = File::Spec->catfile($tmp, 'numeric-transferred.bam');
my $numeric_stats = File::Spec->catfile($tmp, 'numeric-transfer.tsv');
run_ok('natural-equivalence QNAME classes transfer without lexical merge errors',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $numeric_donor_bam, '--acceptor', $numeric_acceptor_bam,
	'--output', $numeric_output, '--stats', $numeric_stats);
my (undef, $numeric_view, undef) = capture($samtools, 'view', $numeric_output);
is_deeply([sort map { (split /\t/)[0] } grep { $_ ne '' } split /\n/, $numeric_view],
	[qw(read02 read10 read2)],
	'exact QNAMEs are joined correctly inside a zero-padded natural-name class');
unlike($numeric_view, qr/\tSA:Z:/, 'primary output cannot retain an SA link to a dropped supplementary');
like(slurp($numeric_stats), qr/^dropped_no_primary_groups\t1$/m,
	'an orphan supplementary-only exact QNAME group is counted and dropped');
like(slurp($numeric_stats), qr/^dropped_supplementary_alignments\t1$/m,
	'the safe default drops supplementary records from otherwise valid groups');

my $variant_donor_sam = File::Spec->catfile($tmp, 'variant-donor.sam');
my $variant_donor_bam = File::Spec->catfile($tmp, 'variant-donor.bam');
my $variant_acceptor_sam = File::Spec->catfile($tmp, 'variant-acceptor.sam');
my $variant_acceptor_bam = File::Spec->catfile($tmp, 'variant-acceptor.bam');
write_file($variant_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+mh?,0;A+123.,0;G+Z,0;\tML:B:C,100,100,200,50\tMN:i:12\n");
write_file($variant_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t0\tref\t1\t60\t12M\t*\t0\t0\t$numeric_sequence\t$numeric_quality\n");
run_ok('multi-code and numeric-ChEBI donor is created', $samtools, 'view', '-b',
	'-o', $variant_donor_bam, $variant_donor_sam);
run_ok('multi-code acceptor is created', $samtools, 'view', '-b',
	'-o', $variant_acceptor_bam, $variant_acceptor_sam);
my $variant_output = File::Spec->catfile($tmp, 'variant-transferred.bam');
run_ok('formal MM alphabetic, numeric, question, and dot forms transfer intact',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $variant_donor_bam, '--acceptor', $variant_acceptor_bam,
	'--output', $variant_output);
my (undef, $variant_view, undef) = capture($samtools, 'view', $variant_output);
like($variant_view, qr/MM:Z:C\+mh\?,0;A\+123\.,0;G\+Z,0;/,
	'valid modern MM group syntax is preserved byte-for-byte');

my $probability_donor_sam = File::Spec->catfile($tmp, 'probability-donor.sam');
my $probability_donor_bam = File::Spec->catfile($tmp, 'probability-donor.bam');
write_file($probability_donor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+mh,0;\tML:B:C,200,200\tMN:i:12\n");
run_ok('invalid multi-code probability donor is created', $samtools, 'view', '-b',
	'-o', $probability_donor_bam, $probability_donor_sam);
my ($probability_status, undef, $probability_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $probability_donor_bam, '--acceptor', $variant_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'probability-rejected.bam'));
ok($probability_status != 0, 'mutually exclusive multi-code probabilities above one fail closed');
like($probability_error, qr/probabilities above 1/,
	'probability-sum failure identifies the inconsistent original read position');

my $missing_mn_sam = File::Spec->catfile($tmp, 'missing-mn.sam');
my $missing_mn_bam = File::Spec->catfile($tmp, 'missing-mn.bam');
write_file($missing_mn_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "variants\t4\t*\t0\t0\t*\t*\t0\t0\t$numeric_sequence\t$numeric_quality"
	. "\tMM:Z:C+m,0;\tML:B:C,100\n");
run_ok('legacy donor without MN is created', $samtools, 'view', '-b',
	'-o', $missing_mn_bam, $missing_mn_sam);
my ($missing_mn_status, undef, $missing_mn_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $missing_mn_bam, '--acceptor', $variant_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'missing-mn-rejected.bam'));
ok($missing_mn_status != 0, 'MN is required by the fail-closed default');
like($missing_mn_error, qr/lacks MN/, 'missing-MN failure points to the legacy override');
my $legacy_output = File::Spec->catfile($tmp, 'missing-mn-allowed.bam');
my $legacy_stats = File::Spec->catfile($tmp, 'missing-mn-allowed.tsv');
run_ok('explicit legacy override accepts a reviewed donor without MN',
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $missing_mn_bam, '--acceptor', $variant_acceptor_bam,
	'--output', $legacy_output, '--stats', $legacy_stats, '--allow-missing-mn');
my (undef, $legacy_view, undef) = capture($samtools, 'view', $legacy_output);
like($legacy_view, qr/\tMN:i:12(?:\t|$)/,
	'legacy override regenerates a current MN on the exact-sequence acceptor');
like(slurp($legacy_stats), qr/^legacy_missing_mn_reads\t1$/m,
	'legacy use is explicit in machine-readable transfer statistics');

my $overlap_acceptor_sam = File::Spec->catfile($tmp, 'overlap-acceptor.sam');
my $overlap_acceptor_bam = File::Spec->catfile($tmp, 'overlap-acceptor.bam');
write_file($overlap_acceptor_sam,
	"\@HD\tVN:1.6\tSO:queryname\n\@SQ\tSN:ref\tLN:100\n"
	. "clip\t0\tref\t1\t60\t5S13M5S\t*\t0\t0\t$clip_sequence\t$clip_quality\n"
	. "clip\t2048\tref\t30\t60\t7S12M4S\t*\t0\t0\t$clip_sequence\t$clip_quality\n");
run_ok('overlapping split-alignment fixture is created', $samtools, 'view', '-b',
	'-o', $overlap_acceptor_bam, $overlap_acceptor_sam);
my ($overlap_status, undef, $overlap_error) = capture(
	$^X, $transfer_script, '--samtools', $samtools,
	'--donor', $clip_donor_bam, '--acceptor', $overlap_acceptor_bam,
	'--output', File::Spec->catfile($tmp, 'overlap-rejected.bam'),
	'--supplementary', 'keep');
ok($overlap_status != 0, 'supplementary retention rejects overlapping native query spans');
like($overlap_error, qr/overlapping primary\/supplementary query spans/,
	'overlap failure explains the double-counting risk');

my ($split_left, $split_right) = ('', '');
for (1 .. 1800) {
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$split_left .= $bases[($seed >> 16) % 4];
	$seed = (1103515245 * $seed + 12345) % 2147483648;
	$split_right .= $bases[($seed >> 16) % 4];
}
my $split_read = substr($split_left, 300, 700) . substr($split_right, 500, 700);
my $split_reference = File::Spec->catfile($tmp, 'split-reference.fa');
my $split_fastq = File::Spec->catfile($tmp, 'split-read.fastq');
write_file($split_reference, ">left\n$split_left\n>right\n$split_right\n");
write_file($split_fastq,
	"\@split\n$split_read\n+\n" . ('I' x length($split_read)) . "\n");
my ($split_status, $split_sam, $split_error) = capture(
	$minimap2, '-a', '-Y', '--secondary=no', '-x', 'map-ont',
	$split_reference, $split_fastq,
);
is($split_status, 0, 'minimap2 creates the split-read clipping fixture')
	or diag($split_error);
my @split_records = map { [split /\t/, $_, -1] }
	grep { $_ ne '' && $_ !~ /^\@/ } split /\n/, $split_sam;
my @supplementary = grep { $_->[1] & 0x800 } @split_records;
ok(@supplementary, 'the chimeric read produces a supplementary alignment');
ok(!(grep { $_->[5] =~ /H/ } @supplementary),
	'-Y prevents hard clipping on every supplementary alignment');
ok(!(grep { length($_->[9]) != length($split_read) } @supplementary),
	'every supplementary alignment retains the complete read sequence for exact transfer');

done_testing();
