use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Meth2Rep qw(
	read_target_mgs read_modbam_manifest read_mgs_report
	build_meth2rep_plan representative_fasta_path
);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub write_gzip {
	my ($path, $contents) = @_;
	gzip(\$contents => $path) or die "Cannot write $path: $GzipError";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $selection_file = File::Spec->catfile($tmp, 'mgs.txt');
write_file($selection_file, "# selected targets\nMGS.2\nMGS.1 # retained\n");
my $selected = read_target_mgs(mgs => 'MGS.3,MGS.1', mgs_file => $selection_file);
is_deeply([sort keys %{$selected}], [qw(MGS.1 MGS.2 MGS.3)],
	'comma and file MGS selections form a deterministic union');

my (%map, %groups);
$map{opt}{smpl_order} = [qw(S1 S2)];
$map{altNms} = {};
for my $spec (
	[qw(S1 G1 ONT)],
	[qw(S2 G2 PB)],
) {
	my ($sample, $group, $technology) = @{$spec};
	my $work = File::Spec->catdir($tmp, $sample);
	my $assembly = File::Spec->catdir($tmp, "assembly-$sample");
	make_path(
		File::Spec->catdir($work, 'assemblies', 'metag'),
		File::Spec->catdir($work, 'mapping'),
		File::Spec->catdir($assembly, 'Binning', 'SB'),
	);
	write_file(File::Spec->catfile($work, 'assemblies', 'metag', 'assembly.txt'), "$assembly\n");
	write_file(File::Spec->catfile($assembly, 'scaffolds.fasta.filt'), ">$sample-contig\nACGT\n");
	write_file(File::Spec->catfile($assembly, 'Binning', 'SB', $sample),
		"Sequence ID\tBin\n$sample-contig\t" . ($sample eq 'S1' ? '1.fa.gz' : '2.fa.gz') . "\n");
	my $cram = File::Spec->catfile($work, 'mapping', "$sample-smd.cram");
	write_file($cram, 'fixture');
	write_file("$cram.sto", 'done');
	my $modbam = File::Spec->catfile($tmp, "$sample.mod.bam");
	write_file($modbam, 'fixture');
	$map{$sample} = {
		SmplID => $sample, wrdir => "$work/", AssGroup => $group,
		SeqTech => $technology, SupportReads => '', hasPrimaryRds => 1,
	};
	$groups{$group} = { SmplID => [$sample], wrdir => ["$work/"] };
}

my $manifest_file = File::Spec->catfile($tmp, 'modbams.tsv');
write_file($manifest_file,
	"sample\tscope\ttechnology\tmodbam\n"
	. "S1\tprimary\tONT\tS1.mod.bam\n"
	. "S2\tprimary\tPB\tS2.mod.bam\n"
	. "S9\tprimary\tONT\tarchived-S9.mod.bam\n");
my $manifest = read_modbam_manifest($manifest_file);
is($manifest->{"S1\tprimary"}[0]{technology}, 'ONT',
	'modBAM manifest records explicit technology and scope');
like($manifest->{"S2\tprimary"}[0]{modbam}, qr/\Q$tmp\E\/S2\.mod\.bam\z/,
	'relative donor paths resolve against the manifest directory');

my $duplicate_manifest = File::Spec->catfile($tmp, 'duplicate-modbams.tsv');
write_file($duplicate_manifest,
	"sample\tscope\ttechnology\tmodbam\n"
	. "S1\tprimary\tONT\tS1.mod.bam\n"
	. "S1\tprimary\tONT\tS1.mod.bam\n");
eval { read_modbam_manifest($duplicate_manifest) };
like($@, qr/repeats donor/,
	'an exact repeated sample/scope donor fails instead of being scanned twice');
my $multi_manifest_file = File::Spec->catfile($tmp, 'multi-modbams.tsv');
write_file($multi_manifest_file,
	"sample\tscope\ttechnology\tmodbam\n"
	. "S1\tprimary\tONT\tS1.mod.bam\n"
	. "S1\tprimary\tONT\tS1-part2.mod.bam\n");
my $multi_manifest = read_modbam_manifest($multi_manifest_file);
is(scalar(@{$multi_manifest->{"S1\tprimary"}}), 2,
	'multiple distinct original modBAMs are declared for one sample/scope');

my $report_file = File::Spec->catfile($tmp, 'MAGvsGC.txt.gz');
write_gzip($report_file,
	"MAG\tMGS\tRepresentative4MGS\tCompleteness\n"
	. "S1__1.fa.gz\tMGS.1\t*\t95\n"
	. "S2__2.fa.gz\tMGS.1\t\t90\n"
	. "Cano__42\tMGS.2\t*\t80\n");
my $representatives_dir = File::Spec->catdir($tmp, 'representatives');
make_path($representatives_dir);
my $representative_fasta = representative_fasta_path(
	$representatives_dir, 'MGS.1', 'S1__1.fa.gz',
);
write_gzip($representative_fasta, ">S1-contig\nACGT\n");

my $report = read_mgs_report($report_file, { 'MGS.1' => 1, 'MGS.2' => 1 });
is($report->{'MGS.1'}{representative_mag}, 'S1__1.fa.gz',
	'MGS report representative is selected by the named header and star marker');
is($report->{'MGS.2'}{reason}, 'canopy_only_no_mag_reference',
	'Canopy-only targets are represented as an explicit unavailability state');

my $unrelated_work = File::Spec->catdir($tmp, 'S3-incomplete');
make_path($unrelated_work);
$map{S3} = {
	SmplID => 'S3', wrdir => "$unrelated_work/", AssGroup => 'G3',
	SeqTech => 'ONT', SupportReads => '', hasPrimaryRds => 1,
};
$groups{G3} = { SmplID => ['S3'], wrdir => ["$unrelated_work/"] };

my $plan = build_meth2rep_plan(
	map => \%map, assembly_groups => \%groups,
	manifest => $manifest, manifest_file => $manifest_file,
	selected => { 'MGS.1' => 1, 'MGS.2' => 1 },
	mgs_report => $report_file, representatives_dir => $representatives_dir,
	binner => 'SB', modes => [qw(mgs2rep rep2rep)],
);
ok(!exists($plan->{scopes}{"S3\tprimary"}),
	'an explicitly selected MGS does not inspect or require unrelated assembly groups');
is_deeply($plan->{manifest_unused}, ["S9\tprimary"],
	'an unavailable donor in an unrelated manifest row does not block selected work');
is(scalar(@{$plan->{targets}{"mgs2rep\tMGS.1"}{sources}}), 2,
	'mgs2rep uses every non-Canopy MAG in the selected MGS');
is(scalar(@{$plan->{targets}{"rep2rep\tMGS.1"}{sources}}), 1,
	'rep2rep uses only the representative MAG');
is_deeply([sort keys %{$plan->{targets}{"mgs2rep\tMGS.1"}{scope_keys}}],
	["S1\tprimary", "S2\tprimary"],
	'mgs2rep follows member MAG assembly groups to every eligible sample scope');
is_deeply([sort keys %{$plan->{targets}{"rep2rep\tMGS.1"}{scope_keys}}],
	["S1\tprimary"],
	'rep2rep does not pull candidates from non-representative member assemblies');
ok($plan->{contig_targets}{File::Spec->catdir($tmp, 'assembly-S2')}{'S2-contig'}{"mgs2rep\tMGS.1"},
	'member-MAG contigs are connected to the mgs2rep target');
ok(!exists($plan->{contig_targets}{File::Spec->catdir($tmp, 'assembly-S2')}{'S2-contig'}{"rep2rep\tMGS.1"}),
	'non-representative contigs are excluded from rep2rep');
is(scalar(grep { $_->{mgs} eq 'MGS.2' } @{$plan->{unavailable}}), 2,
	'Canopy-only MGS gets one explicit unavailable result per requested mode');

my %incomplete_manifest = %{$manifest};
delete $incomplete_manifest{"S2\tprimary"};
eval {
	build_meth2rep_plan(
		map => \%map, assembly_groups => \%groups,
		manifest => \%incomplete_manifest, manifest_file => $manifest_file,
		selected => { 'MGS.1' => 1 }, mgs_report => $report_file,
		representatives_dir => $representatives_dir, binner => 'SB',
		modes => ['mgs2rep'],
	);
};
like($@, qr/No original modBAM.*S2.*primary/,
	'a missing donor relationship fails before any alignment is attempted');

my %missing_modbam = map { $_ => [map { +{%{$_}} } @{$manifest->{$_}}] } keys %{$manifest};
$missing_modbam{"S2\tprimary"}[0]{modbam} = File::Spec->catfile($tmp, 'missing-S2.mod.bam');
eval {
	build_meth2rep_plan(
		map => \%map, assembly_groups => \%groups,
		manifest => \%missing_modbam, manifest_file => $manifest_file,
		selected => { 'MGS.1' => 1 }, mgs_report => $report_file,
		representatives_dir => $representatives_dir, binner => 'SB',
		modes => ['mgs2rep'],
	);
};
like($@, qr/Original modBAM is missing or empty.*S2.*primary/,
	'a selected donor path must resolve to a nonempty modBAM');

my $mgs_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
unlike($mgs_source, qr/meth2rep|mgs2rep|rep2rep/i,
	'MGS runs independently and has no meth2rep launch path');
my ($mgs_checkpoint_parameters) = $mgs_source =~ /%checkpointParameters\s*=\s*\((.*?)\n\);/s;
ok(defined($mgs_checkpoint_parameters), 'the established MGS checkpoint parameter block is identifiable');
unlike($mgs_checkpoint_parameters, qr/meth2rep/i,
	'meth2rep options do not invalidate the established MGS clustering checkpoint');
my $gene_cat_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'geneCat.pl'));
unlike($gene_cat_source, qr/meth2rep|mgs2rep|rep2rep/i,
	'geneCat runs independently and does not forward meth2rep options');
my $config_source = slurp(File::Spec->catfile($Bin, '..', 'Mods', 'config_internal.txt'));
unlike($config_source, qr/^meth2rep_scr\t/m,
	'the standalone tool does not modify core program configuration');
unlike($config_source, qr/^modkit\t/m,
	'modkit is not a runtime tool dependency');
my $executor_source = slurp(File::Spec->catfile(
	$Bin, '..', 'secScripts', 'MGS', 'meth2rep.pl',
));
like($executor_source,
	qr/\$options\{minimap2\}.*?'-Y'.*?'--secondary=no'/s,
	'fresh representative mapping retains full SEQ on supplementary records for exact tag transfer');
like($executor_source, qr/transfer_mod_tags\.pl/,
	'the executor uses the first-party tag transfer');
unlike($executor_source, qr/modkit/i,
	'the executor never invokes modkit');
like($executor_source,
	qr/for my \$scope_key.*?scan_scope_candidates.*?for my \$modbam \(sort keys %donor_scopes\).*?filter_bam_by_names/s,
	'all CRAM candidate lists are compiled before one extraction pass per distinct donor modBAM');
like($executor_source,
	qr/my \$index_key\s*=\s*join\("\\0", \$reference, \$preset\).*?'-x', \$preset.*?'-d', \$index/s,
	'minimap2 indexes are keyed and built with the same technology-specific preset used for mapping');
like($executor_source, qr/assert_source_donor_identity/,
	'candidate provenance requires native-sequence agreement between assembly CRAM and donor modBAM');
like($executor_source, qr/'output-format=s'/,
	'the standalone module exposes BAM or self-contained CRAM publication');
like($executor_source, qr/'supplementary-alignments=s'/,
	'supplementary retention is an explicit output policy rather than an implicit side effect');
unlike($executor_source, qr/getProgPaths\(/,
	'standalone tool discovery cannot silently select bamFilter from another checkout');

done_testing();
