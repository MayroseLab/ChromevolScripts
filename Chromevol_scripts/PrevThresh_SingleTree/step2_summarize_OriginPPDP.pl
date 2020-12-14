use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use chromevol;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
*STDERR = *STDOUT;

my ($controlFile) = @ARGV;
my @paramsList = @{chromevol::read_input($controlFile)};
my ($countsFile,$inTreesFile, $infTreeNum, $simNum,
	 $paramTemplateDir, $runModels, $excludeModels, $baseNum, $ploidyCallType,
	 $duplModelsOnly, $workdir, $chromEvolExe, $plot, $cpus, $logFile, $rootFix, $analysisName) = @paramsList;

my $inferDir = "$workdir/infer";
## check if inferences step completed
my $completedPhylogenies = chromevol::check_inferences($inferDir,$infTreeNum,$workdir);
if ($completedPhylogenies == 0){
	die "Error while running\n";
}
else{
	my @completedTreesNums = @{$completedPhylogenies};
	if (scalar(@completedTreesNums) < 50){
		## summarize models
		my $modelsSumFile = "$workdir/models_summary";
		my $models_sum = chromevol::summarize_models($inferDir,$modelsSumFile);
		#if ($models_sum == 0){
		#	open(NoDupleMoreThan50,">$workdir/NoDuple_MoreThan50.txt");
		#	print NoDupleMoreThan50 " More than 50% NoDuple models\n";
		#	die "More than 50% NoDuple models\n";
		#}
		#else{
		#	open(NotEnoughTrees,">$workdir/Less_than_50_phylogenies.txt");
		#	print NotEnoughTrees "completedTreesNums: (@completedTreesNums)\n";
		#	die "Less than 50 phylogenies to work with\n";
		#}
	}
}

## summarize models
my $modelsSumFile = "$workdir/models_summary_origin";
chromevol::summarize_models($inferDir,$modelsSumFile);

## Determine threshold for calling diploids/polyploids (using simulations)
my $thresFilePP = "$workdir/thresholds_PP_origin";
my $thresFileDP = "$workdir/thresholds_DP_origin";
my ($thresPP, $thresDP) = @{chromevol::determine_thresholds_for_power_OneTree($inferDir,$thresFilePP, $thresFileDP, $ploidyCallType)}
	or die "ERROR: Failed to determine thresholds\n";

### Compute reliability from real inference
my $relInfer = "$workdir/reliability_infer.txt";
chromevol::compute_reliability_from_real_for_power($inferDir, $thresPP, $thresDP, $relInfer, $ploidyCallType)
	or die "ERROR: Failed to compute reliability from real inferences\n";
my $inferReliabilityRef = chromevol::parse_reliability_file($relInfer);	
#
#
### Compute reliability from simulations inference
my $relSim = "$workdir/reliability_sims.txt_origin";
#chromevol::compute_reliability_from_sim_for_power_OneTree($inferDir,$thresDP, $thresPP, $ploidyCallType,$relSim,$inferReliabilityRef)
#	or die "ERROR: Failed to compute reliability from simulations\n";

## Summarize reliability
# find counts edit file
opendir(DIR,$workdir);
my @countsEdit = grep {/counts_max_edit$/} readdir(DIR);	#MD changed from: my @countsEdit = grep {/counts_edit$/} readdir(DIR);
my $countsEditFile = "$workdir/$countsEdit[0]";
my $finalOut = "$workdir/ploidy.csv";
chromevol::summarize_reliability_for_power($relSim,$relInfer,$countsEditFile,$finalOut)
	or die "ERROR: Failed to summarize reliability\n";