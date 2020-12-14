use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use chromevol;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
my ($inferDir, $chromEvolExe, $countsEdit, $paramTemplateDir) = @ARGV;

### get list of parameter files to run (and fetch tree file for later use)
my @paramFiles;
my $treeFile;
opendir (DIR, $inferDir) or die "Can't open directory $inferDir!\n$!";
while (my $file = readdir(DIR)) {
	if ($file =~ m/^param_/){
		push(@paramFiles,"$inferDir/$file");
	}
	elsif ($file =~ m/^tree/){
		$treeFile = "$inferDir/$file";
	}
}
### run inferences
foreach my $paramFile (@paramFiles){
	system("$chromEvolExe $paramFile") == 0
			or die "ERROR: chromEvol run on $paramFile failed.\n";
}

### summarize results
my @chromNums = @{chromevol::calc_min_max_chrom_num($countsEdit)};
my $maxChromNum = $chromNums[2];	# (max in data + 10)*2;
my $resultSumFile = "$inferDir/result_sum";
my $rootFreqFile = "$inferDir/root_freq";
chromevol::parse_chromevol_results($inferDir, $resultSumFile, $maxChromNum,$rootFreqFile);
my $bestModel = chromevol::get_best_model($resultSumFile);
if ($bestModel eq "CONST_RATE_NO_DUPL" or $bestModel eq "LINEAR_RATE_NO_DUPL"){
	my $flagFile = "$inferDir/NO_DUPL";
	open(F,"> $flagFile");
	close F;
	die "NO_DUPL model chosen, no ploidy inference will be performed for this phylogeny\n";
}

### simulate
## prepare simulations 
my $simDir = "$inferDir/simulation/";
mkdir $simDir;
my $simParamTemplate = "$paramTemplateDir/param_SIM";
my $minInData = $chromNums[0];
my $maxInData = $chromNums[1];
my $rangeInData = $maxInData - $minInData;
my $simNum = 100;
my $simRunFilesListRef = chromevol::prepare_simulations($inferDir, $rootFreqFile, $simParamTemplate, $simNum,
$maxChromNum, $rangeInData, $simDir,0)
	or die "Error: Failed to prepare simulations\n";		
## run simulations
my @simRunFiles = @{$simRunFilesListRef};
foreach my $simRunFile (@simRunFiles){
	system("$chromEvolExe $simRunFile") == 0
			or die "ERROR: chromEvol simulation from file $simRunFile failed.\n";
}
## test simulations
my $simulationDir = "$inferDir/simulation";
my $goodSim = chromevol::test_simulations($simNum, $maxInData, $minInData, $simulationDir);
while($goodSim == 0){	# while a good simulation was not found
	my $fakeRootFile = "$inferDir/fake_root_freq";
	chromevol::prepare_fake_root_freq($rootFreqFile, $fakeRootFile);
	$rootFreqFile = $fakeRootFile;
	$simRunFilesListRef = chromevol::prepare_simulations($inferDir, $fakeRootFile, $simParamTemplate, $simNum,
	$maxChromNum, $rangeInData, $simDir,1)
		or die "Error: Failed to prepare fake simulations\n";
	my @simRunFiles = @{$simRunFilesListRef};
	foreach my $simRunFile (@simRunFiles){
		system("$chromEvolExe $simRunFile") == 0
				or die "ERROR: chromEvol simulation from file $simRunFile failed.\n";
	}

	$goodSim = chromevol::test_simulations($simNum, $maxInData, $minInData, $simulationDir);
}

if ($goodSim){	# if a good simulation was found
	## remove all but one good simulation
	for (my $n = 1; $n <= $simNum*2; $n++){
		if ($n != $goodSim){
			#system("rm -rf $simulationDir/dir_sim_$n");
			#system("rm -rf $simulationDir/param_sim_$n");
		}
	}
	
	## prepare infer from sim
	my $goodSimDir = "$simulationDir/dir_sim_$goodSim/0";
	my $simCountsFile = "$goodSimDir/simCounts.txt";
	my $simInferDir = "$simulationDir/sim_infer";
	mkdir $simInferDir;
	my $simDataFile = "$simInferDir/sim_data.counts"; # same as simCounts, but with missing species replaced with x
	chromevol::parse_sim_results($countsEdit,$simCountsFile,$simDataFile);	
	chromevol::prepare_initial_run($treeFile, $simDataFile, $paramTemplateDir, $bestModel, 0, 0, 0, $simInferDir);
	
	## run inference from sim
	my $runFile = "$simInferDir/param_$bestModel\.txt";
	system("$chromEvolExe $runFile") == 0
			or die "ERROR: chromEvol run on file $runFile failed.\n";
}
else{
	warn "A good simulation could not be acquired for $treeFile.\n";
}
# create flag file to signal end of step
my $flagFile = "$inferDir/inference_step_complete";
open(F,"> $flagFile");
close F;




	