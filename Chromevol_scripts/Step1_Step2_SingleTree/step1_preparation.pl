#!/usr/bin/perl
#-------------------------------------------------------------------------#
=begin

Program name: 	Ploidy Inference Power Pipeline (for use in queue)
Version:		1.0
Last change:	Feb 2015

Author: 		Lior Glick
E-mail: 		liorglic@mail.tau.ac.il

Description:	This is the source code for the ploidy inference power pipeline.
				Similar to the normal pipeline, this script is used to infer
				ploidy from phylogenies and chromosome numbers. However,
				the power pipeline performs a model choice step on each of the
				given phylogenies, instead of on the guide-tree alone.
				This version is for use on the SGE queue. The pipeline consist
				of three steps:
				1. preparation - step1_preparation.pl -
					reads the pipeline control (parameters) file
					and prepares the files and directories for inference. Then
					sends inferences as jobs to queue, so that each phylogeny is
					analyzed in a single job.
				2. running inferences - run_inferences.pl -
					runs chromEvol several times, sequentially, for each phylogeny.
					First, all required models are tested, then, using the best
					fitting model, simulations are created. One good simulation
					(if such exists) is chosen, and an inference is run on simulated
					data using best fitting model.
				3. summary - step2_summarize.pl - 
					summarizes results: selected models, thresholds, reliability, ploidy.
					Terminates if less than 50 phylogenies have good simulations.
				
				* This script is not for distribution *
				
Usage:			For usage instructions and examples, see the manual:
				http://www.tau.ac.il/~itaymay/cp/chromEvol/chromEvol_v2.0_manual.pdf
				
Requirements:	The script is intended for Unix\Linux machines. Make sure
				you have a working chromEvol executable, as well as the
				chromevol.pm package and the Perl packages detailed in the manual.

=cut
#-------------------------------------------------------------------------#

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use chromevol;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use Cwd 'abs_path';


my $usage = "
Usage: Ploidy Inference Pipeline
	Use control file to specify parameters
\n";
die $usage unless (scalar(@ARGV) == 1);



# Global parameters
my ($controlFile) = @ARGV;
my @paramsList = @{chromevol::read_input($controlFile)};
my ($countsFile,$inTreesFile, $infTreeNum, $simNum,
	 $paramTemplateDir, $runModels, $excludeModels, $baseNum, $ploidyCallType,
	 $duplModelsOnly, $workdir, $chromEvolExe, $plot, $cpus, $logFile, $rootFix, $analysisName) = @paramsList;
	 
if (!-d $workdir){
	system("mkdir -p $workdir") == 0
		or die "ERROR: Could not create output directory $workdir";
}
chromevol::print_params_to_log(\@paramsList,$logFile);
open (LOG, ">> $logFile") or die "ERROR: Can't write to log file '$logFile'\n";

### Run proccess
main();

sub main{
	### 1. Edit counts file
	print LOG "Editing chromosome counts file...\n";
	my $countsEdit = "$workdir/".basename($countsFile)."_edit";
	chromevol::edit_counts_file($inTreesFile, $countsFile, $countsEdit) or die "ERROR: Failed to edit chromosome counts file\n$!";
	if (chromevol::check_counts($countsEdit)){
		print LOG "Error in edited counts file!\n";
		die "Error in edited counts file!   $inTreesFile $countsFile $countsEdit \n";
	}

	### 2. Prepare for ChromEvol inference from real trees
	print LOG "Preparing for chromEvol inference from $infTreeNum phylogenies...\n";
	## choose trees
	my $treesListRef = chromevol::choose_trees($inTreesFile,$infTreeNum);
	if (not $treesListRef){	# not enough trees in input file
		print LOG "Number of required trees is larger than number of trees in input file\n";
		die;
	}
	my $chosenTreesFile = "$workdir/chosen_trees";
	open (my $ofh, "> $chosenTreesFile") or die "ERROR: Can't write chosen trees to file\n";
	print $ofh join("\n",@{$treesListRef});
	close $ofh;
	## prepare inference
	my $realInferDir = "$workdir/infer";
	system("mkdir -p $realInferDir") == 0
		or die "ERROR: Failed to create inference directory $realInferDir";	
	chromevol::prepare_infer_for_power($inTreesFile,$countsEdit, $paramTemplateDir,
		$treesListRef, $rootFix, $realInferDir, $runModels, $excludeModels, $baseNum)
		or die "ERROR: Failed to prepare inference from trees in file $inTreesFile";
		
	### 3. Send inferences to queue
	my $script = $FindBin::Bin."/run_inference_moreSim.pl";
	for (my $n =1; $n <= $infTreeNum; $n++){
		#my $dir = "$realInferDir/infer_tree_special_$n";
		my $dir = "$realInferDir/infer_tree_$n";
		my $command = "perl $script $dir $chromEvolExe $countsEdit $paramTemplateDir";
		my $jobName = "chromEvol_$analysisName\_tree_$n";
		send_to_queue($command, $dir, $jobName);
	}
	# create flag file to signal end of step
	my $flagFile = "$workdir/preparation_step_complete";
	open(F,"> $flagFile");
	close F;		
}

sub send_to_queue{
	
	my ($command, $outDir, $jobName) = @_;
	my $shFile = "$outDir/$jobName.sh";
	my $stderrFile = "$outDir/$jobName.ER";
	my $stdoutFile = "$outDir/$jobName.OU";
	
	print($shFile);
	open (SH_SCRIPT,"> $shFile");
	print SH_SCRIPT '#!/bin/tcsh',"\n";
	print SH_SCRIPT '#$ -N ',"$jobName\n";
	print SH_SCRIPT '#$ -S /bin/tcsh',"\n";
	print SH_SCRIPT '#$ -cwd',"\n";
	print SH_SCRIPT '#$ -e ',$stderrFile,"\n";
	print SH_SCRIPT '#$ -o ',$stdoutFile,"\n";
	print SH_SCRIPT "module load  perl/perl-5.16.3\n";
	print SH_SCRIPT "$command\n";
	close(SH_SCRIPT);
			
	system ("qsub -l itaym -p -2 $shFile") == 0
			or die "ERROR: Can't send $shFile to queue\n";
			
}