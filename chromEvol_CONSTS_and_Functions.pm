package chromEvol_CONSTS_and_Functions;

use strict;
use warnings;
use File::Slurp;
use File::Path;# qw(make_path);
use File::Basename;


use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(WriteToFile WriteToFileWithTimeStamp ReadFromFile CheckInputTypeValid ReadParam ReadToArrayParam CheckIfFileExists ReadParamLineToArray);



######### Consts ###########
use constant CHROMEVOL_SCRIPTS_PATH		=> "/bioseq/chromEvol/Chromevol_scripts";
use constant CHROMEVOL_PARAM_PATH		=> "/bioseq/chromEvol/power_PARAM_templates/";
use constant CHROMEVOL_EXE_PATH			=> "/bioseq/chromEvol";
use constant SERVER_NAME 				=> "chromevol";
use constant RESULTS_DIR_ABSOLUTE_PATH 	=> "/bioseq/data/results/chromEvol";
use constant RESULTS_LINK 				=> "chromEvol_results";
use constant LOG_FILE 					=> "PIP_control";		#"output.log";
use constant ERROR_STATUS_LOG_FILE 		=> "Status.txt";		#"error.OU";
use constant OUTPUT_FILES_LIST 			=> "output_files.dat";
use constant LOG_DIR_ABSOLUTE_PATH 		=> "/bioseq/data/results/chromEvol/logs/";
use constant MAIN_PIPELINE_LOG 			=> "chromEvol.log";
#use constant PYTHON_MODULE_TO_LOAD 		=> 'python/python-3.3.0';
#use constant PYTHON_MODULE_TO_LOAD 		=> 'python/anaconda3-5.0.0';
use constant PYTHON_MODULE_TO_LOAD 		=> 'python/python-anaconda3.6.5-michaldrori';
#use constant PERL_MODULE_TO_LOAD 		=> 'perl/perl-5.16.3';
use constant PERL_MODULE_TO_LOAD 		=> 'perl/perl-5.28.1';
use constant QSUB_JOB_NUM_FILE			=> "qsub_job_num.dat";
use constant RESULTS_PAGE_URL			=> "/results.html";
use constant RESULTS_COUNTS_PAGE_URL	=> "/counts_results.html";
use constant FILENAME_PARAMS_TXT		=> "params.txt";
use constant FILENAME_TREE_FILE			=> "tree.tre";
use constant FILENAME_COUNTS_FILE		=> "counts.txt";
use constant FILENAME_PARAMS_HTML		=> "params.html";
use constant FILENAME_LOCUS_TREENUM_TXT	=> "NumOfClusterTrees.txt";
use constant FILENAME_EMAIL				=> "email.txt";
use constant FILENAME_SELECTED_MODELS	=> "selected_models.txt";
use constant FILENAME_STATUS_FLAGS		=> "status_flags.txt";
use constant FILENAME_INPUT_REPORT		=> "taxa_input_report.txt";
use constant FILENAME_MA_PLOIDY_FLAGS	=> "ma_ploidy_flags.txt";
use constant FILENAME_DONE_FILE			=> "DONE.txt";
use constant FILENAME_CHROM_COUNTS		=> "chrom_counts_data.txt";
use constant FILENAME_ADEQUACY			=> "Adeq_Results.csv";
use constant FILENAME_BEST_MODELS		=> "BestModel_per_Tree.txt";
use constant FILENAME_PIP	=> "PIP_control";


######### Functions ############
sub WriteToFile
{
	my ($file, $message, $shouldOverwrite) = @_;
	
	# creating file dir if doesnt exist
	my ($fileName, $fileDir) = fileparse($file);
	#make_path($fileDir);
	mkpath($fileDir);
	
	#$message =~ s/^\s+//;
	
	$message =~ s/\R/\n/g;

	if (defined $shouldOverwrite && $shouldOverwrite){
		write_file($file, "$message\n");
	}
	else {
		append_file($file, "$message\n");
	}	
}


sub ReadParam
{
	my ($file, $paramName) = @_;
	my $NameValue = "NONE";
	#Check if file exists:
	if (-e $file)
	{
		open my $fh, '<', $file or die "Could not open '$file' $!\n";
		while (my $line = <$fh>) {
			chomp $line;
			if (index($line, $paramName) != -1)
			{
				my @fields = split /:/, $line;
				$NameValue = $fields[1];
				return $NameValue;	
			}
		}
	}
	$NameValue  = "NoData";
	return $NameValue;
}

sub ReadToArrayParam  
{
	my ($file, $paramName) = @_;
	my @values = ();
	
	#Check if file exists:
	if (-e $file)
	{
		open my $fh, '<', $file or die "Could not open '$file' $!\n";
		while (my $line = <$fh>) {
			chomp $line;
			#@values = split(',', $line);
			@values = split /,/, $line;
			return \@values;
		}
	}
	return \@values;
}

sub ReadParamLineToArray  
{
	my ($str_param) = @_;
	my @values = ();
	
	if ($str_param eq ""){
		chomp $str_param;
		@values = split /,/, $str_param;
		return \@values;
	} else {
		return \@values;
	}

}

sub CheckIfFileExists  
{
	my ($file) = @_;
	my $flag_exists = 'NO';
	
	#Check if file exists:
	if (-e $file)
	{
		$flag_exists = 'YES';
		return $flag_exists;
	}
	return $flag_exists;
}

sub ReadFromFile
{
	my ($file, $defaultValue) = @_;
	
	if (defined $defaultValue)
	{
		if (-e $file)
		{
			my $line = read_file($file);
			return $line;
		}
		else
		{
			# this is ugly. delete this after renaming all relevant files to dat.
			my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	if ($ext eq ".dat")	{
		 		my $newFile = substr($file, 0, -4);
		 		
		 		if (-e $newFile) {
					my $line = read_file($newFile);
					return $line;
				}
		 	}
			
			return $defaultValue;
		}
	}
	else
	{
		my @lines;
		
		if (-e $file)
		{
			@lines = read_file($file);
		}
		else {
			return 'abc975'
			# this is ugly. delete this after renaming all relevant files to dat.
			#my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	#if ($ext eq ".dat")	{
		 	#	my $newFile = substr($file, 0, -4);
		 	#	
		 	#	if (-e $newFile) {
			#		@lines = read_file($newFile);
			#	}
		 	#}
		}
		
		return @lines;
	}
}