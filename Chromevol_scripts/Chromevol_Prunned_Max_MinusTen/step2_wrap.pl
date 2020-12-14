use strict;

# get working_dir and list of genera
my ($working_dir,$list,$step2_script,$output_file) = @ARGV;

chdir($working_dir) or die "Can't open directory $working_dir\n";
my $ref_genera = get_genera($working_dir,$list);
my @genera = @{$ref_genera};
chomp @genera;

# try running step2 and get the process status
my $control_file;
my $qstatString = `qstat -r | grep "Full jobname:" | awk '{print $3}'`;	# get full job names, including genus names
open (OUT,">$output_file") or die "Can't open $output_file\n";
my $status;
foreach my $genus (@genera){
	print ("$genus\n");
	if ($qstatString =~ m/_$genus\_/){
		$status = "Still running...";
	}
	else{
		#/groups/itay_mayrose/michaldrori/MSA_JmodelTree/output/Acanthus/Acanthus_Chromevol
		$control_file = $working_dir.$genus."/".$genus."_Chromevol/PIP_control";
		$status = check_status($genus,$working_dir);
		if ($status eq "ready"){
			print("perl $step2_script $control_file\n");
			$status = `perl $step2_script $control_file`;
			chomp $status;
			if ($status eq ""){
				$status = "done";
			}
		}
	}
	print OUT ($genus.",".$status."\n");
}
close (OUT);

sub get_genera{
	my ($working_dir,$list) = @_;
	my @genera;
	if ($list eq ""){ # the list is empty - get all sub-directories from the working dir
		opendir my($dh), $working_dir or die "Can't open dir $working_dir\n";
		@genera = readdir $dh;
		closedir $dh;
	} else { # the list is not empty - read it and go over the genera in it
		open (GENERA, "<$list") or die "Can't open $list\n";
		@genera = <GENERA>;
		close (GENERA)
	}
	return (\@genera);
}

sub check_status{
	my ($genus,$working_dir) = @_;
	my $status;
	my $curr_dir = $working_dir.$genus."/".$genus."_Chromevol/chromevol_out/";
	chdir ($curr_dir);
	if (-e "ploidy.csv"){
		$status = "done";
	} else {
		$status = "ready";
	}
	chdir($working_dir);
	return($status);
}
