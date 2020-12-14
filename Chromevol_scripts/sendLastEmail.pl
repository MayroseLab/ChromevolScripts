use strict;
use warnings;

use Getopt::Long;

my $toEmail;
my $JobStatus;
my $id;
my $jobTitle;

GetOptions(
	"jobTitle=s"        => \$jobTitle,
	"JobStatus=s"        => \$JobStatus,
	"toEmail=s"        => \$toEmail,
	"id=s"     => \$id
);

my $from = 'evolseq@tauex.tau.ac.il';

my $subject;
my $subject_fail = "Your ChromEvol job $id has failed";
my $final_status_file = "/bioseq/data/results/ChromEvol/$id/Status.txt";
my $final_pass_file = "/bioseq/data/results/ChromEvol/$id/JOB_PASS.txt";

my $resultsLink = "http://ChromEvol.tau.ac.il/results.html?jobId=$id"."&jobTitle=".$jobTitle;;
my $message .= "Hello,\n\n";
#if (-e $final_pass_file)
if ($JobStatus eq "PASS")
{
	$subject = "Your ChromEvol results of run $id are ready";
	$message .= "$subject at:\n";
	$message .= "Your job was completed successfully\n";
} else {
	$subject = "Your ChromEvol run $id has failed";
	$message .= "$subject_fail at:\n";
	$message .= "Please contact us for more info\n";
}





#my $message .= "Hello,\n\n".$subject. " at:\n";
$message .= "$resultsLink\n\n";


$message .= "Please note: the results will be kept on the server for one month.\n\n";
$message .= "Thank you for using ChromEvol,\n";
$message .= "ChromEvol Team\n";

 
#open(MAIL, "|/usr/sbin/sendmail -t");
open(MAIL, "|/usr/sbin/sendmail -t") or die "can't open mail";
 
# Email Header
print MAIL "To: $toEmail\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n\n";

# Email Body
print MAIL $message;
my $return_code = $? >> 8;
print "print MAIL returned: $return_code\n";

close(MAIL);
#print "Email Sent Successfully to: $toEmail\n";

#use Sys::Hostname;
#$host = hostname;
my $host = `hostname`;

print "hostname: $host\n";


