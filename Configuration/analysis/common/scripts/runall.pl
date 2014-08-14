#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %arg;
getopts('d:c:r:sj:m:f:R:x:v:', \%arg);

unless ($arg{R} || $arg{r} || !-e $arg{d}) {
    print "********\n$arg{d} exists! - please remove it before running runall!\n" .
        " as an alternative, you can use the -r and/or -R parameter to write into $arg{d}\n";
    exit;
}
unless ($arg{d} && $arg{c} && -f $arg{c}) {
    print <<'USAGE';
Syntax:
 $ runall2.pl -d directoryName -c configFile.py [-r regexp] [-s]
        [-j NJobs] [-x factor] [-m maxEventsPerJob] [-v verbosityLevel] [-f files.txt]

Run runall2.pl to run over all data samples given in files.txt using the
configuration file configFile.py.

-r regexp
If a regexp is given only jobs matching regexp will be submitted (case sensitive).

-R regexp
If a regexp is given only jobs NOT matching the regexp will be submitted (case sensitive).

-s
really submit jobs using nafJobSplitter
otherwise dry run only

-j NJobs
submit NJobs for each dataset. Overwrites default value from files.txt

-x number
multiply NJobs with number (use a number < 1 to submit less jobs, > 1 for more jobs)

-f files.txt (default value: files.txt)
files.txt must contain lines of
cff_file root_output number_of_jobs
# comments with # and empty lines are allowed
Variables can be added via ${VAR} and will be replaced with environment variable $VAR

-v verbosityLevel (default value: 0)
Verbose mode, 0 means no verbosity


USAGE
    exit;
}

my $source = $arg{f} ? $arg{f} : 'files.txt';
open my $IN, '<', $source or die "Cannot open input file: $!\n";
my $ds = Dataset->new($arg{d}, $arg{c}, $arg{s});
my $totalJobs = 0;
while(<$IN>) {
    chomp;
    next if /^\s*#/; #skip comments
    next if exists $arg{'R'} && /$arg{'R'}/;
    next unless /\w/; #skip empty lines
    if (!exists $arg{'r'} || /$arg{'r'}/) {
        my @params = split ' ';
        s!\${(\w+)}!$ENV{$1}!g for @params;
        $params[0] = $arg{'j'} if $arg{'j'};
        $params[0] = multiplyWithJobNumberFactor($params[0]);
        $totalJobs += $params[0];
        $ds->setMaxEventsPerJob($arg{'m'}) if $arg{'m'};
        $ds->setVerbosity($arg{'v'}) if $arg{'v'};
        $ds->prepare(@params);
    }
}
print "Total jobs: $totalJobs\n";
$ds->submit();

sub multiplyWithJobNumberFactor {
    my $no = shift;
    $no *= $arg{'x'} if $arg{'x'};
    $no = int($no);
    return $no < 1 ? 1 : $no;
}

package Dataset;

use strict;
use warnings;
use File::Path;
use File::Copy;
use Cwd qw(abs_path getcwd);

sub new {
    my ($class, $path, $configFile, $reallySubmit) = @_;
    my $self = {};
    $self->{path} = abs_path($path || '.');
    $self->{submit} = $reallySubmit;
    $self->{maxEventsPerJob} = -1;
    $self->{verbosityLevel} = 0;
    $self->{configFile} = $configFile or die "No config file given";
    $self->{configFile} =~ s/-|\./_/g;
    $self->{commands} = [];
    -d $self->{path} or mkpath($self->{path}) or die "Invalid path: >>$self->{path}<<\n$!";
    open my $in, '<', $configFile or die "Cannot read config file $!";
    $self->{fullConfigFile} = do { local $/; <$in> };
    #copy($configFile, "$path/$self->{configFile}") or die "Cannot copy config file to $path: $!";
    bless $self, $class;
}

sub setMaxEventsPerJob {
    my $self = shift;
    $self->{maxEventsPerJob} = shift;
}

sub setVerbosity {
    my $self = shift;
    $self->{verbosityLevel} = shift;
}

sub prepare {
    my ($self, $NJobs, $inputSample, $outputFile, $samplename) = @_;
    $samplename ||= 'standard';
    $inputSample =~ s/\.py$//;
    (my $outputFileWithoutRoot = $outputFile) =~ s/\.root$//;
    my $cmsRunCmdLine = "outputFile=$outputFile,".
        "inputScript=TopAnalysis.Configuration.$inputSample,".
        "samplename=$samplename";
    $self->createConfig($inputSample, $outputFile, $cmsRunCmdLine);
    push @{$self->{commands}}, [$self->{path},
        "-m $self->{maxEventsPerJob} ".
        "-v $self->{verbosityLevel} ".
        "-d $outputFileWithoutRoot ".
        "-M ".($cmsRunCmdLine =~ /PDFWeight/i ? 8000 : 3700)." ".
        "-c '$cmsRunCmdLine' ".
        "$NJobs $self->{configFile}_for_$outputFileWithoutRoot.py"];
}

sub submit {
    my $self = shift;
    my $oldDir = getcwd();
    for (@{$self->{commands}}) {
        my ($path, $cmdline) = @$_;
        chdir($path) or die "Cant chdir to $path: $!";
        if ($self->{submit}) {
            system("nafJobSplitter.pl -W 0 $cmdline"); #nafJobSplitter is in the search path
        } else {
            print "DRY RUN: nafJobSplitter.pl -W 0 $cmdline\n";
        }
    }
    chdir($oldDir) or die "cant chdir back to $oldDir: $!\n";

}

sub createConfig {
    my ($self, $inputSample, $outputFile, $cmsRunCmdLine) = @_;
    my $configFileWithoutPy;
    for ($configFileWithoutPy = "$self->{configFile}_for_$outputFile") {
        s/\.root//g;
        s/\.py//g;
        s/-|\./_/g;
    }

    my $newPY = "$self->{path}/$configFileWithoutPy.py";
    open my $MODIFIED, '>', $newPY or die "$newPY: $!";
    print $MODIFIED "#!cmsRun $configFileWithoutPy.py $cmsRunCmdLine\n";
    print $MODIFIED $self->{fullConfigFile};
    close $MODIFIED;
}

1;

