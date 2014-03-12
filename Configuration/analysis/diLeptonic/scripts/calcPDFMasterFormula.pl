#!/usr/bin/perl -w

# % cat Plots_temp/PDF_CENTRAL/combined/HypToppTLaTeX.txt
# Variable: HypToppT   Channel: Dilepton Combined
# BinCenter & LowXbinEdge  &  HighXbinEdge  &   DiffXSec  &  StatError(\%)  & SystError(\%)  & TotalError(\%) \\
# \hline
# $40$ & $0$ to $80$   &  0.0051167  &   1.33 &    0 &    1.33 \\
# $105$ & $80$ to $130$   &  0.0058073  &   1.5 &    0 &    1.5 \\
# $165$ & $130$ to $200$   &  0.0030112  &   1.45 &    0 &    1.45 \\
# $250$ & $200$ to $300$   &  0.00084305  &   1.76 &    0 &    1.76 \\
# $350$ & $300$ to $400$   &  0.00016227  &   3.4 &    0 &    3.4 \\

use strict;
use warnings;
use List::Util qw(max);
use File::Path qw(mkpath);

use Data::Dumper;

my $outputPath = "UnfoldingResults/PDF_";

sub pdfSum {
    my ($quantity, $channel) = @_;
    my $nominal = UnfoldedResult->new("UnfoldingResults/PDF_CENTRAL/Nominal/$channel/${quantity}Results.txt");
    #print Dumper $result;
    my (@up, @down);
    for my $var_no (1..22) { 
        push @up, UnfoldedResult->new("UnfoldingResults/PDF_${var_no}_UP/Nominal/$channel/${quantity}Results.txt");
        push @down, UnfoldedResult->new("UnfoldingResults/PDF_${var_no}_DOWN/Nominal/$channel/${quantity}Results.txt");
    }
    
    my $filename = "$outputPath/$channel/${quantity}Results.txt";
    open my $OUTFH, '>', $filename or die $!;
    print "Writing to $filename\n";
    
    for my $bin (0..@{$nominal->{bins}}-1) {
        my $binNom = $nominal->{bins}->[$bin];
        my $nominalXsec = $binNom->{-xsec};
        my ($upSUM2, $downSUM2) = (0,0);
        for my $no (0..21) {
            my $upMnom = $up[$no]{bins}[$bin]{-xsec} - $nominalXsec;
            my $downMnom = $down[$no]{bins}[$bin]{-xsec} - $nominalXsec;
            $upSUM2 += max(0, $upMnom, $downMnom)**2;
            $downSUM2 += max(0, -$upMnom, -$downMnom)**2;
        }
        my $relUncUp = sqrt($upSUM2)/$nominalXsec;
        my $relUncDown = sqrt($downSUM2)/$nominalXsec;
#         my $unc = 0.5 * ($relUncUp + $relUncDown);
        my $unc = max($relUncUp, $relUncDown);
        #my $line = "XAxisbinCenters[bin]: $binNom->{-center} bin: $binNom->{-from} to $binNom->{-to} SystematicRelError: $unc (+$relUncUp -$relUncDown)\n";
        my $line = "XAxisbinCenters[bin]: $binNom->{-center} bin: $binNom->{-from} to $binNom->{-to} SystematicRelError: $unc\n";
        print $line;
        print $OUTFH $line;
    }

}

sub pdfSumIncl {
    my ($channel) = @_;
#    my $nominalIncl = InclusiveResult->new("Plots/PDF_CENTRAL/Nominal/$channel/InclusiveXSec.txt");
    my $nominalIncl = InclusiveResult->new("Plots/Nominal/$channel/InclusiveXSec.txt");
    my (@upIncl, @downIncl);
    for my $var_no (1..22) { 
        push @upIncl, InclusiveResult->new("Plots/PDF_${var_no}_UP/Nominal/$channel/InclusiveXSec.txt");
        push @downIncl, InclusiveResult->new("Plots/PDF_${var_no}_DOWN/Nominal/$channel/InclusiveXSec.txt");
    }
    
    my ($upInclSUM2, $downInclSUM2) = (0,0);
    for my $no (0..21) {
        my $upMnom = $upIncl[$no]->{-xsec} - $nominalIncl->{-xsec}; 
        my $downMnom = $downIncl[$no]->{-xsec} - $nominalIncl->{-xsec};
        $upInclSUM2 += max(0, $upMnom, $downMnom)**2;
        $downInclSUM2 += max(0, -$upMnom, -$downMnom)**2;
    }

    storeInclXsec("PDF_UP", $channel, $nominalIncl->{-xsec} + sqrt($upInclSUM2));
    storeInclXsec("PDF_DOWN", $channel, $nominalIncl->{-xsec} - sqrt($downInclSUM2));    
}

sub storeInclXsec {
    my ($syst, $channel, $xsec) = @_;
    mkpath("Plots/$syst/$channel"); 
    my $result = "Systematic: $syst Channel: /$channel InclXSection: $xsec AbsStatError: 0\n";
    
    print "Inclusive xsec:\n$result";
    
    open my $FH, '>', "Plots/$syst/$channel/InclusiveXSec.txt" or die $!;
    print $FH $result;
}

for my $channel qw(ee emu mumu combined) {
    mkpath("$outputPath/$channel");
    for my $plot qw(
    HypLeptonpT
    HypLeptonpTLead
    HypLeptonpTNLead
    HypLeptonEta
    HypLeptonEtaLead
    HypLeptonEtaNLead
    HypLeptonBjetMass
    HypLLBarpT
    HypLLBarMass
    HypBJetpT
    HypBJetpTLead
    HypBJetpTNLead
    HypBJetEta
    HypBJetpTNLead
    HypBJetEta
    HypBJetEtaLead
    HypBJetEtaNLead
    HypTopRapidity
    HypTopRapidityLead  
    HypToppT
    HypToppTLead
    HypToppTNLead
    HypTTBarRapidity
    HypTTBarpT
    HypTTBarMass
    HypTTBarDeltaPhi
    HypBBBarpT
    HypBBBarMass
    HypToppTTTRestFrame
    HypTTBarDeltaRapidity
    HypJetMultpt30
    HypJetMultpt60
    HypJetMultpt100
    ) 
    {
        print "Unc for $plot in channel $channel:\n";
        eval{ pdfSum($plot, $channel); };
        print $@ if $@;
        print "\n";
    }
    pdfSumIncl($channel);
}

package UnfoldedResult;

sub new {
    my ($class, $filename) = @_;
    my $self = {};
    $self->{filename} = $filename;
    bless $self, $class;
    $self->{bins} = $self->readFile();
    return $self;
}

sub readFile {
    my $self = shift;
    open my $f, '<', $self->{filename} or die "Opening $self->{filename}\n$!";
    my @bins;
    while (my $line = <$f>) {
        $line =~ s/\$//g;
        my @cols = split / +/, $line;
        push @bins, {-center => $cols[1],
                     -from => $cols[3], -to => $cols[5], 
                     -xsec => $cols[7], -xsecstatunc => $cols[7]};

    }
    return \@bins;
}

package InclusiveResult;

sub new {
    my ($class, $filename) = @_;
    my $self = {};
    $self->{filename} = $filename;
    bless $self, $class;
    $self->{-xsec} = $self->readFile();
    return $self;
}

sub readFile {
    my $self = shift;
    open my $f, '<', $self->{filename} or die "Opening $self->{filename}\n$!";
    my (undef, $syst, undef, $ch, undef, $xsec, undef, $xsecunc) = split / +/, <$f>;
    return $xsec;
}
