#!/usr/bin/perl -w

# % cat UnfoldingResults/PDF_CENTRAL/Nominal/combined/HypToppTResults.txt
# XAxisbinCenters[bin]: 35 bin: 0 to 65 DiffXsec: 0.00403064 StatError: 5.86689e-05 GenDiffXsec: 0.980369
#XAxisbinCenters[bin]: 100 bin: 65 to 125 DiffXsec: 0.00663165 StatError: 6.96881e-05 GenDiffXsec: 1.53973
# XAxisbinCenters[bin]: 165 bin: 125 to 200 DiffXsec: 0.00367398 StatError: 2.49495e-05 GenDiffXsec: 0.834946
# XAxisbinCenters[bin]: 245 bin: 200 to 290 DiffXsec: 0.00065152 StatError: 6.52835e-06 GenDiffXsec: 0.256277
# XAxisbinCenters[bin]: 350 bin: 290 to 400 DiffXsec: 5.23712e-05 StatError: 6.98898e-07 GenDiffXsec: 0.0563997


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
    my $nominalIncl = InclusiveResult->new("Plots/PDF_CENTRAL/Nominal/$channel/InclusiveXSec.txt");
    my $nominal = InclusiveResult->new("Plots/Nominal/$channel/InclusiveXSec.txt");
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

#    storeInclXsec("PDF_UP", $channel, $nominalIncl->{-xsec} + sqrt($upInclSUM2));
#    storeInclXsec("PDF_DOWN", $channel, $nominalIncl->{-xsec} - sqrt($downInclSUM2));    
    storeInclXsec("PDF_UP", $channel, $nominal->{-xsec} * (1.0 + sqrt($upInclSUM2)/$nominalIncl->{-xsec}));
    storeInclXsec("PDF_DOWN", $channel, $nominal->{-xsec} * (1.0  - sqrt($downInclSUM2)/$nominalIncl->{-xsec}));
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
