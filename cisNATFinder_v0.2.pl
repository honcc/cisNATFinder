#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a revamped-from-scratch perl script to discover and characterize cis natural antisense transcripts, mainly based on antisenseTranscriptDiscoverer_v0.11
#	the main improvement is it will use perl storables extensively.
#
#	Input
#
#		--trnsfrgInfoPlsPath=				path[compulsory]; path of the trnsfrgInfoHsh.pls generated from transfragDoscoverer;
#		--gffPath=							path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=							output directory
#	v0.1
#	[Thu  5 Sep 2013 19:07:54 CEST] cleaned;
#
#	V0.2
#	[11/10/2013 13:47] will take transfrgInfo storable identified using transfragDiscoverer instead of from queryCovPlsIndexPath;
#	[11/10/2013 13:54] read5End25OffsetPlsIndexPath option removed;
#	[11/10/2013 13:55] polyAClusterStorablePath option removed;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-10-11 20:27]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/antisenseTranscripts/cisNATFinder/v0.2/cisNATFinder_v0.2.pl --trnsfrgInfoPlsPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/transfragDiscoverer/storable/trnsfrgInfoHsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff --outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/antisenseTranscripts/cisNATFinder/v0.2/EHI_Standard_CL0.95/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/antisenseTranscripts/cisNATFinder/v0.2/cisNATFinder_v0.2.pl
#	--trnsfrgInfoPlsPath=/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.95/transfragDiscoverer/storable/trnsfrgInfoHsh.pls
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/gff/forPileupCounter.gff
#	--outDir=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/antisenseTranscripts/cisNATFinder/v0.2/EHI_Standard_CL0.95/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1218, readParameters|1452
#	secondaryDependOnSub: currentTime|634
#
#<section ID="startingTasks" num="0">
#----------Read parameters ----------#
&printCMDLogOrFinishMessage("CMDLog");#->1218
my ($trnsfrgInfoPlsPath, $gffPath, $outDir) = &readParameters();#->1452
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedPath" num="1">
my @fileToCheckAry = ();
my $GCPath = "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.GC.tdf"; push @fileToCheckAry, $GCPath;
my $reptvPath = "/Volumes/A_MPro2TB/NGS/IGVInfo/infoTracks/EHI/GC.Repetitiveness/win.100_step.10.Reptv.tdf"; push @fileToCheckAry, $reptvPath;
my $plusCovTDFPath = "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.99/correctedWig.upr/corrected.plus.wig.tdf"; push @fileToCheckAry, $plusCovTDFPath;
my $minusCovTDFPath = "/Volumes/A_MPro2TB/NGS/analyses/EHI_NAT/artefactCorrectedPileup/EHI_Standard/modelFit/W200.S100.I1.R0.10000000/modelFitData.CL0.99/correctedWig.upr/corrected.minus.wig.tdf"; push @fileToCheckAry, $minusCovTDFPath;
my $IGVGenomeID = "EHI_v3.0";
foreach (@fileToCheckAry) {die "Can't read $_" if not -s $_;}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="2">
my $maxThread = 12;
my $minOvrlpSize = 50; #---minimum to be qualified as an overlapping transfrag
my $trnsfrgGffGeneLineOnly = 'yes'; #---only print the gene line in the GFF for trnsfrg
my $minIsolatedDist = 50; #---minimum isolated distance to be regarded as lncRNA

################---[16/10/2013 12:53] not recommended to change, 5 is the best so far for calculating bias index in amoeba, and some other script depends on it, including polyATailFinder ################################
my $numRegion = 5;#---[11/10/2013 16:08] num of the region to be divided in a gene for coverage hit string
###########################################################################################################################################################################################################################

my $minGeneLenRegPlot = 1000; #----will only include the gene longer than this length in the coverage plot
my $resultDirTag = "ID$minIsolatedDist.OS$minOvrlpSize.RG$numRegion.ML$minGeneLenRegPlot";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="3">
my @mkDirAry  = ();
my $resultGFFDir = "$outDir/$resultDirTag/gff/"; push @mkDirAry, $resultGFFDir;
my $resultLogDir = "$outDir/$resultDirTag/log/"; push @mkDirAry, $resultLogDir;
my $resultWigDir = "$outDir/$resultDirTag/wig/"; push @mkDirAry, $resultWigDir;
my $resultXMLDir = "$outDir/$resultDirTag/xml/"; push @mkDirAry, $resultXMLDir;
my $resultStorableDir = "$outDir/$resultDirTag/storable/"; push @mkDirAry, $resultStorableDir;
my $ggplotDirHsh_ref = {};
my @ggplotFileTypeAry = qw /dat pdf R log/;
foreach my $fileType (@ggplotFileTypeAry) {$ggplotDirHsh_ref->{$fileType} = "$outDir/$resultDirTag/ggplot/$fileType"; push @mkDirAry, $ggplotDirHsh_ref->{$fileType};}
system ("mkdir -pm 777 $_") foreach (@mkDirAry);
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="4">
my $allTrnsfrgGFFPrefix = "$resultGFFDir/all.transfrag";
my $lncRNAGFFPrefix = "$resultGFFDir/lncRNA.transfrag";
my $antisenseTrnfrgGFFPrefix = "$resultGFFDir/antisense.transfrag";
my $senseTrnfrgGFFPrefix = "$resultGFFDir/sense.transfrag";
my $XMLPath = "$resultXMLDir/cisNAT.igv.xml";
my $mRNATrnsfrgOvrlpInfoLogPath = "$resultLogDir/mRNATrnsfrgOvrlpInfo.log.xls";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_processInputData
#	primaryDependOnSub: checkGeneInfo|256, getCtgryGeneInfo|680, printGFF_oneRNAPerGene|1251, readGFF_oneRNAPerGene|1353, reportStatus|1485
#	secondaryDependOnSub: currentTime|634, reportStatus|1485
#
#<section ID="processInputData" num="5">
#----------Read Gff
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->1353
&checkGeneInfo($geneInfoHsh_ref);#->256

&reportStatus("retrieving trnsfrgInfoPlsPath", 10, "\n");#->1485
my $trnsfrgInfoHsh_ref = retrieve($trnsfrgInfoPlsPath);
&printGFF_oneRNAPerGene($trnsfrgInfoHsh_ref, $allTrnsfrgGFFPrefix, $trnsfrgGffGeneLineOnly);#->1251

#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/bfmRNA/;
my $mRNAInfoHsh_ref = &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->680
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_defineOverlapping
#	primaryDependOnSub: checkTrnsfrgProximityAndDefinelncRNA|518, checkmRNATrnsfrgOverlap|586, pairTransfragAndmRNA|954, plotRegionHitAndBiasIndex|1156, printmRNATrnsfrgOvrlpInfo|1312
#	secondaryDependOnSub: checkOverlapAndProximity_withMargin|282, ggplotBarChart|748, ggplotHistogram|782, reportStatus|1485
#
#<section ID="defineOverlapping" num="6">
#----------Find overlapping and proximity between mRNA and trnsfrg
my ($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref) = &checkmRNATrnsfrgOverlap($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $resultStorableDir);#->586
my ($mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref, $regionHitHsh_ref) = &pairTransfragAndmRNA($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref, $numRegion, $minOvrlpSize, $minGeneLenRegPlot, $resultStorableDir);#->954
&printmRNATrnsfrgOvrlpInfo($mRNATrnsfrgOvrlpInfoHsh_ref, $mRNATrnsfrgOvrlpInfoLogPath);#->1312
&plotRegionHitAndBiasIndex($regionHitHsh_ref, $ggplotDirHsh_ref);#->1156
my ($hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref, $lncRNAInfoHsh_ref) = &checkTrnsfrgProximityAndDefinelncRNA($geneInfoHsh_ref, $trnsfrgInfoHsh_ref, $minIsolatedDist, $resultStorableDir);#->518
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_outputFiles
#	primaryDependOnSub: outputIGVXML|853, printGFF_oneRNAPerGene|1251
#	secondaryDependOnSub: >none
#
#<section ID="outputFiles" num="7">
&printGFF_oneRNAPerGene($antisenseTrnsfrgInfoHsh_ref, $antisenseTrnfrgGFFPrefix, $trnsfrgGffGeneLineOnly);#->1251
&printGFF_oneRNAPerGene($senseTrnsfrgInfoHsh_ref, $senseTrnfrgGFFPrefix, $trnsfrgGffGeneLineOnly);#->1251
&printGFF_oneRNAPerGene($lncRNAInfoHsh_ref, $lncRNAGFFPrefix, $trnsfrgGffGeneLineOnly);#->1251
&outputIGVXML($XMLPath, $GCPath, $reptvPath, $gffPath, $plusCovTDFPath, $minusCovTDFPath, $IGVGenomeID, $lncRNAGFFPrefix, $antisenseTrnfrgGFFPrefix, $senseTrnfrgGFFPrefix, $allTrnsfrgGFFPrefix);#->853
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|1218
#	secondaryDependOnSub: currentTime|634
#
#<section ID="finishingTasks" num="8">
#system ("open $resultXMLDir");
close DEBUGLOG;
&printCMDLogOrFinishMessage("finishMessage");#->1218
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	XML [n=1]:
#		outputIGVXML
#
#	general [n=7]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printCMDLogOrFinishMessage, readGFF_oneRNAPerGene, readParameters
#		reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, getCtgryGeneInfo, readGFF_oneRNAPerGene
#
#	ggplot [n=2]:
#		ggplotBarChart, ggplotHistogram
#
#	multithread [n=1]:
#		generateThreadHshWithRandomCntg
#
#	plot [n=1]:
#		plotRegionHitAndBiasIndex
#
#	plotInR [n=2]:
#		ggplotBarChart, ggplotHistogram
#
#	range [n=2]:
#		checkOverlapAndProximity_withMargin, checkmRNATrnsfrgOverlap
#
#	reporting [n=1]:
#		currentTime
#
#	specific [n=2]:
#		pairTransfragAndmRNA, plotRegionHitAndBiasIndex
#
#	storable [n=1]:
#		getIndivCntgCovPlsPath
#
#	unassigned [n=3]:
#		checkTrnsfrgProximityAndDefinelncRNA, printGFF_oneRNAPerGene, printmRNATrnsfrgOvrlpInfo
#
#====================================================================================================================================================#

sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|1485
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|149
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: none
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 156
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 20, "\n");#->1485
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 20, "\n");#->1485
	}
}
sub checkOverlapAndProximity_withMargin {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: currentTime|634, generateThreadHshWithRandomCntg|653, reportStatus|1485
#	appearInSub: checkTrnsfrgProximityAndDefinelncRNA|518, checkmRNATrnsfrgOverlap|586
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_defineOverlapping|169
#	input: $checkPrxmty, $maxThread, $qryInfoHsh_ref, $qryMargin, $qryRngType, $refInfoHsh_ref, $refMargin, $refRngType, $reportExactMatch
#	output: $hitAndPrxmtyByQryHsh_ref, $hitAndPrxmtyByRefHsh_ref
#	toCall: my ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref) = &checkOverlapAndProximity_withMargin($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);
#	calledInLine: 558, 626
#....................................................................................................................................................#
	#---incoming variables
	my ($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin) = @_;

	#---outgoing variables
	my $refGeneNumTotal = 0;
	
	#---make a tmpHsh to contain all cntgs the have either ref and qry
	my $tmpCntgHsh_ref = {};
	my $refCntgHsh_ref = {};
	my $qryCntgHsh_ref = {};
	
	foreach my $refGeneID (keys %{$refInfoHsh_ref}) {
		if ($refInfoHsh_ref->{$refGeneID}{$refRngType}) {
			@{$refInfoHsh_ref->{$refGeneID}{$refRngType}} = sort {$a <=> $b} @{$refInfoHsh_ref->{$refGeneID}{$refRngType}};
			$refGeneNumTotal++;
			$tmpCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}++;
			$refCntgHsh_ref->{$refInfoHsh_ref->{$refGeneID}{'cntg'}}{$refGeneID}++;
		}
	}
	
	foreach my $qryGeneID (keys %{$qryInfoHsh_ref}) {
		if ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}) {
			@{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}} = sort {$a <=> $b} @{$qryInfoHsh_ref->{$qryGeneID}{$qryRngType}};
			$tmpCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}++;
			$qryCntgHsh_ref->{$qryInfoHsh_ref->{$qryGeneID}{'cntg'}}{$qryGeneID}++;
		}
	}
	
	my $totalCntgNum = keys %{$tmpCntgHsh_ref};
	my @cntgAry = keys %{$tmpCntgHsh_ref};

	my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->653
	my $refGeneNumProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->1485

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
	
				my $hitAndPrxmtyByRefHsh_InThr_ref = {};
				my $hitAndPrxmtyByQryHsh_InThr_ref = {};

				foreach my $cntg (@{$cntgAry_ref}) {

					#---update the on-screen progress
					#&reportStatus("Finding overlaping on $cntg", 20, "\r");#->1485
		
					my $tmpPrxmtyByRefHsh_ref = {};
					my $tmpPrxmtyByQryHsh_ref = {};

					if ((exists $qryCntgHsh_ref->{$cntg}) and (exists $refCntgHsh_ref->{$cntg})) {#---if there are both ref and qry can both be found on cntg
						foreach my $refGeneID (keys %{$refCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of refGff
							my ($refStart, $refEnd) = ($refInfoHsh_ref->{$refGeneID}{$refRngType}->[0]-$refMargin, $refInfoHsh_ref->{$refGeneID}{$refRngType}->[-1]+$refMargin);
				
							$refGeneNumProc++;
				
							&reportStatus("$refGeneNumProc of $refGeneNumTotal reference genes checked", 20, "\r");#->1485

							foreach my $qryGeneID (keys %{$qryCntgHsh_ref->{$cntg}}) {#--- all ftur on the $strnd of $cntg of QryGtf
	
								my $samestrnd = "no";
								$samestrnd = "yes" if ($refInfoHsh_ref->{$refGeneID}{'strnd'} eq $qryInfoHsh_ref->{$qryGeneID}{'strnd'});
								my ($qryStart, $qryEnd) = ($qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[0]-$qryMargin, $qryInfoHsh_ref->{$qryGeneID}{$qryRngType}->[-1]+$qryMargin);

								my $scene;
								my $ovrlpSize;

								if (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0
									$scene = 'exactMatch';
									$ovrlpSize = $qryEnd - $qryStart;
						
								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1
									$scene = 'overlapTail';
									$ovrlpSize = $refEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2
									$scene = 'overlapHead';
									$ovrlpSize = $qryEnd - $refStart;

								} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3
									$scene = 'cover';
									$ovrlpSize = $qryEnd - $qryStart;

								} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4
									$scene = 'within';
									$ovrlpSize = $refEnd - $refStart;

								#------Proximity with ref's tail proximal to qry's head
								} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

									$scene = 'prxmtyTail';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $qryStart - $refEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"T"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"H"}{$refGeneID} = $tmpPrmxty;
										}
									}

								#------Proximity with ref's head proximal to qry's tail
								} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

									$scene = 'prxmtyHead';

									if ($checkPrxmty eq "yes") {
										my $tmpPrmxty = $refStart - $qryEnd;
										$tmpPrxmtyByRefHsh_ref->{'XS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
										$tmpPrxmtyByQryHsh_ref->{'XS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;

										if ($samestrnd eq "yes") {
											$tmpPrxmtyByRefHsh_ref->{'SS'}{$refGeneID}{"H"}{$qryGeneID} = $tmpPrmxty;
											$tmpPrxmtyByQryHsh_ref->{'SS'}{$qryGeneID}{"T"}{$refGeneID} = $tmpPrmxty;
										}
									}

								} else {#---BUG! possibly other scene?
									#print "[".&currentTime()."] refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";#->634
									die "Unexpected overlapping scene between $refGeneID and $qryGeneID. It's a Bug. Program qutting.\n";
								}
					
								if ($scene ne 'prxmtyTail' and $scene ne 'prxmtyHead' and not ($reportExactMatch eq 'no' and $scene eq 'exactMatch')) {

									@{$hitAndPrxmtyByRefHsh_InThr_ref->{'XS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
									@{$hitAndPrxmtyByQryHsh_InThr_ref->{'XS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);

									if ($samestrnd eq "yes") {
										@{$hitAndPrxmtyByRefHsh_InThr_ref->{'SS'}{'hit'}{$refGeneID}{$qryGeneID}} = ($scene, $ovrlpSize);
										@{$hitAndPrxmtyByQryHsh_InThr_ref->{'SS'}{'hit'}{$qryGeneID}{$refGeneID}} = ($scene, $ovrlpSize);
									}
								}
							}
						}
					}

					#---find the closest proximity for all refs
					if ($checkPrxmty eq "yes") {
						my $refQryRefHsh_ref = {};

						$refQryRefHsh_ref->{'ref'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByRefHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'cntgHsh_ref'} = $refCntgHsh_ref;
						$refQryRefHsh_ref->{'ref'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByRefHsh_InThr_ref;

						$refQryRefHsh_ref->{'qry'}{'tmpPrxmtyHsh_ref'} = $tmpPrxmtyByQryHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'cntgHsh_ref'} = $qryCntgHsh_ref;
						$refQryRefHsh_ref->{'qry'}{'hitAndPrxmtyHsh_ref'} = $hitAndPrxmtyByQryHsh_InThr_ref;
			
						foreach my $refOrQry ('ref', 'qry') {

							my $cntgHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'cntgHsh_ref'};
							my $tmpPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'tmpPrxmtyHsh_ref'};
							my $hitAndPrxmtyHsh_ref = $refQryRefHsh_ref->{$refOrQry}{'hitAndPrxmtyHsh_ref'};

							foreach my $ftur (keys %{$cntgHsh_ref->{$cntg}}) {
								foreach my $XSOrSS ('XS', 'SS') {
									foreach my $HOrT ('H', 'T') {
										$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{"edge"} = -999 if (not exists $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT});
										foreach my $otherFtur (sort {$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$a} <=> $tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$b}} keys %{$tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}}) {
											@{$hitAndPrxmtyHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = ($tmpPrxmtyHsh_ref->{$XSOrSS}{$ftur}{$HOrT}{$otherFtur}, $otherFtur);
											last; #---sample the smallest only
										}
									}
								}
							}
						}
					}
				}
				
				return ($hitAndPrxmtyByRefHsh_InThr_ref, $hitAndPrxmtyByQryHsh_InThr_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my %tmpTransferThrDataHsh = ();
	$tmpTransferThrDataHsh{'ref'}{'all'} = {};
	$tmpTransferThrDataHsh{'qry'}{'all'} = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				($tmpTransferThrDataHsh{'ref'}{'thr'}, $tmpTransferThrDataHsh{'qry'}{'thr'}) = $thr->join;
				foreach my $refOrQry (keys %tmpTransferThrDataHsh) {
					my ($allHsh_ref, $thrHsh_ref) = ($tmpTransferThrDataHsh{$refOrQry}{'all'}, $tmpTransferThrDataHsh{$refOrQry}{'thr'});
					foreach my $XSOrSS ('XS', 'SS') {
						foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}}) {
							foreach my $hitftur (keys %{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}}) {
								@{$allHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}} = @{$thrHsh_ref->{$XSOrSS}{'hit'}{$ftur}{$hitftur}};
							}
						}
						
						if ($thrHsh_ref->{$XSOrSS}{'prxmty'}) {
							foreach my $ftur (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}}) {
								foreach my $HOrT (keys %{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}}) {
									@{$allHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} = @{$thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT}} if $thrHsh_ref->{$XSOrSS}{'prxmty'}{$ftur}{$HOrT};
								}
							}
						}
					}
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	print "\n";

	my $hitAndPrxmtyByRefHsh_ref = $tmpTransferThrDataHsh{'ref'}{'all'};
	my $hitAndPrxmtyByQryHsh_ref = $tmpTransferThrDataHsh{'qry'}{'all'};

	return ($hitAndPrxmtyByRefHsh_ref, $hitAndPrxmtyByQryHsh_ref);
}
sub checkTrnsfrgProximityAndDefinelncRNA {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: checkOverlapAndProximity_withMargin|282, reportStatus|1485
#	appearInSub: >none
#	primaryAppearInSection: 6_defineOverlapping|169
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref, $minIsolatedDist, $resultStorableDir, $trnsfrgInfoHsh_ref
#	output: $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref, $hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $lncRNAInfoHsh_ref
#	toCall: my ($hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref, $lncRNAInfoHsh_ref) = &checkTrnsfrgProximityAndDefinelncRNA($geneInfoHsh_ref, $trnsfrgInfoHsh_ref, $minIsolatedDist, $resultStorableDir);
#	calledInLine: 179
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $trnsfrgInfoHsh_ref, $minIsolatedDist, $resultStorableDir) = @_;
	
	my $lncRNAInfoHsh_ref = {};
	
	#################Check the overlapping or retrieve from the storables
	my $hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref = {};
	my $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref = {};
	my $hitAndPrxmtyToTrnsfrgByAllGeneHshStorablePath = "$resultStorableDir/hitAndPrxmtyToTrnsfrgByAllGeneHsh.pls";
	my $hitAndPrxmtyToAllGeneByTrnsfrgHshStorablePath = "$resultStorableDir/hitAndPrxmtyToAllGeneByTrnsfrgHsh.pls";

	if (-s $hitAndPrxmtyToTrnsfrgByAllGeneHshStorablePath and -s $hitAndPrxmtyToAllGeneByTrnsfrgHshStorablePath) {

		&reportStatus("Retrieving hitAndPrxmtyToTrnsfrgByAllGeneHshStorablePath and hitAndPrxmtyToAllGeneByTrnsfrgHshStorablePath", 20, "\n");#->1485
		
		$hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref = retrieve($hitAndPrxmtyToTrnsfrgByAllGeneHshStorablePath);
		$hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref = retrieve($hitAndPrxmtyToAllGeneByTrnsfrgHshStorablePath);

	} else {
		&reportStatus("Checking overlapping and proximty between mRNA and trnsfrgs", 20, "\n");#->1485
		my $refInfoHsh_ref = $geneInfoHsh_ref;
		my $qryInfoHsh_ref = $trnsfrgInfoHsh_ref;
		my $checkPrxmty = 'yes';
		my $reportExactMatch = 'yes';
		my $maxThread = 4;
		my $refRngType = 'geneRng';
		my $qryRngType = 'geneRng';
		my $refMargin = 0;
		my $qryMargin = 0;
		($hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref) = &checkOverlapAndProximity_withMargin($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);#->282
		store($hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $hitAndPrxmtyToTrnsfrgByAllGeneHshStorablePath);
		store($hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref, $hitAndPrxmtyToAllGeneByTrnsfrgHshStorablePath);
	}
	
	#################Pick up lncRNA
	foreach my $trnsfrgID (keys %{$trnsfrgInfoHsh_ref}) {
		
		#---not hitting anything
		if (not exists $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref->{'XS'}{'hit'}{$trnsfrgID}) {
			my $closeToAllGene = 'no';
			my $tmpDistHsh_ref = {};
			
			foreach my $HOrT ('H', 'T') {
				my ($prmxtyDist, $allGeneID) = @{$hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref->{'XS'}{'prxmty'}{$trnsfrgID}{$HOrT}};
				$tmpDistHsh_ref->{$HOrT} = $prmxtyDist.".".$allGeneID;
				$closeToAllGene = 'yes' if $prmxtyDist < $minIsolatedDist;
			}

			if ($closeToAllGene eq 'no') {
				%{$lncRNAInfoHsh_ref->{$trnsfrgID}} = %{$trnsfrgInfoHsh_ref->{$trnsfrgID}};
			}
		}
	}
	
	return ($hitAndPrxmtyToTrnsfrgByAllGeneHsh_ref, $hitAndPrxmtyToAllGeneByTrnsfrgHsh_ref, $lncRNAInfoHsh_ref);
}
sub checkmRNATrnsfrgOverlap {
#....................................................................................................................................................#
#	subroutineCategory: range
#	dependOnSub: checkOverlapAndProximity_withMargin|282, reportStatus|1485
#	appearInSub: >none
#	primaryAppearInSection: 6_defineOverlapping|169
#	secondaryAppearInSection: >none
#	input: $mRNAInfoHsh_ref, $resultStorableDir, $trnsfrgInfoHsh_ref
#	output: $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref
#	toCall: my ($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref) = &checkmRNATrnsfrgOverlap($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $resultStorableDir);
#	calledInLine: 175
#....................................................................................................................................................#

	my ($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $resultStorableDir) = @_;

	my $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref = {};
	my $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref = {};
	
	my $hitAndPrxmtyToTrnsfrgBymRNAHshStorablePath = "$resultStorableDir/hitAndPrxmtyToTrnsfrgBymRNAHsh.pls";
	my $hitAndPrxmtyTomRNAByTrnsfrgHshStorablePath = "$resultStorableDir/hitAndPrxmtyTomRNAByTrnsfrgHsh.pls";
	
	
	if (-s $hitAndPrxmtyToTrnsfrgBymRNAHshStorablePath and -s $hitAndPrxmtyTomRNAByTrnsfrgHshStorablePath) {

		&reportStatus("Retrieving hitAndPrxmtyToTrnsfrgBymRNAHshStorablePath and hitAndPrxmtyTomRNAByTrnsfrgHshStorablePath", 20, "\n");#->1485
		$hitAndPrxmtyTomRNAByTrnsfrgHsh_ref = retrieve($hitAndPrxmtyTomRNAByTrnsfrgHshStorablePath);
		$hitAndPrxmtyToTrnsfrgBymRNAHsh_ref = retrieve($hitAndPrxmtyToTrnsfrgBymRNAHshStorablePath);
		
	} else {

		&reportStatus("Checking overlapping and proximty between mRNA and trnsfrgs", 20, "\n");#->1485

		my $refInfoHsh_ref = $mRNAInfoHsh_ref;
		my $qryInfoHsh_ref = $trnsfrgInfoHsh_ref;
		my $checkPrxmty = 'no';
		my $reportExactMatch = 'yes';
		my $maxThread = 4;
		my $refRngType = 'geneRng';
		my $qryRngType = 'geneRng';
		my $refMargin = 0;
		my $qryMargin = 0;
		($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref) = &checkOverlapAndProximity_withMargin($refInfoHsh_ref, $qryInfoHsh_ref, $checkPrxmty, $reportExactMatch, $maxThread, $refRngType, $qryRngType, $refMargin, $qryMargin);#->282
		store($hitAndPrxmtyTomRNAByTrnsfrgHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHshStorablePath);
		store($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyToTrnsfrgBymRNAHshStorablePath);
	}

	return ($hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref);
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity_withMargin|282, getCtgryGeneInfo|680, getIndivCntgCovPlsPath|714, printCMDLogOrFinishMessage|1218, readGFF_oneRNAPerGene|1353, reportStatus|1485
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|71, 5_processInputData|149, 8_finishingTasks|197
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 420, 698, 709, 743, 1238, 1241, 1246, 1373, 1501
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;

}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: checkOverlapAndProximity_withMargin|282
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 324
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return $randCntgInThreadHsh_ref;

}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|634
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|149
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 164
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	print "[".&currentTime()."] Filtering GFF on cgtry $ctgryStr.\n";#->634
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	print "[".&currentTime()."] $numGene gene filtered on cgtry $ctgryStr.\n";#->634
	
	return $geneCtgryInfoHsh_ref;
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: currentTime|634
#	appearInSub: >none
#	primaryAppearInSection: >none
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 726
#....................................................................................................................................................#
	
	#my $cntgCovInPlsPathHsh_ref = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);#->714
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	print "[".&currentTime()."] pls path of $numCntg contig stored.\n";#->634
	
	return $cntgCovInPlsPathHsh_ref;
}
sub ggplotBarChart {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotRegionHitAndBiasIndex|1156
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_defineOverlapping|169
#	input: $RScriptPath, $XAXis, $YAxis, $dataPath, $extraArg, $extraStatment, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: 
#	toCall: &ggplotBarChart($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width);
#	calledInLine: 1189
#....................................................................................................................................................#
	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($XAXis, $YAxis)), "\n";
	foreach my $XVal (sort keys %{$plotDataHsh_ref}) {
		my $YVal = $plotDataHsh_ref->{$XVal};
		print PLOTDATA join "", (join "\t", ($XVal, $YVal)), "\n";
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "$extraStatment"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, fill=$YAxis)) + geom_bar(stat=\"identity\") $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

	return ();
}
sub ggplotHistogram {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotRegionHitAndBiasIndex|1156
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_defineOverlapping|169
#	input: $RScriptPath, $XAxis, $binWidth, $dataPath, $dataPtMax, $extraArg, $height, $leftXAxisPercentileLimit, $log2OrLinear, $logPath, $pdfPath, $plotAry_ref, $rightXAxisPercentileLimit, $width
#	output: $plotValueAry_ref
#	toCall: my ($plotValueAry_ref) = &ggplotHistogram($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $XAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);
#	calledInLine: 1211
#....................................................................................................................................................#
	
	my ($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $XAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width) = @_;
	
	my $valueStatObj = Statistics::Descriptive::Full->new();
	$valueStatObj->add_data(@{$plotAry_ref});

	my $leftXAxisLimitValue;
	if ($leftXAxisPercentileLimit eq 'min') {
		$leftXAxisLimitValue = $valueStatObj->min();
	} else {
		$leftXAxisLimitValue = $valueStatObj->percentile($leftXAxisPercentileLimit);
	}

	my $rightXAxisLimitValue;
	if ($rightXAxisPercentileLimit eq 'max') {
		$rightXAxisLimitValue = $valueStatObj->max();
	} else {
		$rightXAxisLimitValue = $valueStatObj->percentile($rightXAxisPercentileLimit);
	}
	
	#---trim the end values
	my @trimmedAry = ();
	foreach my $value (@{$plotAry_ref}) {
		my $transformedValue = $value;
		if ($log2OrLinear eq 'log2') {
			$transformedValue = 0;
			eval {$transformedValue = log($value)/log(2)};
		}
		push @trimmedAry, $transformedValue if $value <= $rightXAxisLimitValue and $value >= $leftXAxisLimitValue;
	}
	
	$dataPtMax = @trimmedAry if $dataPtMax > @trimmedAry;

	#---down sample the data point number
	my @shuffleIndexAry = shuffle(0..$#trimmedAry);
	my $plotValueAry_ref = ();
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA $XAxis."\n";
	foreach my $i (0..$dataPtMax-1) {
		my $shuffleValue = $trimmedAry[$shuffleIndexAry[$i]];
		push @{$plotValueAry_ref}, $shuffleValue;
		print PLOTDATA $shuffleValue."\n";
		
	}
	close PLOTDATA;
	
	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAxis)) + ggtitle(\"Distribution of $XAxis $log2OrLinear scale [n=$dataPtMax]\") + geom_histogram(binwidth=$binWidth, aes(y = ..density.., fill = ..count..)) + geom_density() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath ");

	return $plotValueAry_ref;
	
}
sub outputIGVXML {
#....................................................................................................................................................#
#	subroutineCategory: XML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_outputFiles|184
#	secondaryAppearInSection: >none
#	input: $GCPath, $IGVGenomeID, $XMLPath, $allTrnsfrgGFFPrefix, $antisenseTrnfrgGFFPrefix, $gffPath, $lncRNAGFFPrefix, $minusCovTDFPath, $plusCovTDFPath, $reptvPath, $senseTrnfrgGFFPrefix
#	output: none
#	toCall: &outputIGVXML($XMLPath, $GCPath, $reptvPath, $gffPath, $plusCovTDFPath, $minusCovTDFPath, $IGVGenomeID, $lncRNAGFFPrefix, $antisenseTrnfrgGFFPrefix, $senseTrnfrgGFFPrefix, $allTrnsfrgGFFPrefix);
#	calledInLine: 192
#....................................................................................................................................................#

	my ($XMLPath, $GCPath, $reptvPath, $gffPath, $plusCovTDFPath, $minusCovTDFPath, $IGVGenomeID, $lncRNAGFFPrefix, $antisenseTrnfrgGFFPrefix, $senseTrnfrgGFFPrefix, $allTrnsfrgGFFPrefix) = @_;

	my $colorHsh_ref = {};
	$colorHsh_ref->{'plus'} = '255,153,153';
	$colorHsh_ref->{'minus'} = '153,153,255';

	my $wigPrefixHsh_ref = {};
	foreach my $name (keys %{$wigPrefixHsh_ref}) {$wigPrefixHsh_ref->{$name}{'ext'} = 'wig';}

	my $GFFPrefixHsh_ref = {};
	$GFFPrefixHsh_ref->{'antisenseTrnsfrg'}{'prefix'} = $antisenseTrnfrgGFFPrefix;
	$GFFPrefixHsh_ref->{'senseTrnfrg'}{'prefix'} = $senseTrnfrgGFFPrefix;
	$GFFPrefixHsh_ref->{'lncRNA'}{'prefix'} = $lncRNAGFFPrefix;
	$GFFPrefixHsh_ref->{'allTrnsfrg'}{'prefix'} = $allTrnsfrgGFFPrefix;
	foreach my $name (keys %{$GFFPrefixHsh_ref}) {$GFFPrefixHsh_ref->{$name}{'ext'} = 'gff';}
	
	open (XML, ">", $XMLPath);
	printf XML "%0s", "<\?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"\?>\n";
	printf XML "%0s", "<Session genome=\"$IGVGenomeID\" version=\"5\">\n";
	printf XML "%4s", "<Resources>\n";
	
	foreach my $filePath ($GCPath, $reptvPath, $gffPath, $plusCovTDFPath, $minusCovTDFPath) {
		printf XML "%8s", "<Resource path=\"$filePath\"\/>\n";
	}
	
	foreach my $hsh_ref ($wigPrefixHsh_ref, $GFFPrefixHsh_ref) {
		foreach my $name (sort keys %{$hsh_ref}) {
			foreach my $strnd ('plus', 'minus') {
				my $filePath = "$hsh_ref->{$name}{'prefix'}.$strnd.$hsh_ref->{$name}{'ext'}";
				printf XML "%8s", "<Resource path=\"$filePath\"\/>\n";
			}
		}
	}
	
	printf XML "%4s", "</Resources>\n";
	printf XML "%4s", "<Panel height=\"900\" name=\"DataPanel\" width=\"1200\">\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;20.0;10.0;25.0;35.0;153,153,255;255,255,255;255,153,153\" fontSize=\"10\" id=\"$GCPath\" name=\"GC\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track colorScale=\"ContinuousColorScale;2.0;1.0;3.0;15.0;255,255,255;153,153,255;255,153,153\" fontSize=\"10\" id=\"$reptvPath\" name=\"Repetitiveness\" normalize=\"false\" renderer=\"HEATMAP\" height=\"20\" sortable=\"true\" visible=\"true\" windowFunction=\"none\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;17.0;255,255,255;0,0,178\" displayMode=\"SQUISHED\"  fontSize=\"10\" height=\"20\" id=\"EHI_v13_genes\" name=\"bonaFide mRNA\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;40.0;255,255,255;0,0,178\" displayMode=\"COLLAPSED\"  fontSize=\"10\" height=\"20\" id=\"$gffPath\" name=\"All Genomic Features\" renderer=\"GENE_TRACK\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">\n";
	printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"255,153,153\" displayMode=\"COLLAPSED\" fontSize=\"10\" height=\"100\" id=\"$plusCovTDFPath\" name=\"Corrected Plus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	printf XML "%8s", "<Track autoScale=\"true\" color=\"153,153,255\" displayMode=\"COLLAPSED\" fontSize=\"10\" height=\"100\" id=\"$minusCovTDFPath\" name=\"Corrected Minus Cov\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">\n";
	printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
	printf XML "%8s", "</Track>\n";

	foreach my $name (sort keys %{$GFFPrefixHsh_ref}) {
		foreach my $strnd ('plus', 'minus') {
			my $filePath = "$GFFPrefixHsh_ref->{$name}{'prefix'}.$strnd.$GFFPrefixHsh_ref->{$name}{'ext'}";
			my $trackName = "$name.$strnd";
			printf XML "%8s", "<Track autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"$colorHsh_ref->{$strnd}\" displayMode=\"SQUISHED\"  fontSize=\"10\" height=\"20\" id=\"$filePath\" name=\"$trackName\" renderer=\"GENE_TRACK\" sortable=\"false\" visible=\"true\" windowFunction=\"none\">\n";
			printf XML "%12s", "<DataRange type=\"LINEAR\"/>\n";
			printf XML "%8s", "</Track>\n";
		}
	}

	foreach my $name (sort keys %{$wigPrefixHsh_ref}) {
		foreach my $strnd ('plus', 'minus') {
			my $filePath = "$wigPrefixHsh_ref->{$name}{'prefix'}.$strnd.$wigPrefixHsh_ref->{$name}{'ext'}";
			my $trackName = "$name.$strnd";
			printf XML "%8s", "<Track autoScale=\"true\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"$colorHsh_ref->{$strnd}\" fontSize=\"10\" height=\"50\" id=\"$filePath\" name=\"$trackName\" renderer=\"BAR_CHART\" sortable=\"false\" visible=\"true\" windowFunction=\"none\">\n";
			printf XML "%12s", "<DataRange type=\"LOG\"/>\n";
			printf XML "%8s", "</Track>\n";
		}
	}
	
	printf XML "%4s", "</Panel>\n";

	printf XML "%0s", "</Session>\n";

	close XML;

}
sub pairTransfragAndmRNA {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|1485
#	appearInSub: >none
#	primaryAppearInSection: 6_defineOverlapping|169
#	secondaryAppearInSection: >none
#	input: $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref, $mRNAInfoHsh_ref, $minGeneLenRegPlot, $minOvrlpSize, $numRegion, $resultStorableDir, $trnsfrgInfoHsh_ref
#	output: $antisenseTrnsfrgInfoHsh_ref, $mRNATrnsfrgOvrlpInfoHsh_ref, $regionHitHsh_ref, $senseTrnsfrgInfoHsh_ref
#	toCall: my ($mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref, $regionHitHsh_ref) = &pairTransfragAndmRNA($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref, $numRegion, $minOvrlpSize, $minGeneLenRegPlot, $resultStorableDir);
#	calledInLine: 176
#....................................................................................................................................................#
	my ($mRNAInfoHsh_ref, $trnsfrgInfoHsh_ref, $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref, $hitAndPrxmtyTomRNAByTrnsfrgHsh_ref, $numRegion, $minOvrlpSize, $minGeneLenRegPlot, $resultStorableDir) = @_;
	
	&reportStatus("Pairing mRNA and trnsfrgs", 20, "\n");#->1485

	my $mRNATrnsfrgOvrlpInfoHsh_ref = {};
	my $antisenseTrnsfrgInfoHsh_ref = {};
	my $senseTrnsfrgInfoHsh_ref = {};
	my $regionHitHsh_ref = {};

	my $mRNATrnsfrgOvrlpInfoHshPlsPath = "$resultStorableDir/mRNATrnsfrgOvrlpInfoHsh.pls";
	my $antisenseTrnsfrgInfoHshPlsPath = "$resultStorableDir/antisenseTrnsfrgInfoHsh.pls";
	my $senseTrnsfrgInfoHshPlsPath = "$resultStorableDir/senseTrnsfrgInfoHsh.pls";
	my $regionHitHshPlsPath = "$resultStorableDir/regionHitHsh.pls";

	if (-s $mRNATrnsfrgOvrlpInfoHshPlsPath and -s $antisenseTrnsfrgInfoHshPlsPath and -s $senseTrnsfrgInfoHshPlsPath and -s $regionHitHshPlsPath) {
	
		$mRNATrnsfrgOvrlpInfoHsh_ref = retrieve($mRNATrnsfrgOvrlpInfoHshPlsPath);
		$antisenseTrnsfrgInfoHsh_ref = retrieve($antisenseTrnsfrgInfoHshPlsPath);
		$senseTrnsfrgInfoHsh_ref = retrieve($senseTrnsfrgInfoHshPlsPath);
		$regionHitHsh_ref = retrieve($regionHitHshPlsPath);

	} else {

		foreach  my $dirtn (qw/a s/) {
			foreach my $region (1..$numRegion) {
				$regionHitHsh_ref->{$dirtn}{$region} = 0;
			}
		}
		my @paramAry = qw /ovrlpSize covPerNt trnsfrgLen/;
		my $mRNAProc = 0;
		foreach my $geneID (sort keys %{$mRNAInfoHsh_ref}) {

			$mRNAProc++;
			&reportStatus("$mRNAProc mRNA processed", 10, "\r");#->1485

			#---get gene length
			my @CDSRngAry = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$geneID}{'CDSRng'}};
			my ($mRNAStart, $mRNAEnd) = ($CDSRngAry[0], $CDSRngAry[-1]);
			my $mRNALength = $mRNAEnd - $mRNAStart + 1;
			my $tmpInfoHsh_ref = {};
			my $neighborOvrlpHsh_ref = {};

			#---[11/10/2013 19:26] set an empty array reference
			$tmpInfoHsh_ref->{$_}{'trnsfrgID'} = [] foreach (qw/a s/);

			#---trnsfrg exists
			if (exists $hitAndPrxmtyToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$geneID}) {
				foreach my $trnsfrgID (keys %{$hitAndPrxmtyToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$geneID}}) {
					my $tmpParamValueHsh_ref = {};
					(undef, $tmpParamValueHsh_ref->{'ovrlpSize'}) = @{$hitAndPrxmtyToTrnsfrgBymRNAHsh_ref->{'XS'}{'hit'}{$geneID}{$trnsfrgID}};
				
					#---skip the trnsfrg if ovrlpSize < $minOvrlpSize;
					next if $tmpParamValueHsh_ref->{'ovrlpSize'} < $minOvrlpSize;
				
					$tmpParamValueHsh_ref->{'covPerNt'} = $trnsfrgInfoHsh_ref->{$trnsfrgID}{'covPerNt'};
					$tmpParamValueHsh_ref->{'trnsfrgLen'} = $trnsfrgInfoHsh_ref->{$trnsfrgID}{'trnsfrgLen'};
				
					my $dirtn;
					if ($mRNAInfoHsh_ref->{$geneID}{'strnd'} ne $trnsfrgInfoHsh_ref->{$trnsfrgID}{'strnd'}) {
						$dirtn = 'a';
						%{$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}} = %{$trnsfrgInfoHsh_ref->{$trnsfrgID}};
						$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'no';
					} else {
						$dirtn = 's';
						%{$senseTrnsfrgInfoHsh_ref->{$trnsfrgID}} = %{$trnsfrgInfoHsh_ref->{$trnsfrgID}};
						$senseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'no';
					}

					push @{$tmpInfoHsh_ref->{$dirtn}{'trnsfrgID'}}, $trnsfrgID;
				
					#----push the trsnfrgInfo
					foreach my $param (@paramAry) {
						push @{$tmpInfoHsh_ref->{$dirtn}{$param}}, $tmpParamValueHsh_ref->{$param};
					}

					#----check overlapping
					foreach my $otherGeneID (keys %{$hitAndPrxmtyTomRNAByTrnsfrgHsh_ref->{'XS'}{'hit'}{$trnsfrgID}}) {
						if ($otherGeneID ne $geneID) {
							$neighborOvrlpHsh_ref->{$dirtn}++;
							$antisenseTrnsfrgInfoHsh_ref->{$trnsfrgID}{'neighborOvrlp'} = 'yes' if $dirtn eq 'a';
						}
					}
				}
			}
		
			#---summary trnsfrg
			foreach my $dirtn ('a', 's') {

				#---zero the ovrlpPct
				my $trnsfrgNum = @{$tmpInfoHsh_ref->{$dirtn}{'trnsfrgID'}};
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'ovrlpPct'} = 0;
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{"total_ovrlpSize"} = 0;
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'trnsfrgNum'} = $trnsfrgNum;

				#---record the neighborOvrlp
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'neighborOvrlp'} = 0; 
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'neighborOvrlp'} = $neighborOvrlpHsh_ref->{$dirtn} if $neighborOvrlpHsh_ref->{$dirtn};
			
				#---averge the transfrg param
				foreach my $param (@paramAry) { #---@paramAry = (ovrlpSize covPerNt trnsfrgLen)
					$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{"avg_".$param} = 0;
					if ($tmpInfoHsh_ref->{$dirtn}{$param}) {
						#---calculate the average for all param
						$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{"avg_".$param} = sprintf "%.2f", sum(@{$tmpInfoHsh_ref->{$dirtn}{$param}})/@{$tmpInfoHsh_ref->{$dirtn}{$param}};
						#---calculate the ovrlpPct with mRNA
						if ($param eq 'ovrlpSize') {
							$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'ovrlpPct'} = sprintf "%.2f", 100*sum(@{$tmpInfoHsh_ref->{$dirtn}{'ovrlpSize'}})/$mRNALength; 
							$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{"total_ovrlpSize"} = sprintf "%.2f", sum(@{$tmpInfoHsh_ref->{$dirtn}{'ovrlpSize'}});
						}
					}
				}
			
				#---[11/10/2013 15:38] calculate the coverage string
				my @headToTailRngAry = ($mRNAStart..$mRNAEnd);
				@headToTailRngAry = reverse(@headToTailRngAry) if $mRNAInfoHsh_ref->{$geneID}{'strnd'} eq '-';
				my $rltvPos = 0;
				my %posRegConversionHsh = ();
				foreach my $pos (@headToTailRngAry) {
					#---[11/10/2013 16:28] minium 1, maximum $numRegion since $rltvPos is alway < $mRNALength;
					my $region = int($numRegion*($rltvPos/$mRNALength))+1;
					$posRegConversionHsh{$pos} = $region;
					$rltvPos++;
				}
			
				#---[11/10/2013 16:01] push - to every region
				my @regionHitAry = ();
				push @regionHitAry, '-' foreach (1..$numRegion);
			
				#---[11/10/2013 16:14] assing + to the hit region
				foreach my $trnsfrgID (@{$tmpInfoHsh_ref->{$dirtn}{'trnsfrgID'}}) {
					my ($trnsfrgStart, $trnsfrgEnd) = sort {$a <=> $b} @{$trnsfrgInfoHsh_ref->{$trnsfrgID}{'geneRng'}};
					foreach my $pos ($trnsfrgStart..$trnsfrgEnd) {
						if ($posRegConversionHsh{$pos}) {
							my $region = $posRegConversionHsh{$pos}; #---[11/10/2013 16:14] covert to region number
							my $i = $region - 1;
							$regionHitAry[$i] = '+';
						}
					}
				}
				my $regHitStr = join "", @regionHitAry;
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'regHitStr'} = $regHitStr;

				#---[11/10/2013 18:28] calculate the bias index
				my $biasIndex = 0;
				my $plusRegCount = 0;
				my $minusRegCount = 0;
				foreach  my $i (0..$#regionHitAry) {
					my $score = $i - int($numRegion/2);
					$score++ if ($score >= 0 and $numRegion % 2 == 0);#---[11/10/2013 18:43] skip the zero score if numRegion is a even number
					if ($regionHitAry[$i] eq '+') {
						$biasIndex += $score;
						$plusRegCount++;
					} else {
						$minusRegCount++;
					}
				}
				
				$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{'biasIndex'} = $biasIndex;
			
				#---[11/10/2013 16:16] put the region into regionHitHsh_ref 
				if ($mRNALength >= $minGeneLenRegPlot) {
					
					#---[11/10/2013 19:01] take genes that are not fully covered
					$regionHitHsh_ref->{$dirtn}{'totalGeneNum'}++;
					if ($plusRegCount == $numRegion) {
						$regionHitHsh_ref->{$dirtn}{'fullyCover'}++;
					} elsif ($minusRegCount == $numRegion) {
						$regionHitHsh_ref->{$dirtn}{'noTransfrg'}++;
					} else {
						$regionHitHsh_ref->{$dirtn}{'partiallyCover'}++;
						push @{$regionHitHsh_ref->{$dirtn}{'biasIndex'}}, $biasIndex;
					}
					
					foreach  my $i (0..$#regionHitAry) {
						my $region = $i+1;
						$regionHitHsh_ref->{$dirtn}{'count'}{$region}++ if $regionHitAry[$i] eq '+';
					}
				}
			}
		}
		
		store($mRNATrnsfrgOvrlpInfoHsh_ref, $mRNATrnsfrgOvrlpInfoHshPlsPath);
		store($antisenseTrnsfrgInfoHsh_ref, $antisenseTrnsfrgInfoHshPlsPath);
		store($senseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHshPlsPath);
		store($regionHitHsh_ref, $regionHitHshPlsPath);
	
	}
	
	return ($mRNATrnsfrgOvrlpInfoHsh_ref, $antisenseTrnsfrgInfoHsh_ref, $senseTrnsfrgInfoHsh_ref, $regionHitHsh_ref);
}
sub plotRegionHitAndBiasIndex {
#....................................................................................................................................................#
#	subroutineCategory: specific, plot
#	dependOnSub: ggplotBarChart|748, ggplotHistogram|782
#	appearInSub: >none
#	primaryAppearInSection: 6_defineOverlapping|169
#	secondaryAppearInSection: >none
#	input: $ggplotDirHsh_ref, $regionHitHsh_ref
#	output: 
#	toCall: &plotRegionHitAndBiasIndex($regionHitHsh_ref, $ggplotDirHsh_ref);
#	calledInLine: 178
#....................................................................................................................................................#
	my ($regionHitHsh_ref, $ggplotDirHsh_ref) = @_;

	foreach my $dirtn (keys %{$regionHitHsh_ref}) {
		my $totalGeneNum = $regionHitHsh_ref->{$dirtn}{'totalGeneNum'};
		my $regHitDataHsh_ref = {};
		foreach my $region (keys %{$regionHitHsh_ref->{$dirtn}{'count'}}) {
			$regHitDataHsh_ref->{$region} = 100*$regionHitHsh_ref->{$dirtn}{'count'}{$region}/$totalGeneNum;
		}
		
		{
			my $plotDataHsh_ref = $regHitDataHsh_ref;
			my $nameTag = "region.hit.dirtn.$dirtn";
			my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $ggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $XAXis = 'region';
			my $YAxis = 'percentageOfmRNA';
			my $extraStatment = "";
			my $extraArg = "+ ylim(0,100) + ggtitle(\"Region covered by transfrag in $dirtn dirtn N=$totalGeneNum\")";
			my $height = 5;
			my $width = 10;
			&ggplotBarChart($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width);#->748
		}
		
		{
			my $plotAry_ref = $regionHitHsh_ref->{$dirtn}{'biasIndex'};
			my $partiallyCoverNum = $regionHitHsh_ref->{$dirtn}{'partiallyCover'};
			my $fullyCoverNum = $regionHitHsh_ref->{$dirtn}{'fullyCover'};
			my $noTransfrgNum = $regionHitHsh_ref->{$dirtn}{'noTransfrg'};
			my $nameTag = "biasIndex.dirtn.$dirtn";
			my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $ggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $leftXAxisPercentileLimit = 'min';
			my $rightXAxisPercentileLimit = 'max';
			my $XAxis = "biasIndex";
			my $binWidth = 1;
			my $extraArg = " + ggtitle(\"bias index of partially covered genes N=$partiallyCoverNum \[fullyCover=$fullyCoverNum, noTransfrg=$noTransfrgNum\]\") + scale_fill_gradient(\"count\", low = \"grey\", high = \"red\")";
			my $log2OrLinear = 'linear';#-----choose the value to be plot in linear or log2 scale here
			my $dataPtMax = 99999999;
			my $height = 8;
			my $width = 24;
			&ggplotHistogram($plotAry_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $leftXAxisPercentileLimit, $rightXAxisPercentileLimit, $XAxis, $binWidth, $dataPtMax, $extraArg, $log2OrLinear, $height, $width);#->782
		}
	}

	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|634
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|71, 8_finishingTasks|197
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 77, 204
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->634
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->634
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->634
		print "=========================================================================\n\n";
	}
}
sub printGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|149, 7_outputFiles|184
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref, $gffGeneLineOnly, $outGFFPrefix
#	output: none
#	toCall: &printGFF_oneRNAPerGene($geneInfoHsh_ref, $outGFFPrefix, $gffGeneLineOnly);
#	calledInLine: 160, 189, 190, 191
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $outGFFPrefix, $gffGeneLineOnly) = @_;
	
	$outGFFPrefix =~ s/\/+/\//g;
	
	open (GFFOUTALL, ">$outGFFPrefix.both.gff");
	open (GFFOUTPLUS, ">$outGFFPrefix.plus.gff");
	open (GFFOUTMINUS, ">$outGFFPrefix.minus.gff");
	print GFFOUTALL "##gff-version\t3\n";
	print GFFOUTPLUS "##gff-version\t3\n";
	print GFFOUTMINUS "##gff-version\t3\n";

	foreach my $geneID (sort {$a cmp $b} keys %{$geneInfoHsh_ref}) {
		my $cntg = $geneInfoHsh_ref->{$geneID}{"cntg"};
		my $strnd = $geneInfoHsh_ref->{$geneID}{"strnd"};
		my $ctgry = $geneInfoHsh_ref->{$geneID}{"ctgry"};
		my $description = $geneInfoHsh_ref->{$geneID}{"description"};

		my ($geneStart, $geneEnd) = @{$geneInfoHsh_ref->{$geneID}{'geneRng'}};
		print GFFOUTALL join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";
		if ($strnd eq "+") {
			print GFFOUTPLUS join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";
		}
		if ($strnd eq "-") {
			print GFFOUTMINUS join "", (join "\t", ($cntg, 'BCP', 'gene', $geneStart, $geneEnd, ".", $strnd, ".", "ID=$geneID;Name=$geneID;description=$description")), "\n";
		}
		if ($gffGeneLineOnly eq 'no') {#-----will print the RNA and exon lines also, aim to avoid display annoyance on IGV
			my $RNAID = $geneInfoHsh_ref->{$geneID}{"RNAID"};
			my ($RNAStart, $RNAEnd) = @{$geneInfoHsh_ref->{$geneID}{'RNARng'}};
			print GFFOUTALL join "", (join "\t", ($cntg, 'BCP', $ctgry, $RNAStart, $RNAEnd, ".", $strnd, ".", "ID=$RNAID;Name=$geneID;Parent=$geneID;description=$description;")), "\n";
			
			foreach my $rngType ('exon', 'CDS') {
				if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
					my $num = 1;
					my @rng = @{$geneInfoHsh_ref->{$geneID}{$rngType}};
					for (my $i=0; $i < $#rng; $i += 2) {
						$num++;
						my ($start, $end) = ($rng[$i], $rng[$i+1]);
						my $ID = "$rngType\_$num\_$RNAID";
						print GFFOUTALL join "", (join "\t", ($cntg, "BCP", $rngType, $start, $end, ".", $strnd, ".", "ID=$ID;Name=$ID;Parent=$RNAID;description=.;")), "\n";
					}
				}
			}
		}
	}
	close GFFOUTALL;
	close GFFOUTPLUS;
	close GFFOUTMINUS;
}
sub printmRNATrnsfrgOvrlpInfo {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 6_defineOverlapping|169
#	secondaryAppearInSection: >none
#	input: $mRNATrnsfrgOvrlpInfoHsh_ref, $mRNATrnsfrgOvrlpInfoLogPath
#	output: none
#	toCall: &printmRNATrnsfrgOvrlpInfo($mRNATrnsfrgOvrlpInfoHsh_ref, $mRNATrnsfrgOvrlpInfoLogPath);
#	calledInLine: 177
#....................................................................................................................................................#

	my ($mRNATrnsfrgOvrlpInfoHsh_ref, $mRNATrnsfrgOvrlpInfoLogPath) = @_;
	
	open (OUTFILE, ">", $mRNATrnsfrgOvrlpInfoLogPath);
	foreach my $geneID (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref}) {
		my @outLineAry = ();
		push @outLineAry, 'geneID';
		foreach my $dirtn (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}}) {
			foreach my $param (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}}) {
				my $column = join '_', ($dirtn, $param);
				push @outLineAry, $column;
			}
		}
		print OUTFILE join "", ((join "\t", @outLineAry), "\n");
		last;
	}
	foreach my $geneID (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref}) {
		my @outLineAry = ();
		push @outLineAry, $geneID;
		foreach my $dirtn (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}}) {
			foreach my $param (sort keys %{$mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}}) {
				push @outLineAry, $mRNATrnsfrgOvrlpInfoHsh_ref->{$geneID}{$dirtn}{$param};
			}
		}
		print OUTFILE join "", ((join "\t", @outLineAry), "\n");
	}
	close OUTFILE;

}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|634
#	appearInSub: >none
#	primaryAppearInSection: 5_processInputData|149
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 155
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->634
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|71
#	secondaryAppearInSection: >none
#	input: none
#	output: $gffPath, $outDir, $trnsfrgInfoPlsPath
#	toCall: my ($trnsfrgInfoPlsPath, $gffPath, $outDir) = &readParameters();
#	calledInLine: 78
#....................................................................................................................................................#
	
	my ($trnsfrgInfoPlsPath, $gffPath, $outDir);
	my $dirPath = dirname(rel2abs($0));

	$outDir = "$dirPath/cisNATFinder/";

	GetOptions 	("trnsfrgInfoPlsPath=s" => \$trnsfrgInfoPlsPath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($trnsfrgInfoPlsPath, $gffPath) {
		die "Can't read $fileToCheck" if (not -s $fileToCheck and not -s "$fileToCheck.gz");
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($trnsfrgInfoPlsPath, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|634
#	appearInSub: checkGeneInfo|256, checkOverlapAndProximity_withMargin|282, checkTrnsfrgProximityAndDefinelncRNA|518, checkmRNATrnsfrgOverlap|586, pairTransfragAndmRNA|954
#	primaryAppearInSection: 5_processInputData|149
#	secondaryAppearInSection: 5_processInputData|149, 6_defineOverlapping|169
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 158, 270, 278, 331, 345, 356, 542, 548, 609, 615, 967, 998
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->634

	return ();
}

exit;
