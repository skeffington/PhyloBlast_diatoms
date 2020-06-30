#For blasting against the MMETSP transcriptomes. First make sure appropriate PATH and BLASTDB variables are set and that the blast module and perl modules are loaded.
#Note that the evalue threshold chosen will effect the qcov of the output. qcov reflects all HSPs found for that hit.
#[input taxa: MMETSP] is of the form: "MMETSP1471	Pycnococcus provasolii	RCC733"
#[input taxa: taxonomy] is of the form: "Pycnococcus provasolii	131567	2759	33090	3041	3152	41878	41879	41880" Note goes from top level downwards.
#[input: phyla name to number] is of the form: "Bacillariophyta	2836"

use strict;
use warnings;

my $usage = "perl mmetspblast.pl [input fasta] [input taxa: MMETSP] [input taxa: taxonomy] [input: phyla name to number] [evalue threshold] [qcov cutoff] [output prefix]\n";

my $fastain = shift or die $usage;
my $mmetspin = shift or die $usage;
my $taxonomyin = shift or die $usage;
my $phylain = shift or die $usage;
my $evalue = shift or die $usage;
my $qcovco = shift or die $usage;
my $outpre = shift or die $usage;

open (FAIN, '<', $fastain);
open (MMIN, '<', $mmetspin);
open (TAIN, '<', $taxonomyin);
open (PYIN, '<', $phylain);

my %TAXA;

#Read in species name - MMETSP id pairs into %TAXA hash: $TAXA{$mmid}{$spname}. Include counts of transcriptomes for each species. DOESN'T

while (my $line = <MMIN>){
	if ($line =~ /^([^\t]*)\t([^\t]*)/){
		chomp $line;
		my $mmid = $1;
		my $spname = $2;
		$TAXA{$mmid}{$spname} = '1';
		#print "mmid is $mmid and species is $spname\n";
		#if (exists $TAXA{$mmid}{$spname} ){
		#	my $count = $TAXA{$mmid}{$spname};
		#	$count++;
		#	$TAXA{$mmid}{$spname} = $count;
		#	}else{
		#	$TAXA{$mmid}{$spname} = '1';
		#}
	}
}

#Read in the taxonomy for each species: $PHYLA{species name} = phyla. The value is first the NCBI taxonomy number, and then converted to a name in the second loop

my %PHYLA;

while (my $line = <TAIN>){
	if ($line =~ /^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/){
		chomp $line;
		my $name = $1;
		my $l1 = $2;
		my $l2 = $3;
		my $l3 = $4;
		my $l4 = $5;
		$PHYLA{$name}= $l4;
		}
}

while (my $line = <PYIN>){
	if ($line =~ /^([^\t]*)\t([^\t]*)/){
		my $phyla = $1;
		my $number = $2;
		foreach my $species (keys %PHYLA){
			my $id = $PHYLA{$species};
			if ($id =~ /\d/){ # to prevent next loop throwing warnings
				if ($id == $number){
					$PHYLA{$species} = $phyla;
				}
			}
		}
	}
}

#Generate hash with the correct structure for ordering the blast searches

my %TAXON;
my %TOMEC; #transcriptome counts

foreach my $species (keys %PHYLA){
	#print "species is $species\n";
	my $phylum = $PHYLA{$species};
	$TAXON{$phylum}{$species} = ();#avoiding void context
	my $count = '0';
	foreach my $mmetid (keys %TAXA){
		my $spname;
		foreach my $sp (keys %{ $TAXA{$mmetid} }){
			$spname = $sp;
			}
		#print "\t spname is $spname\n";
		if ($spname eq $species){
			$TAXON{$phylum}{$species}{$mmetid} = ();
			$count++;
		}
	}
	$TOMEC{$phylum}{$species}=$count;
	#print "TOMEC entry for $phylum $species is $TOMEC{$phylum}{$species} \n";
}

#Generate output files. summary contains the number of hits per transcriptome for each diatom and for the broader classes. 
#top contains the top hit for each transcriptome. All contains all hits for each transcriptome

my $allout = $outpre.'_allhits.txt';
my $topout = $outpre.'_tophits.txt';
my $sumout = $outpre.'_summary.txt';

open (AO, '>>', $allout);
open (TO, '>>', $topout);
open (SO, '>>', $sumout);

#print headers
print AO "Protein query\tPhylum\tSpecies\tTranscriptome\tsseqid\tpident\tlength\tevalue\tqcovs\tqcovhsp";
print TO "Protein query\tPhylum\tSpecies\tTranscriptome\tsseqid\tpident\tlength\tevalue\tqcovs\tqcovhsp";

#Read in the fasta file of queries for blasting and write the 'all' and 'top' files. 

####
my $debug = 'Debug.txt';
open (DEBUG, '>>', $debug);

{
local $/ = ">";
my $protcount = '0';#first record is simply '>', so $protcount must start at zero.
while (my $record = <FAIN>){
	my %SUMMDIAT; #Record information for the 'summary' file in the two hashes below
	my %SUMMOTH;
	if ($record =~ /^([^\n]*)\n(.*)/gs ){
		my $name = $1;		
		my $seq = $2;		
		$seq =~ s/>//g;
		
		#write temp blast input
		my $blastin = 'tempblastin.fasta';
		open (BLIN, '>>', $blastin);
		print BLIN ">$name\n$seq";
		close BLIN;
		
		foreach my $phylum (sort keys %TAXON){
			
		
			#For diatoms perform blast routine
		
			if ($phylum eq 'Bacillariophyta'){
				####
				print DEBUG 'Bacillariophyta\n';
				foreach my $species (sort keys %{ $TAXON{$phylum} } ){
					
					#species level
					my $hitcount = '0'; #counts all hits for the species, accross all transcriptomes
					foreach my $tomes (sort keys %{ $TAXON{$phylum}{$species} } ){
						#transcriptome level
						#my $mmetsp = $tomes.'.trinity_out_2.2.0.Trinity.fasta';
						system (qq(tblastn -db $tomes -query tempblastin.fasta -evalue $evalue -outfmt "6 qseqid sseqid pident length evalue qcovs qcovhsp" -num_threads 6 -out tempblastout.txt));
						my $blastout = 'tempblastout.txt';
						open (BLO, '<', $blastout);
						my $previous = '0';
						my $top = '0'; #for selecting top hit for the transcriptome
						{
						local $/ = "\n";
							while (my $line = <BLO>){
								if ($line =~ /^[^\t]*\t(([^\t]*)\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)\t[^\t]*)/){
									my $outline = $1;
									my $transcript = $2;
									my $qcov = $3;
									if ($qcov >= $qcovco){ #qcov threshold implemented
										if ($transcript ne $previous){
											if ($top == '0'){										
												print TO "\n$name\t$phylum\t$species\t$tomes\t$outline";
												print AO "\n$name\t$phylum\t$species\t$tomes\t$outline";
												$hitcount++;
												$top++;
											}else{
												print AO "\n$name\t$phylum\t$species\t$tomes\t$outline";
												$hitcount++;
											}
										}							
										$previous = $transcript;
									}
								}
								
							}
						}
						close BLO;
						system ('rm tempblastout.txt');
					}
					$SUMMDIAT{$species} = $hitcount;
				}
			}	
			#For other taxa, perform blast routine
			
			
			if ($phylum ne 'Bacillariophyta'){
				my $hitcount = '0';
				print DEBUG "\n\n\nNOT Bacillariophyta\n hitcount is $hitcount\n";
				print DEBUG "\t\t\tPHYLUM is $phylum\n";
				foreach my $species (sort keys %{ $TAXON{$phylum} } ){
					print DEBUG "\nSpecies is $species ****\n";
					foreach my $tomes (sort keys %{ $TAXON{$phylum}{$species} }){
						print DEBUG "\n\tTranscriptome is $tomes ****\n";
						#my $mmetsp = $tomes.'.trinity_out_2.2.0.Trinity.fasta';
							system (qq(tblastn -db $tomes -query tempblastin.fasta -evalue $evalue -outfmt "6 qseqid sseqid pident length evalue qcovs qcovhsp" -num_threads 6 -out tempblastout.txt));
							my $blastout = 'tempblastout.txt';
							open (BLO, '<', $blastout);
							my $previous = '0';
							my $top = '0'; #for selecting top hit for the transcriptome
							{
							local $/ = "\n";
								while (my $line = <BLO>){
									if ($line =~ /^[^\t]*\t(([^\t]*)\t[^\t]*\t[^\t]*\t[^\t]*\t([^\t]*)\t[^\t]*)/){
										my $outline = $1;
										my $transcript = $2;
										my $qcov = $3;
										if ($qcov >= $qcovco){ #qcov threshold implemented
											if ($transcript ne $previous){ #ensures only works with top HSP
												if ($top == '0'){										
													print TO "\n\t\t$name\t$phylum\t$species\t$tomes\t$outline";
													print AO "\n\t\t$name\t$phylum\t$species\t$tomes\t$outline";
													print DEBUG "\n\t\t$name\t$phylum\t$species\t$tomes\t$outline";
													$hitcount++;
													print DEBUG "\n\t\t hitcount from top =0 loop is: $hitcount";
													$top++;
												}else{
													print AO "\n\t\t\t$name\t$phylum\t$species\t$tomes\t$outline";
													print DEBUG "\n\t\t\t$name\t$phylum\t$species\t$tomes\t$outline";
													$hitcount++;
													print DEBUG "\n\t\t\t hitcount from top > 0 loop is: $hitcount";
												}
											}							
											$previous = $transcript;
										}
									}
									
								}
							}
							close BLO;
							system ('rm tempblastout.txt');
						}
						$SUMMOTH{$phylum}{$species} = $hitcount;
						print DEBUG "\n\tWRITING SUMMOTH hash: SUMMOTH : $phylum : $species = $hitcount\n";
						$hitcount = '0'; ###
				}
			}
		}
		system ('rm tempblastin.fasta');
		
		#Now generate the summary file from the count information stored in hashes
		print "protcount is $protcount\n";
		#if the first protein, generate header and first line (transcriptome counts):
		if ($protcount == '1'){
			print SO "Protein query";
			foreach my $species (sort keys %SUMMDIAT){
				print SO "\t$species";
			}
			foreach my $phyla (sort keys %SUMMOTH){
				print SO "\t$phyla";
			}
			print SO "\nnumber of transcriptomes";
			foreach my $species (sort keys %SUMMDIAT){
				print SO "\t$TOMEC{'Bacillariophyta'}{$species}";
			}
			foreach my $phyla (sort keys %SUMMOTH){
				my $phylatranscount = '0';
				foreach my $species (sort keys %{ $SUMMOTH{$phyla} } ){
					my $transcount = $TOMEC{$phyla}{$species};
					$phylatranscount += $transcount;
				}
				print SO "\t$phylatranscount";
			}
		}
		
		#Print summary info for this protein
		
		print SO "\n$name";
		print DEBUG "\n\n @@@@@@@@ Generating output from SUM hash @@@@@@@@@\n\n";
		foreach my $species (sort keys %SUMMDIAT){
			my $hits = $SUMMDIAT{$species};
			#my $transcount = $TOMEC{'Bacillariophyta'}{$species};
			#print "species really is $species\n";
			# $normhits = $hits.';'.$transcount;
			#print "\tnormhits is $normhits\n";
			print SO "\t$hits";
		}
		
		foreach my $phyla (sort keys %SUMMOTH){
			#print "running SUMMOTH loop\n";
			#print "\tphyla is $phyla\n";
			my $phylahits = '0';
			#my $phylatranscount = '0';
			print DEBUG "phyla is $phyla\n";
			foreach my $species (sort keys %{ $SUMMOTH{$phyla} } ){
				#print "\tspecies is $species\n";
				print DEBUG "\tspecies is $species \n";
				my $hits = $SUMMOTH{$phyla}{$species};
				print DEBUG "\t\tnumber of hits is: $hits\n";
				$phylahits += $hits;
				print DEBUG "\t\t\tValue of phylahits is $phylahits\n";
				#my $transcount = $TOMEC{$phyla}{$species};
				#$phylatranscount += $transcount;
			}
			#my $normhits = $phylahits.';'.$phylatranscount;
			print SO "\t$phylahits";
			print DEBUG"\tEnd value of phylahits is $phylahits\n";
		}
			
		
	}
	$protcount++;
}		
}

close FAIN;
close MMIN;
close TAIN;
close PYIN;
close AO;
close TO;
close SO;
close DEBUG;

my $endtime = time();
my $runtime =  $endtime - $^T;
print "run time was $runtime\n";
	