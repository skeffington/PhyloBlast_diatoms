
#mmetsp_evoprot.pl
#requires EMBOSS.
use warnings;
use strict;
#use Cwd;

my $usage = "perl -S mmetsp_evoprot.pl [top hits out] [fasta of queries] [output prefix]\n";
#my $dir = getcwd();
my $tophit = shift or die $usage;
my $fastain = shift or die $usage;
my $outpre = shift or die $usage;

#Read in top hits file and store information in %TOP hash:
open(TH, '<', $tophit);

my %TOP;
my $index = '0';

while (my $line = <TH>){
	$line =~ s/\R//;
	if ($line =~ /Protein query/){
		next;
	}
	if ($line =~ /^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t[^\t]*\t([^\t]*)\t([^\t]*)\t[^\t]*/){
		my $prot = $1;
		my $phylum = $2;
		my $species = $3;
		
		my $transcriptome = $4;
		my $transcript = $5;
		#print "transcript is $transcript\n";
		my $pcid = $6;
		my $eval = $7;
		my $qcov = $8;
		
		$TOP{$prot}{$transcriptome}{'phylum'} = $phylum ;
		$TOP{$prot}{$transcriptome}{'species'} = $species;
		$TOP{$prot}{$transcriptome}{'transcript'} = $transcript;
		$TOP{$prot}{$transcriptome}{'pcid'} = $pcid;
		$TOP{$prot}{$transcriptome}{'eval'} = $eval;
		$TOP{$prot}{$transcriptome}{'qcov'} = $qcov;
		$TOP{$prot}{$transcriptome}{'index'} = $index;
		$index++;
		#print "index is $TOP{$prot}{$transcriptome}{'index'} ; transcript is  $TOP{$prot}{$transcriptome}{'transcript'}\n";
	}
}

###Database of groups hardcoded. Note not all MMETSPs are in these groups. Some are too insignificant for inclusion.

#The names stored in the hash are the 'phylum' column of the tophits file for non-diatoms, but the 'species' column of the top hits file for diatoms.

my %GROUPS;

$GROUPS{'notdia'}{'Miozoa'} = ['Apicomplexa', 'Dinophyceae']; 
$GROUPS{'notdia'}{'Rhodophyta'} = ['Bangiophyceae', 'Compsopogonophyceae', 'Rhodellophyceae', 'Stylonematophyceae'];
$GROUPS{'notdia'}{'Ochrophyta'} = ['Bolidophyceae', 'Chrysophyceae', 'Dictyochophyceae', 'Pelagophyceae', 'Pinguiophyceae', 'Raphidophyceae', 'Synchromophyceae', 'Synurophyceae'];
$GROUPS{'notdia'}{'Cercozoa'} = ['Cercozoa'];
$GROUPS{'notdia'}{'Chlorophyta'} = ['Chlorophyta'];
$GROUPS{'notdia'}{'Choanoflagellida'} = ['Choanoflagellida'];
$GROUPS{'notdia'}{'Ciliophora'} = ['Ciliophora'];
$GROUPS{'notdia'}{'Haptophyta'} = ['Coccolithales', 'Coccosphaerales', 'Haptophyta incertae sedis', 'Isochrysidales', 'Pavlovales', 'Phaeocystales', 'Prymnesiales'];
$GROUPS{'notdia'}{'Cryptophyta'} = ['Cryptomonadales'];
$GROUPS{'notdia'}{'Glaucophyta'} = ['Cyanoptyche', 'Gloeochaetales'];
#Diatoms. Should perhaps exclude T. oceanica??
#Note Bacillariophyceae includes Bacillariophyta classis incertae sedis
$GROUPS{'diatom'}{'Bacillariophyceae'} = ['Cylindrotheca closterium','Fragilariopsis kerguelensis','Nitzschia','Pseudo-nitzschia arenysensis','Pseudo-nitzschia australis','Pseudo-nitzschia delicatissima','Pseudo-nitzschia fraudulenta','Pseudo-nitzschia heimii','Pseudo-nitzschia pungens','Tryblionella compressa','Cyclophora tenuis','Staurosira','Synedropsis cf. recta','Licmophora paradoxa','Craspedostauros australis','Amphiprora','Amphiprora paludosa','Stauroneis constricta','Grammatophora oceanica','Asterionellopsis glacialis','Striatella unipunctata','Entomoneis sp.','Thalassionema frauenfeldii','Thalassionema nitzschioides','Triceratium dubium','Amphora coffeiformis','Astrosyne radiata'];
$GROUPS{'diatom'}{'Coscinodiscophyceae'} = ['Aulacoseira subarctica','Corethron hystrix','Corethron pennatum','Coscinodiscus wailesii','Dactyliosolen fragilissimus','Proboscia alata','Proboscia inermis','Rhizosolenia setigera','Stephanopyxis turris'];
$GROUPS{'diatom'}{'Mediophyceae_non_Thalassiosirales'} = ['Attheya septentrionalis','Eucampia antarctica','Helicotheca tamesis','Chaetoceros affinis','Chaetoceros brevis','Chaetoceros cf. neogracilis','Chaetoceros curvisetus','Chaetoceros debilis','Chaetoceros dichaeta','Chaetoceros neogracilis','Chaetoceros sp.','Leptocylindrus aporus','Leptocylindrus danicus','Extubocellulus spinifer','Minutocellus polymorphus','Odontella aurita','Trieres chinensis','Ditylum brightwellii','Cyclotella meneghiniana','Thalassiothrix antarctica'];
$GROUPS{'diatom'}{'Skeletonema'} = ['Skeletonema costatum','Skeletonema dohrnii','Skeletonema grethae','Skeletonema japonicum','Skeletonema marinoi','Skeletonema menzellii'];
$GROUPS{'diatom'}{'Detonula'} = ['Detonula confervacea'];
$GROUPS{'diatom'}{'Thalassiosira'} = ['Thalassiosira','Thalassiosira antarctica','Thalassiosira gravida','Thalassiosira minuscula','Thalassiosira oceanica','Thalassiosira punctigera','Thalassiosira rotula','Thalassiosira weissflogii'];

###



#Add group information to the TOP hash. Compare hits and find the best hit for each group for each protein - collect the indicies of these records in a hash.

my %GRHITS; #stores a list of index values for records which are the best hit in the group.

foreach my $prot (keys %TOP){
	#print "prot is $prot\n";
	my %GHIT; #stores the hits for each of the groups for comparison for a given protein
	#print "###############\n";
	foreach my $ttome (keys %{ $TOP{$prot} } ){
			my $phylum = $TOP{$prot}{$ttome}{'phylum'};
			my $species = $TOP{$prot}{$ttome}{'species'};
			my $index = $TOP{$prot}{$ttome}{'index'};
			my $eval = $TOP{$prot}{$ttome}{'eval'};
			#print "transcriptome is $ttome ; phylum is $phylum; species is $species; index is $index; eval is $eval\n";
			#non-diatoms
			
			if ($phylum ne 'Bacillariophyta'){
				foreach my $group ( keys %{ $GROUPS{'notdia'} } ){
					#print "array is @{ $GROUPS{'notdia'}{$group} }\n";
					if( $phylum ~~ @{ $GROUPS{'notdia'}{$group} } ){
						$TOP{$prot}{$ttome}{'group'} = $group;
						#print "\tphylum match is $phylum\n";
						$GHIT{$group}{$index} = $eval;
						#print "GIT: group is $group ; index is $index ; eval is $eval\n";
					}
				}
			}
			#Diatoms
			if ($phylum eq 'Bacillariophyta'){
				foreach my $group ( keys %{ $GROUPS{'diatom'} } ){
					if( $species ~~ @{ $GROUPS{'diatom'}{$group} } ){
						$TOP{$prot}{$ttome}{'group'} = $group; 
						$GHIT{$group}{$index} = $eval;
						#print "Diatom GIT: group is $group ; index is $index ; eval is $eval\n";
					
					}
				}			
			}
	}
	#print "###############\n";
	#Compare hits for each protein and group and choose the best.
	foreach my $group ( keys %GHIT ){
		#print "\t group is $group\n";
		foreach my $index ( keys %{ $GHIT{$group} } ){
			#print "index is $index; eval is $GHIT{$group}{$index}\n";
		}
		my @sortindex = sort { $GHIT{$group}{$a} <=> $GHIT{$group}{$b} } keys %{$GHIT{$group}};
		my $lowindex = shift @sortindex;
		#print "lowindex is $lowindex\n";
		$GRHITS{$lowindex}=();
	}
	undef %GHIT;
}

#foreach my $index (keys %GRHITS){ # problem = not every group has a lowindex value for each protein...
#	print "$index\n";
#}

####Generate outputs:

my $evalout = $outpre . '_eval.txt';
my $detailout = $outpre . '_details.txt';
open(EOUT, '>>', $evalout);
open(DOUT, '>>', $detailout);

#print headers
print EOUT "protein\teval\tqcov\tpcid\tgroup";
#my $stopvar = '0';
my $count = '0';

foreach my $prot (keys %TOP){
	print DOUT "###########################################\n$prot\n###########################################\n";
	my $temp1 = 'prot_tmp1.fa';
	system(qq(awk 'BEGIN{RS=">";FS="\\n"}NR>1{if (\$1~/$prot/) print ">"\$0}' $fastain >  $temp1));
			foreach my $ttome (keys %{ $TOP{$prot} } ){
				my $index = $TOP{$prot}{$ttome}{'index'};
				if ( exists $GRHITS{$index} ){ # This is one of the best in group hits. 
						#print "###########index found in list: $index\n";
					
						my $eval = $TOP{$prot}{$ttome}{'eval'};
						my $qcov = $TOP{$prot}{$ttome}{'qcov'};
						my $pcid = $TOP{$prot}{$ttome}{'pcid'};
						my $group = $TOP{$prot}{$ttome}{'group'};
						print EOUT "\n$prot\t$eval\t$qcov\t$pcid\t$group";
						print DOUT "$group\n\n";
						my $transcript = $TOP{$prot}{$ttome}{'transcript'};
						$transcript =~ s/\|/\\|/g; #necessary for awk search to function properly, because of pipes in name.
						$transcript = $transcript . '$';
						#print "transcript line 156 is $transcript\n";
						my $mmetspdb = "/home/blast/MMETSP/MMETSP_2019_2/$ttome";
						#print "mmetspdb line 159 is $mmetspdb\n";
						my $temp2 = 'trans_tmp2.fa';
						system(qq(awk 'BEGIN{RS=">";FS="\\n"}NR>1{if (\$1~/$transcript/) print ">"\$0}' $mmetspdb >  $temp2));
						#die;
						#change name of sixpack input so that transcript name is included.
						my $spin = 'spin_temp.txt';
						open (SPIN, '>>', $spin);
						open(TMP2, '<', $temp2);
						my $short_tran_id = 'na';
						
						{
						local $/ = ">";
							my $short_tran_id = 'na';
							while (my $rec = <TMP2>){ #print the correct translation of the transcript
								if ($rec =~ /^([^\n]*)\n(.*)/gs ){
									my $name = $1;		
									my $seq = $2;		
									$seq =~ s/>//g;
									if ($name =~ /.*(MMETSP[^\|]*)\|(\d*)/){
										my $mmet = $1;
										my $trid = $2;
										$short_tran_id = "$mmet" . '-' . "$trid";
									}
									print SPIN ">$short_tran_id\n$seq";
								}
							}
						}
					
						close TMP2;
						close SPIN;
						
						#run sixpack and blasp the result against the query protein
						system(qq(sixpack -sequence $spin -outfile temp3.txt -outseq sp_temp4.txt));
						#sleep(5);
						
						my $temp5 = 'smtemp5' . "$count";
						system(qq(blastp -query $temp1 -subject sp_temp4.txt -out $temp5));  ##doesn't really work for poor alignments - need blast...
						
						#select best matching translation from blasp output:
						
						open(SMA, '<', $temp5);
						my $topmatch;
						while (my $line = <SMA>){
							$line =~ s/\R//;
							if ($line =~ /^\s*(MMETSP[^\s]*)/){
								$topmatch = $1; #note this is short name form eg MMETSP1059-20121228-47836_3_ORF1
								last;
							}
						}
						
						#now generate detailed output for manual curation of the results
						open(TMP2, '<', $temp2);
						while (my $line = <TMP2>){ #print the transcript sequence of the best match for the group
							$line =~ s/\R//g;
							print DOUT "$line\n";
						}
						close TMP2;

						my $sixpack = 'sp_temp4.txt';
						open(SPAC, '<', $sixpack);
						my $translated = 'temptranl.fasta';
						open (TRANS, '>>', $translated);
						#print "topmatch is $topmatch\n";
						{
						local $/ = ">";
							while (my $rec = <SPAC>){ #print the correct translation of the transcript
								if ($rec =~ /^([^\n]*)\n(.*)/gs ){
									my $name = $1;
									#print "\nname is $name\n";
									my $seq = $2;
									#print "seq is \n$seq\n";
									$seq =~ s/>//g;
									if ($name =~ /$topmatch\s/){
										#print "print outputs\n";
										print DOUT ">$name\n$seq\n";
										print TRANS ">$name\n$seq\n";
									}
								}
							}
						}
						close SPAC;
						close TRANS;
						
						#run blastp for the translated transcript vs the query protein and print to output
						my $blastout = 'tempblastout.txt';
						#sleep(1);
						system(qq(blastp -query $temp1 -subject $translated -out tempblastout.txt));
						#sleep(1);
						open (BLAOUT, '<', $blastout);
						my $switch = '0';
						print DOUT "Blastp output:\n";
						while (my $line = <BLAOUT>){
							if ($switch == 1){
								print DOUT "$line\n";
							}
							if ($line =~ /Query=/){
								print DOUT "$line\n";
								$switch = 1;
							}
							if ($line =~ /Effective/){
								$switch = 0;
							}
						}
						close BLAOUT;
								
						#remove temp files
						$count++;
						#$stopvar++;
						#if($stopvar == 4){ 
						#		die;
						#}
						#if($group eq 'Chlorophyta'){die};
						system(qq(rm tempblastout.txt temptranl.fasta spin_temp.txt sp_temp4.txt $temp5 temp3.txt trans_tmp2.fa));
						
					
				}
			}
	system('rm prot_tmp1.fa');
}

close EOUT;
close DOUT;
close TH;

#END
