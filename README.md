# PhyloBlast_diatoms
Tools for phylogenetic analysis using the MMETSP transcriptomes

At the moment this tool is a set of scripts for blasting against the diatom transcriptomes within the MMETSP databases. Basic bioinformatic skills are required to run it. I'm working on a Dockerized version incorporating all databases which will work 'out of the box' and be available as an App with a GUI on the Cyverse Discovery Environment. I'm also working on a more generalized version to allow focus on any taxonomic range desired and that uses plast instead of blast to make it possible to run on a standard desktop in a reasonable time frame. 

If you're interested in beta testing these tools then please let me know!

Assuming you have access to a Linux machine or server, this is how to proceed.

1. Download the version of the MMETSP transcriptomes that have been cleaned of contaminants.
2. Make blast databases for all transcriptomes
3. Set BLASTDB and PATH variables to the location of the databases
4. Run the scripts as described below:

The first script runs the blast searches and collates the results. The best hit per transcriptome is provided and a phylogenetic summary of the results.

```{}
perl mmetspblast3.pl [input fasta] [input taxa: MMETSP] [input taxa: taxonomy] [input: phyla name to number] [evalue threshold] [qcov cutoff] [output prefix]
```
Where the input files are as follows:

Input fasta: The sequences you want to blast
Input taxa MMETSP: The file '

The number of cores in hard coded in the blast command in the main script. Please alter as neccessary.


```{}
	
perl -S mmetsp_evoprot.pl [top hits out] [fasta of queries] [output prefix]

```



