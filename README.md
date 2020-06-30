# PhyloBlast_diatoms
Tools for phylogenetic analysis using the MMETSP transcriptomes

At the moment this tool is a set of scripts for blasting against the diatom transcriptomes within the MMETSP databases. Basic bioinformatic skills are required to run it. I'm working on a Dockerized version incorporating all databases which will work 'out of the box' and be available as an App with a GUI on the Cyverse Discovery Environment. I'm also working on a more generalized version to allow focus on any taxonomic range desired and that uses plast instead of blast to make it possible to run on a standard desktop in a reasonable time frame. 

If you're interested in beta testing these tools then please let me know!

Assuming you have access to a Linux machine or server, this is how to proceed.

1. Download the version of the MMETSP transcriptomes that have been cleaned of contaminants.
2. Make blast databases for all transcriptomes
3. Set BLASTDB and PATH variables to the location of the databases
4. Run the scripts as follows described below. 





