# NTE
Network Reconstruction by Transfer Entropy. For license information, see 
[LICENSE](./LICENSE).

## Overview

NTE uses transfer entropy to estimate the an edge list based on expression values for different sets of genes that cary in time.

INPUT: An expression matrix that has valiues of expression at different time points. The matrix structure should be configured so that genes are represented in columns and time points in rows. It should include a header with the names of the genes. The time points sould be in increasing order.

OUTPUT: An edge list of the interactions of the matrix. The edge list will include all the pairs of genes (self interactions are ignored. And a values of significance for each pair.
##System requirements

A working version of R with the following packages [https://www.r-project.org/](https://www.r-project.org/).

	optparse
	entropy
	gdata

## Installation

Clone the git repository

   ```bash
   $> git clone https://github.com/jccastrog/NTE
   ``` 

You can also download the zip file from the GitHub site [https://github.com/jccastrog/NTE](https://github.com/jccastrog/NTE).

## Running NTE

### Input

The input is a table that includes the time series expression data for a series of transcription factors, columns are genes and rows are time points For each time point replicates have to be consecutive if replicates are missing NA should be added.

Example:

	Gene1	Gene2
	Time1_Rep1 Time1_Rep1
	Time1_Rep2 Time2_Rep1
	Time2_Rep1 Time2_Rep1
	Time2_Rep2 Time2_Rep1

### Output

The output is a table with tab separated values, the first column is the regulator node and the second column is the regulated gene. Third column is the mutual information value for the two expression profiles, without lag. And the fourth column is the pValue for the interaction.

Example:

	gene1	gene2	MI	pVal
	Gene1	Gene2	0.284118818391665	0.007
	Gene2	Gene1	0.284118818391665	0.304

### Options

	num_reps[INTEGER]: The number of replicates of the experiment.
	num_rand[INTEGER]: Number of mutual information randomizations. Used to calculate the p value of mutual information.
	dynamical[LOGICAL]: Wheter or not the optimal lag step is chosen for each gene. Might increase the computation time significantly.
	max_step[INTEGER]: Maximum step to use as lag.
	step_size[INTEGER]: The step size to be used by default for all genes. Invalid if dynamical=TRUE.
