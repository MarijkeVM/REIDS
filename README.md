# REIDS
The software repository for the paper

1) PreProcessing\n
This folder contains the .CEL Files and .cdf files necessary to preprocess the HTA data.
The result is an .RData file containing the measured probe intensities for every sample.
Further, exon and gene level summarized values can be obtained.
	
Note: the folder structure is important for the aroma.affymetrix package to work.
	
	
2) REIDS\n
This folder contains the REIDS software to either run it on your computer (not recommended) or on a supercomputer infrastructure.
After running the REIDS model, the junction information is used to rank the probe sets identified as alternatively spliced.
	
	
3) Analysis\n
The analysis folder contains the scripts used to obtain the content of the paper.
