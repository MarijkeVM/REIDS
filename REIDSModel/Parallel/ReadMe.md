# Running the items in parallel:

* Convert the .RData file after running HTA_PreProcessing.R to a file in which every item (gene) is represented by a single line:
  + run PivotTransformation.R on the R.Data file --> .csv file (batch script is PivotTransformation.pbs)
  + run Line_Indexer.py on the .csv file --> .csv file (batch script is line_indexer.pbs)
  
The latter .csv file contains the line indices and the corresponding start and end points.
Both files are necessary in order for REIDS.R to run on a supercomputer infrastructure.
Further, a geneID.csv has been created and will be used later when concatenating the results.


* Run REIDS_Cluster.R on the items of the Line_Indexer.py. (batch script is REIDS.pbs).
  + The computation time is approximately a day with the provided settings.
  + Review the specified folders for saving and alter to your settings.


* Concatenate the results with CreateOutput.R (batch script is CreateOutput.pbs)


