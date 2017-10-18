Running the items in parallel:

*
Convert the .RData file after running HTA_PreProcessing.R to a file in which every item (gene) is represented by a single line:
1) run PivotTransformation.R on the R.Data file --> .csv file (batch script is PivotTransformation.pbs)
2) run Line_Indexer.py on the .csv file --> .csv file (batch script is line_indexer.pbs)
The latter .csv file contains the line indices and the corresponding start and end points.
Further, a geneID.csv has been created and will be used later when concatenating the results.


*
Run REIDS_Cluster.R on the items of the Line_Indexer.py.
The computation time is approximately a day with the provided settings.


