Running the items in parallel:

Convert the .RData file after running HTA_PreProcessing.R to a file in which every item (gene) is represented by a single line:
1) run PivotTransformation.R on the R.Data file --> .csv file 
2) run LineIndexer.py on the .csv file --> .csv file
The latter .csv file contains the line indices and the corresponding start and end points.


