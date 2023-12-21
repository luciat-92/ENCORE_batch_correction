# Batch Correction for ENCORE data:
## PROJECT ARCHIVED! This step was integrated in the overall pipeline
### Aim:
Correct raw logFC across COLO and BRCA libraries (independentily). 

### How:
Use ComBat strategy for guide pairs in common and approximate the parameters for pair guides not in common via kNN.

### Input:
- Table of raw logFCs guide pairs x CLs
- Libraries annotation

### Output:
- Sanity check plots
- Plots for validation of kNN strategy and parameter selection
- Data Frame of logFCs guide pairs x CLs plus guide pairs info

