# Batch Correction for ENCORE data:

### Aim:
Correct raw logFC across COLO and BRCA libraries (independentily). 

### How:
Use ComBat strategy for guide pairs in common and approximate the parameters for pair guides not in common via kNN.

### Input:
- Table of raw logFCs guide pairs x CLs
- Libraries annotation

### Output:
- Plots for validation of the strategy and parameter selection
- Table of batch corrected logFCs guide pairs x CLs
- Sanity check plots
