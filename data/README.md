# Data

This directory contains the base data files that are used as input to the
analyses in this project.

## model_info.csv

Info about the organisms that are used in the analyses.

## biomass_db.tsv

This file was taken from [CarveMe](https://github.com/cdanielmachado/carveme/blob/master/carveme/data/input/biomass_db.tsv).
The composition of biomass reactions in BiGG template models were based on this table.

## substrate_utilization_data.csv

The data in this folder was pulled from multiple sources, outlined in the table
below:

| Organism        | Reference    | Notes                                   |
| --------------- | :----------: | --------------------------------------- |
| *E. coli*       | [1](iaf1260) |                                         |
| *B. subtilis*   | [2](iyo844)  |                                         |
| *P. putida*     | [3](pputi)   |                                         |
| *C. difficile*  | [4](cdiff)   | Growth > 1.4 (from paper)               |
| *S. aureus*     | [5](iys854)  |                                         |
| *A. baumannii*  | [6](abaum)   | Growth > 101 OmniLog units (from paper) |

## mediadb.tsv

This file was created to mimic the media DB structure needed for CarveMe. It
specifies the compositions of all the media needed for gap-filling and running
substrate utilization predictions with this dataset.

<!-- Links -->
[iaf1260]: https://link.springer.com/article/10.1038/msb4100155
[iyo844]: https://www.sciencedirect.com/science/article/pii/S0021925820586110
[pputi]: https://journals.asm.org/doi/10.1128/msystems.01594-25
[cdiff]: https://www.sciencedirect.com/science/article/pii/S1201971214015781
[iys854]: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006644
[abaum]: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0058628
