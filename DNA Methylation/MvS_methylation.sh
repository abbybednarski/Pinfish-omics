#!/bin/bash

#load 4.5.1 version that has methylKit and genomation
source /apps/profiles/modules_asax.sh.dyn
module load R/4.5.1

#remove DMR output to start from scratch every run
rm /home/aubaxb001/methylation_data/trimmed/DMR_prawn_MvS.Rout

#command to run R
R CMD BATCH --vanilla /home/aubaxb001/DMR_prawn_MvS.R /home/aubaxb001/methylation_data/trimmed/DMR_prawn_MvS.Rout
