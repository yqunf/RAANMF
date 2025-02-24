# RAANMF
The feature representation of sequence based on non-negative matrix factorization and its application in protein thermostability prediction
# Document introduction
> 1.dataset: The Charoenkan.fasta file is used for training the GAA-NMF model, while the Lin.fasta file is divided into a training set and a test set to validate the performance of amino acid grouping.  
> 2.GAA-NMF.py: The GAA-NMF model  
> 3.readfasta: Read the fasta file.  
> 4.Gm-AAC: m groups-amino acid composition.  
> 5.Gm-DPC: m groups-dipeptide composition.  
> 6.Gm-TPC: m groups-tripetide composition.
# 
| Group name | Amino acid grouping |
|------------| -------------------|
| GAA-HP(5) |	**G**AVLMI-**F**YW-**K**RH-**D**E-**S**TCPNQ |
| GAA-PC(7) |	**A**GV-**C**-**D**E-**F**ILP-**H**NQW-**K**R-**M**STY |
| GAA-PAM(6) | **A**GPST-**D**ENQ-**H**RK-**I**LMV-**F**WY-**C** |
| GAA-CS(7) |	**R**EQKACML-**H**WYF-**T**IV-**N**D-**P**-**S**-**G** |
| GAA-KM(8) |	**D**E-**Q**ST-**L**MI-**V**AF-**N**H-**Y**PRK-**G**W-**C**|
| GAA-NMF(8) |**A**-**C**DHNQST-**E**-**F**IMY-**G**PV-**L**W-**R**-**K** |
