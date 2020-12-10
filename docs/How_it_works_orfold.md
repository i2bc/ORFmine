# How works ORFold?

ORFold is a tool developed in python3 which aims at characterizing the foldability potential of any amino acid sequence given. It is principally based on the HCA method by calculating the HCA foldability score while it can also be complemented with the prediction of the disorder propensity (by IUPred) and the aggregation propensity (by Tango).   


## Hydrophobic Clusters Analysis (HCA)
HCA is a method based on the detection of regions enriched in strong hydrophobic residues (hydrophobic clusters) whose patterns are associated to specific regular secondary structures. The distribution of these hydrophobic patterns on the amino acid sequence can be an indication for its ability to fold (complitely or partially) in a more or less well defined 3D conformation. 

HCA estimates the enrichment of a sequence in hydrophobic patters and transforms it into a score indicative of its foldability. This score ranges from -10, for sequences depleated in hydrophobic clustes and expected to be disordered in solution, to +10, for sequences enriched in hydrophobic clusters and expected to generate aggregates in solution. Sequences that can fold uppon the hydrophobic effect present an equilibrium between hydrophobic and hydrophilic resdiues and adopt mostly intermediate HCA score values.

One should notice that the foldability potential of any sequence can also contain the notion of complete disorder (lack of foldabilty) which makes the HCA score    

The HCA score is calculated using the freely available software **pyHCA** which can be downloaded and installed following the instructions of its developers: <https://github.com/T-B-F/pyHCA>

## IUPred
**IUPred2A** is one of the best methods for the prediction of Intrinsically Disordered Proteins (IDPs) and can be used as a complement or a validation for the HCA score prediction. ORFold can calculate and integrate in its output the disorder propensity of a sequence (but this stays optional). However IUPred is not a freely available software and the user of ORFold should first contact the IUPred developers and have access to the source code : <https://iupred2a.elte.hu>

For the disorder propensity estimation, ORFold aims at detecting regions on the protein sequence that present at least five consecutive amino acids with disorder probability more than 0.5. The disorder propensity of each sequence  is defined as the fraction of residues predicted located in a highly disordered segment.    

## Tango
**Tango** is a method which aims at predicting aggregation nucleating regions in protein sequences. ORFold can calculate and integrate the aggregation propensity of a sequence (but this stays optional). However Tango is not a freely available software and the user of ORFold should first contact the Tango developers and have access to the source code: <http://tango.crg.es>

For the aggregation propensity estimation, a residue was considered as participating in an aggregation prone region if it was located in a segment of at least five consecutive residues which were predicted as populating a b-aggregated conformation for more than 5%. Then, the aggregation propensity of each sequence is defined as the fraction of residues predicted in aggregation prone segments. 
