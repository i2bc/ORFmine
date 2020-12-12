# How works ORFold?

ORFold is a tool developed in python3 which aims at 
characterizing the fold potential of a set of amino acid sequences with no 
knowledge of their 3D structures. The fold potential of a sequence is 
calculated with the HCA 
method [REF]. Also ORFold can estimate the disorder, 
and the aggregation propensities of the input sequences with IUPred
and Tango respectively.    


## Hydrophobic Clusters Analysis (HCA)
**HCA**[REF] aims as delineating in an amino acid sequence, regions enriched in 
strong hydrophobic residues (HCA clusters) and regions 
of at least four consecutive non-hydrophobic residues (HCA linkers). 
The patterns of hydrophobic residues can be associated with specific regular 
secondary structures, and the distribution of the HCA clusters and linkers in a protein 
sequence can be used to estimate through the HCA score, its ability to fold (completely or partially). 

This score ranges from -10 to +10 with low HCA scores indicating
sequences depleted in hydrophobic clusters and expected to be disordered in solution, 
while high HCA scores reflect sequences enriched in hydrophobic clusters 
and expected to generate aggregates in solution, though some of them could
fold in lipidic environments. 
Foldable sequences are known to display
an equilibrium between hydrophobic and hydrophilic residues (average of 33% 
of hydrophobic residues in globular proteins). 
They are mostly associated with intermediate HCA score values.
 

The HCA score is calculated using the freely available 
software **pyHCA** which can be downloaded and installed 
following the instructions of its developers: <https://github.com/T-B-F/pyHCA>


## Tango
**Tango**[REF] is a method which aims at predicting aggregation nucleating regions
in protein sequences. 
If specified by the user, ORFold can calculate and add the aggregation propensity 
of a sequence in the output. 
Tango is not freely available software, and the user of ORFold should 
first contact the Tango developers to have access to the source code: <http://tango.crg.es>

For the aggregation propensity estimation, according to the protocol
proposed by XXX et al[REF], a residue is considered as
participating in an aggregation prone region if it is located in a segment 
of at least five consecutive residues which were predicted as populating 
a b-aggregated conformation for more than 5%. 
Then, the aggregation propensity of each sequence is defined as the 
fraction of residues predicted in aggregation prone segments. 

## IUPred
**IUPred2A**[REF] is one of the best methods for the prediction of 
Intrinsically Disordered Proteins (IDPs) and can be used as a 
complement to the HCA score prediction. 
If specified by the user, ORFold can calculate and add the disorder propensity 
of a sequence in the output. 
IUPred is not freely available, and the user of 
ORFold should first contact the IUPred developers to 
have access to the source code : <https://iupred2a.elte.hu>

For the disorder propensity estimation, in order to be consistent with
the estimation of the aggregation propensity, ORFold searches for 
regions on the protein sequence that present at least five consecutive 
residues with a disorder probability higher than 0.5. 
The disorder propensity of each sequence is defined as the fraction 
of residues predicted as located in a highly disordered segment.    
