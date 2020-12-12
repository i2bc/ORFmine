## Running ORFold:


### Inputs
For its most basic run, ORFold requires only a FASTA file containing the amino acid sequences 
to treat (given with the **-fna** label). ORFold can handle several FASTA files at the same
time. In this case, it will treat them independently and will generate as many 
outputs as entered FASTA files.

 <p> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
 FASTA file example:
```{}
>aminoacid_sequence_1
AGNVCFGGRTYMPDFDGMSCVNWQERT
>aminoacid_sequence_2
MPDFMPCNVSDRTEEEPMSPARTYDFGHKLCVSDFTPMLKKPERT
```
</p>

### How to estimate the fold potential and/or disorder and aggregation propensities
By default, ORFold only estimates the fold potential of the input sequences. 
The disorder and aggregation propensities can be however calculated as well.
The user can specify which calculation methods are to be launched with 
the **-options** argument. 

Each method used by ORFold is referred by its initial: 
<pre>
   HCA     : H
   IUPred  : I
   TANGO   : T 
</pre>

The user must specify the combination of methods he wants to apply
on the input sequences giving their initials with the **-options** argument without any space.


<div class="admonition note">
    <p class="first admonition-title">
    </p>
    <p class="last">

<table>
 <tr>
    <th><b>Methods</b></th> 
    <th><b>Options</b></th>
 </tr>
 <tr>
     <th> &nbsp;&nbsp;&nbsp; HCA <b>&</b> IUPred <b>&</b> Tango &nbsp;&nbsp;&nbsp; </th>
     <th> &nbsp;&nbsp;&nbsp; HIT or HTI or IHT or ITH or THI or TIH &nbsp;&nbsp;&nbsp;</th>
 </tr>
 <tr>
     <th> &nbsp;&nbsp;&nbsp; HCA <b>&</b> IUPred &nbsp;&nbsp;&nbsp; </th>
     <th> &nbsp;&nbsp;&nbsp; HI or IH &nbsp;&nbsp;&nbsp;</th>
 </tr>
 <tr>
     <th> &nbsp;&nbsp;&nbsp; HCA <b>&</b> Tango &nbsp;&nbsp;&nbsp; </th>
     <th> &nbsp;&nbsp;&nbsp; HT or TH &nbsp;&nbsp;&nbsp;</th>
 </tr>
 <tr>
     <th> &nbsp;&nbsp;&nbsp; IUPred <b>&</b> Tango &nbsp;&nbsp;&nbsp; </th>
     <th> &nbsp;&nbsp;&nbsp; IT or TI &nbsp;&nbsp;&nbsp;</th>
 </tr>
</table>
</p>
</div>

### Basic run
The following instructions estimates the fold potential, and the disorder and aggregation propensities of
all amino acid sequences contained in the input fASTA file:

```{python}
orfold -fna sequences.fasta -options HIT
```

<p></p>

It's important to mention that **IUPred** and **Tango** are proposed as 
complementary to the **HCA** methods and they will slow down dramatically ORFold. 
If you wish to calculate only the foldability with HCA you can simply type:
```{python}
orfold -fna sequences.fasta -options H
```
<br>
## Output:
The basic output of ORFold is a table named the same as the input fasta with the extention ".tab".
In the example above it will be : **sequences.tab**
<br>This table contains the information asked by the user (in the -options label) per sequence in the fasta file. 
This table file is a tab separated file easy to be imported to any analysis package (R,excel etc.)



 <p> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;

 Output file example with -option HIT:

```{}
Seq_ID                  HCA     Disord  Aggreg
aminoacid_sequence_1	1.340	0.000	0.230	
aminoacid_sequence_2	-0.230	0.120	0.012	
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
 Output file example with -option H:
```{}
Seq_ID                  HCA     Disord  Aggreg
aminoacid_sequence_1    1.340   nan     nan
aminoacid_sequence_2    -0.230  nan     nan
```
</p>
