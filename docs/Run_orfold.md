# Running ORFold:

## Basic run:

For its most basic run, ORFold requires only the amino acid sequence(s) in format 
FASTA (given with the **-fna** label). ORFold can also take as input more than one FASTA 
files and will treat them the one after the other independently.
 <p> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
 FASTA file example:
```{}
>aminoacid_sequence_1
AGNVCFGGRTYMPDFDGMSCVNWQERT
>aminoacid_sequence_2
MPDFMPCNVSDRTEEEPMSPARTYDFGHKLCVSDFTPMLKKPERT
```
</p>
ORFold requires also the calculation method(s) to perform (given with the **-options** label).
There are three methods implemented in ORFold: 
<pre>
   HCA     : H
   IUPred  : I
   TANGO   : T 
</pre>

Use the initial of each method (as given above) without any space (like a single word).     


<div class="admonition note">
    <p class="first admonition-title">
        Note
    </p>
    <p class="last">
Based on the methods you want to use, the <b>-options</b> can be:
<br><br>
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
For calculating the foldability potentail of one given fasta file with all the three methods, type:

```{python}
orfold -fna sequences.fasta -options HIT
```

<p></p>

It's important to mention that **IUPred** and **Tango** are proposed as complementary to the **HCA** methods and they will slow down dramatically ORFold. 
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
