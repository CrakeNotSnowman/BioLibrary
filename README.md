My cheap personal library for research
======================================

Author(s): Keith Murray

Contact: kmurrayis@gmail.com

Certain Functions/snipits from:
* Ufuk Nalbantoglu
* Sam Way




## Personal Notes
(http://python-packaging.readthedocs.org/en/latest/minimal.html
for the legwork of the install) 




## Introduction
This is sort of my main library for all of research. 
Every evaluation I do uses functions wrapped up in here.

I've called this library kmbio because it's my personal library from bioinformatic research,
prefix km refers to my name, and bio is the subject the library is focused on. This lets me 
install any library I'm interested in from the web without major worry of conflicting namespaces.




## Requirements
#### Fasta and Multifasta requriements:
* Only use this library for nucleic acid sequences, it does not support protein sequences.
* Only "A", "T", "G", "C" are supported, and all other characters are treated as "N"

#### Language
The Code is only python 2.7 supported, not other versions have been tested, and it is 
not compatible with with python 3.x

*Eventually it will be ported to 3, but it will only support 3.6 and up
    Please note that eventually is likely to be almost never unless there's a strong force
    to migrate this hobby project. Sorry 2.x bashers, I don't care, and enough of this code
    is dictionary based that I want to optimize the migration, not just do a simple migrate.

#### Library Dependencies
* biopython
* ete3
* numpy
* scipy
* scikit-learn

## Basic Usage
To explore the usage of the library, we'll focus on sequences outlined in 
*Towards defining the chloroviruses: a genomic journey through a genus of large DNA viruses* <cite>[1]</cite>
with the exception of the Ostreococcus virus set. 

Here is ML tree "based on a concatenated alignment of 32 core protein families" for reference:
![bio_trees](https://static-content.springer.com/image/art%3A10.1186%2F1471-2164-14-158/MediaObjects/12864_2012_Article_4824_Fig1_HTML.jpg)

```python
# Import Libraries
import kmbio
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt

# Sample Set Prep (Specific to this set)
# Get labels
fl = open("tests/nSpectLabel.txt")
short_labels = [line.strip() for line in fl]
fl.close()
# Get Display values
fl = open("tests/nSpectDisplay.txt")
display = [float(line[0]) for line in fl]
fl.close()

# Grab Sample Set
fl = "tests/Virus_test_set.fasta"

# Assign header
a, b = kmbio.multifna_read(fl)

# Group into setSequences data structure
seqSet = kmbio.setSequences(short_labels, b)

# Generate profiles
seqSet.getKmer(7)
seqSet.getAMI(64)
seqSet.getDBP1()
seqSet.getDict()


# Get Representations and Visualizations
seqSet.distMat("kmer", "c")

#   Build a tree
kmbio.getTreeUPGMASet(seqSet)
```
<object data="https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/kmer7_Corr_UPGMA_Tree.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/kmer7_Corr_UPGMA_Tree.pdf">
        Rendering PDFs is difficult. Please download the PDF to view sample Tree: <a href="https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/kmer7_Corr_UPGMA_Tree.pdf">Download PDF</a>.</p>
    </embed>
</object>
Here is the pdf rendered as a .png:
![kmer7_Tree](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/kmer7_Corr_UPGMA_Tree.png)



For the plots, we will adhere to the following representations, 

* Green:  SAG
* Blue:   Pbi
* Yellow: NC64A

```python
#   Plot Prep
display = np.array(display)
#   For plot color
blue = display == 1
yellow = display == 3
green = display == 2

# Get Distance Matrix
seqSet.distMat("kmer", "c")
pca = PCA()
kmer_pca = pca.fit_transform(seqSet.dm)

# Lets get plotting! 
plt.plot(kmer_pca[blue, 0],kmer_pca[blue, 1], 'bo')
plt.plot(kmer_pca[yellow, 0],kmer_pca[yellow, 1], 'yo')
plt.plot(kmer_pca[green, 0],kmer_pca[green, 1], 'go')
plt.title("PCA on Distance Matrix of 7mer Profiles")
plt.show()
```
![kmer7_plot](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/kmer7_Corr_Scatter_Plot.png)

Let's add labels to double check accuracy
```python
plt.plot(kmer_pca[blue, 0],kmer_pca[blue, 1], 'bo')
plt.plot(kmer_pca[yellow, 0],kmer_pca[yellow, 1], 'yo')
plt.plot(kmer_pca[green, 0],kmer_pca[green, 1], 'go')
plt.title("PCA on Distance Matrix of 7mer Profiles")
for i in range(len(short_labels)):
    plt.annotate(short_labels[i], (P[i,0], P[i,1]))
plt.show()
```
![kmer7Label_plot](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/labeled_kmer7_Corr_Scatter_Plot.png)

and lets get all of the other plots shown here as well

![all_plots](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/Full_Scatter_Plots.png)




[1]: Jeanniard, A., Dunigan, D. D., Gurnon, J. R., Agarkova, I. V., Kang, M., Vitek, J., ... & Van Etten, J. L. (2013). Towards defining the chloroviruses: a genomic journey through a genus of large DNA viruses. BMC genomics, 14(1), 158.


