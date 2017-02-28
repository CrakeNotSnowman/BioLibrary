My cheap personal library for research
======================================

Author(s): Keith Murray

Contact: kmurrayis@gmail.com

Certain Functions/snipits from:
* Ufuk Nalbantoglu
* Sam Way


~~ ~~ ~~ ~~ ~~

*"Hofstadter's Law: It always takes longer than you expect,*
*even when you take into account Hofstadter's Law."*

*--Douglas Hofstadter*

~~ ~~ ~~ ~~ ~~


## Personal Notes
(http://python-packaging.readthedocs.org/en/latest/minimal.html
for the legwork of the install) 




## Introduction
This is sort of my main library for all of research. 

I've called this library kmbio because it's my personal library from bioinformatic research,
prefix km refers to my name, and bio is the subject the library is focused on. This lets me 
install most any library I'm interested in from the web without major worry of conflicting 
namespaces.




## Requirements
#### Fasta and Multifasta requirements:
* Only use this library for nucleic acid sequences, it does not support protein sequences.
* Only "A", "T", "G", "C" are supported, and all other characters are treated as "N".

#### Python Version 
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

I'm going to skip the code here, since it's mostly the same as above, and lets get all of 
the other plots shown here as well to see how they look

![all_plots](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/Full_Scatter_Plots.png)

### What are these methods?

[Kmer](https://en.wikipedia.org/wiki/K-mer) 
is treated as a 4^k length vector filled with the counts of each specific permutation.

[AMI](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-48) <cite>[2]</cite>
 short for Average Mutual Information (a profile commonly used in signal analysis). An AMI profile
acts as a signature for genomic sequences loosely based on a normalized joint probability between
nucleotides X and Y which are k bases apart. 

DBP: Dictionary Based Profiles. My own method, roughly similar to a kmer vector in that it is a fixed
vector populated with word counts. The words are not all k nucleotides in length, but are much more
flexible. It maintains a linear time complexity by taking advantage of LZW [3] like codebook
construction. 

Dict: This used the same modification to the LZW codebook construction. The construction of the 
codebook follows the rules of LZW: if the word currently in the window is in the codebook, the window
expands by one character. If it is not, it gets added to the dictionary. Unlike LZW, once the word
is added to the dictionary, the head of the window moves one character forward, not to the tail of the
old window. Each sequence builds a codebook in this manner. Once all sequences have a codebook, the 
distance is calculated following a basic grammar distance: One subtracted from the intersection of 
the two dictionaries divided by their union. There are more complex distance measures outlined in 
[4] which may be used in the future. This is a more memory heavy method. 

### One more way to view things

Looks pretty! But it's only 2D! We can use nSpect to kick it up a dimension, and instead of using 
PCA by way of SVD, we can use an iterative approach to Multidimensional Scaling!

```python
# Assuming the code previously shown as been run

# Get the displays formatted for nSpect
nSpectDisplay = [[x, 0, 4] for x in display] 

# Run nSpect
kmbio.nSpect(seqSet.dm, short_labels, nSpectDisplay)
```

This should pop out a window like 

![nSpect_gif](https://github.com/CrakeNotSnowman/BioLibrary/blob/master/examples/nSpect.gif)

and when you're in the window, you can click on the dots to remove them from the visualization,
this will also tell you what the dot is.

#### Some specifics of nSpect
This portion of the kmbio library rarely works. nSpect is a C++ program that runs from the command 
line. This implementation is a very hackish way to bring it into my python library for my convenience.

I can't seem to get this to work in a Jupyter notebook, but my IP[y]thon Qt Console lets it run.
I have not dug in to figure out why. 

nSpect source code and usage instructions can be found [here](http://bioinfo.unl.edu/nspect.php)
At this time it only runs on Linux and Mac OS' 

The compiled version included in this library is for a 32 bit Linux machine.




## You made it this far. Congrats!
~~ ~~ ~~ ~~ ~~

*Hello, Grandmaster Crake. Enter the passnumber now.*
Crake did so. A new sentence popped up: *Adam named the animals. MaddAddam customizes them.*

--From "Oryx and Crake"
by Margaret Atwood

~~ ~~ ~~ ~~ ~~

I recommend Oryx and Crake to anyone who toys with bioinformatics. 


# References 

[1]: Jeanniard, A., Dunigan, D. D., Gurnon, J. R., Agarkova, I. V., Kang, M., Vitek, J., ... & Van Etten, J. L. (2013). Towards defining the chloroviruses: a genomic journey through a genus of large DNA viruses. BMC genomics, 14(1), 158.

[2]: Bauer, M., Schuster, S. M., & Sayood, K. (2008). The average mutual information profile as a genomic signature. BMC bioinformatics, 9(1), 48.

[3]: Welch, T. A. (1984). A technique for high-performance data compression. Computer, 6(17), 8-19.

[4]: Otu, H. H., & Sayood, K. (2003). A new sequence distance measure for phylogenetic tree construction. Bioinformatics, 19(16), 2122-2130.

