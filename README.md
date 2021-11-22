# Anova

Implementation of the algorithm for inferring gene regulatory networks by ANOVA [1] and its variation applying
a non-parametric Friedman Test [2] instead of a Two-way ANOVA. 

Tested on MATLAB R2016a with the Parallel Computing Toolbox. 

Required Input 
--------------

* genes = cell array (# of genes,1) of genes names.
* transcription factors = cell array (# of TF, 1)  of TF names.
* expressiondata = numeric matrix (# of genes, # of conditions) of the gene expression data.

Output
------

* Net_a = Inferred network write as 3-column cell array 
	* Column 1 = Regulator 
	* Column 2 = Target Gene 
	* Column 3 = non-linear correlation coefficient derived from an analysis of variance (Two-way ANOVA)

* Net_f = Inferred network write as 3-column cell array 
	* Column 1 = Regulator 
	* Column 2 = Target Gene 
	* Column 3 = non-parametric, non-linear correlation coefficient derived from a Friedman Test 
		 
Bibliography		 
------------

[1] Küffner R, Petri T, Tavakkolkhah P, Windhager L, Zimmer R. Inferring gene regulatory networks by ANOVA. 
    Bioinformatics. 2012;28:1376–82.

[2] Hoffman JIE. Chapter 26 - Analysis of Variance II. More Complex Forms. In: Hoffman JIE, editor. 
    Biostatistics for Medical and Biomedical Practitioners. Academic Press; 2015. p. 421–47. doi:10.1016/B978-0-12-802387-7.00026-3.


License
-------

This project is licensed under the GNU General Public License. For the exact terms please see the [LICENSE file]


