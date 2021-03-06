


# RDDSVM

RDDSVM (RNA-DNA differences with support vector machines) is an open source R program for the prediction of candidate A-to-I RNA editing sites using machine learning (support vector machines) utilizing the sequence information flanking the editing sites.

## INTRODUCTION

Adenosine to inosine (A-to-I) editing in RNA is not only involved in various biological processes like gene expression, alternative splicing and mRNA degradation but also associated with carcinogenesis and various human diseases. Therefore, accurate identification of RNA editing sites in transcriptome is valuable for research and medicine. RNA-seq, a high-throughput sequencing of the transcriptome is very useful for the detection of RNA editing events in condition specific cells. However, computational analysis methods of RNA-seq data have considerable false-positive risks due to mapping errors. 

We developed a simple machine learning (SVM) method using support vector machines to train sequence information derived from experimentally verified A-to-I editing sites in order to predict new A-to-I editing sites in RNA. Rediportal database of A-to-I RNA editing events (Picardi, D’Erchia, Giudice, & Pesole, 2017) was used for the generation of the positive dataset. This database has more than 4.5 million A-to-I events detected in 55 body sites of humans from thousands of RNAseq experiments. Reference genomic sequences (hg19) for the coordinate of 5000 randomly selected A-to-I RNA editing events from the Rediportal database with the flanking sequences of 50nt (in both 3’ and 5’ direction) were obtained for the positive dataset. As for the negative set, we randomly selected 5000 sequences of the same length from the reference genome which contained adenine in the center and not present in the Rediportal database. 
In our model, we evaluated composition of the triplet sequence elements (AUU,AAU,... etc) in each nucleotide sequence. For a sequence, the appearance of 64 different kind of triplets were counted while sliding over the sequence by one nucleotide. Total count of triplets was normalized and provided to SVM classifier as a 64-dimensional vector. For the optimization of the model and the SVM, please refer to the original article. 

SVM classifier showed high performance on experimentally verified data providing sensitivity of 92.8%, specificity of 77.1% and accuracy of 90.2%. The methodology is very easy to interpret and computationally low demanding making it a convenient and valuable choice for facilities with low sources.

## PREREQUISITES

RDDSVM was implemented using R version 3.6.3 and works on both Windows and Unix. The performance of the program was not tested with previous versions. Additionally RDDSVM requires R libraries BSgenome.Hsapiens.UCSC.hg19, e1071, optparse. The code will try to download those packages if they are not already installed.

## DOWNLOAD

You can easily install RDDSVM by the following commands:

    git clone https://github.com/huseyintac/RDDSVM.git

RDDSVM is a stand-alone script and can be simply launched by

    cd RDDSVM
    Rscript RDDSVM.R -f sample_input.csv -o output.txt 

## INPUT FILES
RDDSVM requires  input files as simple textual files with three tabulated columns including:

-   chromosome number or name (1,2,3,....,X, Y)
-   coordinate of the position of A-to-I editing (1-based)
-   strand information (+ or -)

A sample input file is provided in the package "sample_input.csv"

Sample input file:

14	33084977	+

15	31229832	-

9	112781854	+

2	65235663	-

2	198168947	-

## OUTPUT FILES
The output file is going to be similar to the input file with the addition of 2 extra columns for the sequence used for the predictions, and the prediction result of the SVM. The sequence contains the editing site and the flanking sequences of 50nt (in both 3’ and 5’ direction). The result 1 and -1 indicates a prediction of a true editing site and false editing site respectively. 

### Sample output file

14	33084977	+	AAA....CUC	1

15	31229832	-	GUU....UAC	1

9	112781854	+	UAG....AUA	-1

2	65235663	-	AGC....UUA	1

2	198168947	-	GGU....GCAC	1


RDDSVM will also print out SVM training results and query predictions results as a whole. 


### SVM training results:

Error estimation of ‘svm’ using 10-fold cross validation: 0.1364

### Query results:

Number of sites predicted as true editing sites -> 89

Number of sites predicted as false editing sites -> 11

     
## USAGE INFORMATION

Options:

-f CHARACTER, --file=CHARACTER
	input file name

-o CHARACTER, --out=CHARACTER
	output file name [default= out.txt]

-h, --help
	Show this help message and exit
	

An input file should be provided with -f option. By default an output file named "out.txt" will be created to the current folder. The name and location of the output file can be provided with -o option. 

### Example:
    Rscript RDDSVM.R -f sample_input.csv -o output.txt 
