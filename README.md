# Comparing_Two_Genotype_Files_(Simplified_Version)
Code Written by M. Joseph Tomlinson IV

Compares to very simplified genotype files for concordance between samples and SNPs.

Code was Developed in the Abasht Laboratory at the University of Delaware under
the supervision of Dr. Behnam Abasht
website: http://canr.udel.edu/faculty/behnam-abasht/

Program Files Provided to Run Code:
comparing_files_code.py
TEST_File_1.txt
TEST_File_2.txt
Matching_Key.txt

Program reads in two genotype files that consist of samples and SNP data. 
Files:
TEST_FILE_1.txt
TEST_FILE_2.txt

The sample files have slightly different names, so the program converts one file's sample names to the correct nomenclature
using a "Matching_Key" file (TEST_FILE_1.txt sample name header is converted to the correct names using Matching_Key.txt)
Example: "TestSample_1" converted to "Sample_1"

The files are then compared for overlapping SNPs, re-ordered and common SNPs between the files printed out in new genotype files. 

Output Files:
Commons_SNPs_File_1.txt (Test file 1's data)
Commons_SNPs_File_2.txt (Test file 2's data)

The sample and SNP genotypes of the common files are compared for concordance, with final results of the
concordance printed out in separate Sample and SNP files for which Samples and SNPs matched and did not match.

Output Files:
Comparing_Samples.txt
Comparing_SNPs.txt

Note: This program was initially developed to test basic functions which were later utilized in a much 
larger 600K genotype dataset that is larger and more complex to analyze. I will admit this overall setup
of "for" loops was later found out to be extremely time consuming when anaylzing such large datasets. However, the code
was utilized in a HPC enviroment, so running for a day was not an issue.    
