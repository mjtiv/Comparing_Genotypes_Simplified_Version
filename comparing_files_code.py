#!/usr/bin/env python3.6

#Note: Code was Developed in the Abasht Laboratory at the University of Delaware under
# the supervision of Dr. Behnam Abasht
# website: http://canr.udel.edu/faculty/behnam-abasht/

#################Comparing Two Simplified Genotype Files for Concordance###############
######################Code Written by M. Joseph Tomlinson IV################

#Program reads in two genotype files with samples and SNP data. The samples have slightly
#different names, so the program converts one file's sample names to the correct nomenclature
#using a "Matching_Key" file. The files are then compared for overlapping SNPs, re-ordered and
#common SNPs between the files printed out in new genotype files. The sample and SNP genotypes
#of the common files are compared for concordance, with final results of the
#concordance printed out in separate Sample and SNP files.

#Basic Packages Required for the Software
import collections

####FUNCTIONS#############################################
#1 id_converter(header_info, ID_dictionary_file)    

#ID CONVERTER FUNCTION,
#Function accepts: the "header" of a file alreadly parsed and "key txt file"
#and converts the header IDs to the values from the key file
#allowing for sample IDs that are slightly different to be matched up
#with universal sample IDs for analysis. The function then ultimately
#returns a new header that can be later utilized in analysis
def id_converter(header_info, ID_dictionary_file):
    file_key = open(ID_dictionary_file, 'r') #opens dictionary file
    ID_dict={}
    for line in file_key: #parses through dictionary file
        if line.startswith("Values"):
            pass
        else:
            data=split_line(line)
            key=(data[1]).strip('\n') #strips new lines
            values=data[0]
            ID_dict.update({key:values})
    file_key.close()
    header_names=split_line(header_info.strip('\n')) #need to stip new line from original header
    for x in range(len(header_names)): #starts parsing through changing IDs
        if (header_names[x]) in ID_dict:
            value=ID_dict[header_names[x]]
            header_names[x]=value
        else:
            pass
    new_header=('\t'.join(header_names)+"\n") #converts header to a string
    return(new_header)

#SPLIT LINE FUNCTION
#Function just parses apart a line into a list
def split_line(line):
    values = line.split("\t")
    return (values)


#Parsing through file for header and SNP data
#Function accepts a file and starts to parse through the data
#returning the results as a retrievable dictionary (header and dictionary of SNPs)
#Dictionary of SNPs can be used to compare common SNPs between files
#Header is original header of data
#Key=SNP Values=Genotypes of Samples
def retrieve_data (file_name):
    file= open(file_name, 'r')
    file_dic={}
    for line in file:
        if line.startswith("IDS"):
            header_info=(line)
        else:
            values = line.split("\t")
            key = (values[0])
            values=line
            file_dic.update({key:values})
    return{'header_info':header_info, 'file_dic':file_dic}


#Definition loads in the header and dictionary files from the parsed datasets and begins ordering
#the SNPs and probes using nested foreloops, extremely slow and in-efficient and takes time to run
#if running on a large dataset
def create_common_SNPs_files(header_info_file_1, header_info_file_2, file_1_dic, file_2_dic):
    sorted_file_2_dic=collections.OrderedDict([(k,file_2_dic[k]) for k in sorted(file_2_dic.keys())])
    
    #Opening New Files To Create Matching SNPs Files
    Common_SNPs_File_1=open("Commons_SNPs_File_1.txt","w")
    Common_SNPs_File_2=open("Commons_SNPs_File_2.txt","w")

    Common_SNPs_File_1.write((header_info_file_1))
    Common_SNPs_File_2.write((header_info_file_2))
    #Nested for loop begins sorting and matching SNPs
    for x in (sorted_file_2_dic):
        for y in (file_1_dic):
            if x == y :
                Common_SNPs_File_1.write(file_1_dic[x])
                Common_SNPs_File_2.write(file_2_dic[x])
            else:
                pass
    Common_SNPs_File_2.close()
    Common_SNPs_File_1.close()

##################Setup another function to sort individuals....###########################


#definition returns a dictionary of samples with corresponding genotype data
#allowing for the samples genotypes to be compared and tabulated
#for what genotypes match and which ones are different
def samples_data (file_name):
    file = open(file_name, 'r')
    file_samples_dic={}
    for line in file:
        if line.startswith("IDS"):
            samples_IDs=line.split("\t")
            samples_IDs_len=len(samples_IDs)
            for i in range(1,samples_IDs_len,1):
                key=(samples_IDs[i])
                sample_values=[]
                file.seek(0)  ###Resets code to begining of file
                for line in file:
                    if line.startswith("IDS"):
                        pass
                    else:
                        values = line.split("\t")
                        #print (values)
                        sample_values.append(values[i])
                file_samples_dic[key]=sample_values
            pass
        return(file_samples_dic)


#Function retrieves the list of SNPs from a SNP list file
#Returns the list with all the SNPs
def SNPs_List(file_name):
    file = open(file_name, 'r')
    SNPs_List=[]
    for line in file:
        if line.startswith("IDS"):
            pass
        else:
            line_values = line.split("\t")
            SNPs_List.append(line_values[0])
            pass
    return(SNPs_List)

#Function accepts sample file dictionaries and tallies the number of matches and non-matches
#Between samples. The final results are outputed into a Comparing Samples file
def matching_Samples(file_1_samples_dic, file_2_samples_dic):
    Matching_Values_File=open("Comparing_Samples.txt","w")
    Matching_Values_File.write("Sample\t Match_Count\t Non_Match_Count \n")
    total_match_counter=0
    total_no_match_counter=0
    for x in (file_1_samples_dic):
        #print (x)
        
        for y in (file_2_samples_dic):
            if x == y :
                list_1=(file_1_samples_dic[x]) 
                list_2=(file_2_samples_dic[y])
                sample_match_counter=0
                sample_non_match_counter=0
                for k in range(len(list_1)):
                    if list_1[k]==list_2[k]:
                        total_match_counter +=1
                        sample_match_counter +=1
                    else:
                        total_no_match_counter +=1
                        sample_non_match_counter+=1
                printing_values=(str(x.strip('\n')),str(sample_match_counter),'\t',str(sample_non_match_counter),'\n')
                Matching_Values_File.write("\t".join(printing_values))
            else:
                pass
    print ("The SAMPLE match counter is:", total_match_counter)
    print ("The SAMPLE non-match counter is:", total_no_match_counter) 
    Matching_Values_File.close()



#Counter counts the number of SNPs that match and do not match between the
#genotyping panels. All results are outputted into a Comparing_SNPs.txt file
#Inputs is a Common_SNPs_file (either 1 or 2 will do, because they are the same)
#The 2nd and 3rd inputs are the file sample dictionaries
def matching_SNP_Counter(Common_SNPs_file, file_1_samples_dic, file_2_samples_dic):
    
    file_1_SNPs=SNPs_List(Common_SNPs_file)
    Matching_SNPs_File=open("Comparing_SNPs.txt","w")
    Matching_SNPs_File.write("SNPs\t Match_Count\t Non_Match_Count\n")
    #Setting up a Dictionary of SNPs with values of zero
    match_SNP_dictionary={key:0 for key in file_1_SNPs}
    #print (match_SNP_dictionary)
    non_match_SNP_dictionary={key:0 for key in file_1_SNPs}
    #Setting up Counters for SNPs
    total_match_counter=0
    total_no_match_counter=0
    for x in (file_1_samples_dic):
        for y in (file_2_samples_dic):
            if x == y :
                list_1=(file_1_samples_dic[x])
                list_2=(file_2_samples_dic[y])
                sample_counter=0
                for k in range(len(list_1)):
                    if list_1[k]==list_2[k]:
                        key_value=file_1_SNPs[k]
                        new_key_value=(match_SNP_dictionary[key_value]+1)
                        match_SNP_dictionary[key_value]=new_key_value
                        total_match_counter +=1
                    else:
                        key_value=file_1_SNPs[k]
                        new_key_value=(non_match_SNP_dictionary[key_value]+1)
                        non_match_SNP_dictionary[key_value]=new_key_value
                        total_no_match_counter +=1
            else:
                pass
    #Printing the Final SNP results that were collected in dictionaries        
    for key,value in match_SNP_dictionary.items():
        printing_line=(str(key),str(value),str(non_match_SNP_dictionary[key]),"\n")
        Matching_SNPs_File.write("\t".join(printing_line))

    Matching_SNPs_File.close()

    print ("The SNP match counter is:", total_match_counter)
    print ("The SNP non-match counter is:", total_no_match_counter)

def main():
    #Opening the files and Retrieving the data from file 1
    file_1_data=retrieve_data('TEST_File_1.txt')
    header_info_file_1=file_1_data['header_info']
    header_info_file_1=id_converter(header_info_file_1, 'Matching_Key.txt')
    
    file_1_dic=file_1_data['file_dic']
    
    file_2_data=retrieve_data('TEST_File_2.txt')
    header_info_file_2=file_2_data['header_info']
    file_2_dic=file_2_data['file_dic']
    
    #This function takes the above parsed data and compares the probes to create a common probes files (probes in different orders)
    create_common_SNPs_files(header_info_file_1, header_info_file_2, file_1_dic, file_2_dic)

    file_1_samples_dic=samples_data('Commons_SNPs_File_1.txt')
    
    file_2_samples_dic=samples_data('Commons_SNPs_File_2.txt')

    #Count How Many of the Samples Match
    matching_Samples(file_1_samples_dic, file_2_samples_dic)

    #Count How Many of the SNPs Match
    matching_SNP_Counter('Commons_SNPs_File_1.txt', file_1_samples_dic, file_2_samples_dic)

main()










