# MohammedA.A_Portfolio
Data Science Portfolio

# Project 1: Average Annual Weather Information for the City of Seattle 

 
1. Analyzed and procured data regarding How many days out of the year it actually rained.
2. Analsis of the average and standard deviaon of precipitaon (in inches) when it does rain.
3. Days of the year it rains per month.
4. Total inches does it rain per month.
5. Minimum and maximum recorded temperature per month.
6. Found the average minimum and maximum recorded temperatures on days that have rained and compared to those that haven’t had rain.




# Project 2: Chicago CTA Ridership Data


 1. Created an outputted file that found the annual ridership for the Loyola stop.
 2. Found the date during the COVID-era when Loyola University closed, by calculating the daily ridership for the Loyola stop from February 1, 2020 through March 30, 2020. 
 3. The overall impact of COVID on CTA ridership. Calculated with the monthly ridership for each of the CTA lines for the last 5 years as comparison, as well as calculated the increase/decrease of ridership for each year relative to the prior year.
 4. Relative to 1/1/2019-6/30/2019, found how much more/less money the CTA made during the first 6 months of 2020.



# Project 3: Visualization and Analysis of Different types Avacado sales in the United States of America during the years of 2015, 2016, & 2017 



1. Visualization of the sales of “conventional” and “organic” avocados for 2015, 2016, and 2017.
2. Visualization of the cities that bought the most avocados for this three year period. 
3. Analysis of the fluctuation of the price of avocados monthly over the three year period and seeing if these trends differ between these three cities and comparing them for any given time within the the time frame.



# Project 4: HCMV Transcriptome Analysis Pipeline
This is a pipeline for analyzing the transcriptome of Human Cytomegalovirus (HCMV) using specific sequencing data. The only neccesary tools to run this is Python.

# Download the Data
The first step is to download the RNA sequencing data for four samples from the Sequence Read Archive (SRA) using wget. We are going to use python scripts
```python
import os #executes system commands
#download SRR5660030, SRR5660033, SRR5660044, SRR5660045 files from the SRA database
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030') #os. converts from python
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')
```

# Preprocess the Data and Gather Test Data

The next step is to preprocess the data using fastq-dump to convert the SRA files to FASTQ format. The neccesary tools to run this are, Python and SRA Toolkit. subsets of the files are used to test, modify with "SRR5660030_1_subset.fastq..." when testing

```python

import os

#split the forward and reverse reads and convert the SRA files into FASTQ format using "-I" and "--split"
os.system('fastq-dump -I --split-files SRR5660030')
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')

# This function to creates small subsets of reads that are used for testing
def create_test_subset(input_file, output_file, num_reads):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for _ in range(num_reads * 4):  # Each read has 4 lines
            line = infile.readline() #files
            outfile.write(line)

# Create small subsets for each sample
samples_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
num_reads = 1000  #the number of reads in the subset, can be changed

for sample in samples_list: # Test Data
    create_test_subset(f'{sample}_1.fastq', f'{sample}_1_subset.fastq', num_reads)
    create_test_subset(f'{sample}_2.fastq', f'{sample}_2_subset.fastq', num_reads)




```

# Build the Reference Genome
Use the HCMV genome as the reference genome for analysis. Download the genome from NCBI and use bowtie2-build to build the index. The neccesary tools to run this are, Python, Biopython, and Bowtie2.
```python
import os  
import Bio
from Bio import Entrez
#fetches HCMV genome sequence from NCBI database
Entrez.email = "mohammed16alsawi@gmail.com" #email for Entrez module
handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype='fasta') #retrieves genome sequence from NCBI 
fasta = handle.read() 
handle.close()
with open('HCMV.fasta', 'w') as f: #save the HCMV genome sequence as a fasta file
    f.write(fasta)
#builds the Bowtie2 index using the HCMV genome sequence
os.system('bowtie2-build HCMV.fasta HCMV')
#run Bowtie2 alignment for all four samples
print('starting bowtie2 .sam') #displays the start of alignment process 'start bowtie2.sam' 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')



```

# Map the Reads to the Reference Genome
Use bowtie2 to map the reads to the reference genome and filter out the unmapped reads. The neccesary tools to run this are, Python and Bowtie2.
```python
import os

# created dictionary to store initial read counts
initial_counts = {}
#creates list of sample names
samples_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

#loop through each sample in the list
for s in samples_list:
    with open(f'{s}_1.fastq') as file1: # opens the file with the sample name and '_1.fastq' ext.
        initial_counts[s] = sum([1 for line in file1]) / 4 #gets the initial read count by counting lines in the file and dividing by 4
    print(f'{s} has {initial_counts[s]:,} read pairs before filtering')

# Running Bowtie2 on samples and keeping only HCMV index mapped reads
for s in samples_list:
    output_sam = f'{s}_mapped.sam'
    os.system(f'bowtie2 --quiet -x HCMV -1 {s}_1.fastq -2 {s}_2.fastq -S {output_sam}')

    # Filtering unmapped reads by reading the SAM file and writing mapped reads to a FASTQ file
    post_filter_counts = 0
    
    #opens the output SAM file and a new FASTQ file for writing
    with open(output_sam) as samfile, open(f'{s}_HCMV.fastq', 'w') as fastqfile: # loops through each line in the SAM file
        for line in samfile:
            if line.startswith('@'): #skips the header lines
                continue

            parts = line.split('\t') #splits lines into columns
            flag_val = int(parts[1]) # gets flag value from the second column

            if flag_val & 4 == 0: # check if the read is mapped ('flag value & 4 is 0')
                post_filter_counts += 1 #runs increments of the count of filtered reads
                qname_val = parts[0] # retrieves the values from the SAM file columns
                seq_val = parts[9]
                qual_val = parts[10]
                fastqfile.write(f"@{qname_val}\n{seq_val}\n+\n{qual_val}\n")  # write the mapped read to the new FASTQ file

    post_filter_counts /= 4 # divides the post-filter count by 4 to get the number of all read pairs
    print(f'{s} has {post_filter_counts:,} read pairs after filtering')

    # appends the read counts to a file
    with open('PipelineProject.log', 'a') as logfile: # appends the read counts to the log file
        logfile.write(f'{s} has {initial_counts[s]:,} read pairs before filtering and {post_filter_counts:,} read pairs after filtering.\n')






```
# SPAdes assembly 
Use SPAdes to perform a genome assembly of four different transcriptomes (SRR5660030, SRR5660033, SRR5660044, SRR5660045) and also append the SPAdes command to the log file. The neccesary tools to run this are, Python and SPAdes.
```python
import os
# define the dictionary to store sample names and their corresponding FASTQ files
samples = {
    'SRR5660030': ('SRR5660030_HCMV.fastq'),
    'SRR5660033': ('SRR5660033_HCMV.fastq'),
    'SRR5660044': ('SRR5660044_HCMV.fastq'),
    'SRR5660045': ('SRR5660045_HCMV.fastq'),
}

#assemble all four transcriptomes together using SPAdes

spades_input = '' #initializes the empty string to store input arguments for the SPAdes command
for index, (sample_name, fq) in enumerate(samples.items(), start=1): # loop through all of the samples by extracting the sample name and FASTQ file
    # appends the current sample's input flag and FASTQ file to the spades_input string
    spades_input += f'--s{index} {fq} '

spades_command = f'spades.py -k 77,99,127 -t 2 --only-assembler {spades_input.strip()} -o HCMV_SRR_assembly'
os.system(spades_command)

# wite the SPAdes command to the log file
with open('PipelineProject.log', 'a') as log_file:
    log_file.write(f'SPAdes command: {spades_command}\n')
    
    
```



# Assembling the Longest Contig and Running BLAST+ Search
Parse the file containing contigs in FASTA format and count the number of contigs longer than 1000 bp and then calculate the total assembly length and write the results to the log file. The neccesary tools to run this are, Python and Biopython. 

 ```python
 
 from Bio import SeqIO

# define the assembly file and makes the input file containing contigs in FASTA format
assembly_file = 'contigs.fasta'

# initialize counters
contigs_gt_1000 = 0 # initialize counter for the number of contigs > 1000 bp
assembly_length = 0 # initialize a variable to store the total assembly length

#Iterate's through contigs in the assembly file
for record in SeqIO.parse(assembly_file, 'fasta'): 
    # calculates the length of the current contig
    contig_length = len(record.seq)
    
    # This checks to see if contig length is greater than 1000 bp
    if contig_length > 1000:
        contigs_gt_1000 += 1  
        assembly_length += contig_length #adds the contig length to the total assembly length

# opens the log file in append mode to write the results
with open("PipelineProject.log", 'a') as log_file:
    log_file.write(f'There are {contigs_gt_1000} contigs > 1000 bp in the assembly.\n') #write's the number of contigs longer than 1000 bp
    log_file.write(f'There are {assembly_length} bp in the assembly.\n') # writes the total assembly length



```
# Analyzing HCMV Transcriptomes Using SPAdes and BLAST

Read the SPAdes assembly output file, find the longest contig, and save to a new file. Then uses the longest contig as a query to perform a BLAST+ search against the nr nucleotide database and saves the results in an XML file. The neccesary tools to run this are, Python, Biopython, and BLAST+.

```python
 
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML #from the BLAST+ ppt. 

#reads the SPAdes assembly output and find the longest contig
assembly_file = "contigs.fasta"  #inputs the assembly file
contigs = SeqIO.parse(assembly_file, "fasta") #parses the input file as FASTA format
longest_contig = max(contigs, key=lambda x: len(x)) #finds the longest contig

#saves the longest contig to a file
with open("longest_contig.fasta", "w") as output_handle: #opens the output file to be appended 
    SeqIO.write(longest_contig, output_handle, "fasta") #writes the longest contig to the output file

#use the longest contig as a BLAST+ input to query the nr nucleotide database
query_sequence = str(longest_contig.seq) # convert the longest contig sequence to a string
blast_program = "blastn" #uses blastn program to use
blast_database = "nr" #searches through nr database

#defines the search parameters
search_parameters = {
    "entrez_query": "txid10357[Organism:exp]",  # this is the betaherpesvirinae taxonomy ID
    "hitlist_size": 10,  # the number of results to return
}

result_handle = NCBIWWW.qblast( #performs BLAST search using the parameters
    blast_program,
    blast_database,
    query_sequence,
    **search_parameters
)

# Save the result to a file
with open("blast_result.xml", "w") as output_file: 
    output_file.write(result_handle.read()) #write the BLAST result to the output file

result_handle.close() #close the result handle 

#Output file in GitHub as "blast_result.xml"
```

# Perform BLAST search on longest contig from SPAdes assembly

Parse the BLAST result XML file using NCBIXML module and extract the top 10 alignments for each. Write into to the log file, including the accession number, percentage identity, alignment length, query start and end, subject start and end, bitscore, e-value, and alignment title. The neccesary tools to run this are, Python, Biopython, and BLAST+.

```python

from Bio.Blast import NCBIXML #biopython library modules

#Opens and parses the BLAST result XML file
with open("blast_result.xml", "r") as result_handle:
    blast_records = list(NCBIXML.parse(result_handle))

#writes the top 10 hits to a log file
with open("PipelineProject.log", "w") as log_file:
    # Write the header row
    log_file.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

    # loop through each records in the BLAST results
    for record in blast_records:
        for alignment in record.alignments[:10]:  # Limit to top 10 hits and loops through the top 10 alignments for each record
            hsp = alignment.hsps[0]  # Only keep the best high scoring segment pair 

            # Extract the information and format the output
            output_line = f"{alignment.accession}\t{hsp.identities * 100 / hsp.align_length:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}\n"

            # Write the output line to the log file
            log_file.write(output_line)


#Output file in GitHub as "PipelineProject.log"
```



# Project 5: Peptide Folding Classifier

We have created here a classifier that will predict if a protein sequence will fold. In this repository we test a variety of non-deep machine learning models and a deep learning model to evaluate the folding propensity of intrinsically disordered proteins (IDP)s. ML Models included in this repository are support vector machines (SVM), decision trees (DT), random forests (RF), k-nearest neighbors (KNN), and a deep learning feedforward neural network. Simulated IDP data was provided by Dr. Peter Kekenes-Huskey for use in training models. Scikit-learn was used to train the ML models on the datasets. 

The data.py script generated training and testing subsets of data for use, which are stored in the 'splitData' folder. The datasets we input into data.py can be found in the 'rawData' folder, but data.py can be run on any other data you wish to use for testing these models.

# Installing Packages

The necessary packages are listed in the 'requirements.txt' file and can be installed together by calling that file with your installer tool. Otherwise each tool can also be installed individually/manually.
1. First, either install locally or create a python virtual environment or conda environment
  - To create a python virtual environemnt in linux/unix (nice way to manage project dependencies versions without mixing with other project dependencies versions)
    - ensure you have python3-venv, and if not install with the following command: **apt-get install python3-venv**
    - create a virtual environment: **python3 -m venv [name of environemnt]**
    - activate the viritual environment: **source venv/bin/activate**
2. Then, install packages with the following command: **pip3 install -r requirements.txt**

# Running the scripts 

Each example below assumes you run from the main project directory; if not please change paths on the commands.

## Prepare Data

To generate splitData from a different dataset, the data.py script can be run to split and properly format datasets for model testing. 

```
nohup python3 scripts/data.py [-h] -i INPUT [-o OUTPUT] [-s SPLITPERCENT] [-t THRESHOLD] [-d [DATA ...]] -p PREDICTINGCOLUMN [-n NEWLABEL] &
```

- options:
  - -h, --help          show this help message and exit
  - -i INPUT, --input INPUT
                        input path to data folder
  - -o OUTPUT, --output OUTPUT
                        output folder path for test/train csvs
  - -s SPLITPERCENT, --splitpercent SPLITPERCENT
                        percentage of data to store as test
  - -t THRESHOLD, --threshold THRESHOLD
                        scoring threshold
  - -d [DATA ...], --data [DATA ...]
                        list of column names of data to input
  - -p PREDICTINGCOLUMN, --predictingColumn PREDICTINGCOLUMN
                        name of column to predict off of
  - -n NEWLABEL, --newLabel NEWLABEL
                        new label for scoring column

Alternatively, the model scripts can be run on the provided sample datasets in 'splitData' by calling the folder as the input. The raw datasets provided can be found in the 'data' folder. These datasets were the input files for the 'data.py' script which generated the training and testing data for model prediction. Replication of the experiment should be done using the files produced, which are the .csv files in the 'splitData' folder.

## Run one model at a time

Example test run of a model script:

```
nohup python3 scripts/svm.py -i splitData/ -j params/svm_params.json -o models/svmTEST.pkl -r results/results.csv -m results/svmTESTCM.png &
```

For more customizable run, please refer to the following parameters:

```
nohup python3 scripts/svm.py [-h] -i INPUT [-c CROSSFOLDS] -j JSON -o OUTPUT -r RESULTS [-n NUMPROCESSORS] -m MATRIX &> logFiles/nohupSVM.out &
```
  - options for non-deep learning scripts (SVM, KNN, DT, RF):
    - -i INPUT, --input INPUT
                          input folder path for data (this is the splitData folder)
    - -c CROSSFOLDS, --crossfolds CROSSFOLDS
                          number of crossfolds (crossfold validation number of random splits for training across entire dataset)
    - -j JSON, --json JSON  
                          json file with parameters for grid search cv (these are in the params folder and depend per model as hyperparameters differ; will test all comibinations)
    - -o OUTPUT, --output OUTPUT
                          output pickle file path (store model object for future use and analysis)
    - -r RESULTS, --results RESULTS
                          path to results csv (this will store f1 score and basic model information)
    - -n NUMPROCESSORS, --numProcessors NUMPROCESSORS
                          number of processers (this is how many processers to use; recommend more because takes long time)
    - -m MATRIX, --matrix MATRIX
                          confusion matrix path (stores confusion matrix for basic evaluation)
```
nohup python3 scripts/feedforward.py [-h] -i INPUT [-o OUTPUT] -r RESULTS [-m MATRIX] [-t TUNINGHYPERPARAMETERS] [-f FOLDS] [-c CURVE] &
```
  - options for deep learning script (feedforward):
    - -h, --help            show this help message and exit
    - -i INPUT, --input INPUT
                        input folder with data (should include scaled/normalized x_train, x_test, y_test, y_train CSVs
    - -o OUTPUT, --output OUTPUT
                        output file for model object; must be pkl
    - -r RESULTS, --results RESULTS
                        path to results.csv (must have Name, Description, Metric, Path)
    - -m MATRIX, --matrix MATRIX
                        path to matrix file; must be png
    - -t TUNINGHYPERPARAMETERS, --tuninghyperparameters TUNINGHYPERPARAMETERS
                        json file with lists of hyperparameters for each option (update with types that can be played
                        around with)
    - -f FOLDS, --folds FOLDS
                        how many folds for crossfold validation selection
    - -c CURVE, --curve CURVE
                        file path to store ROC curve

## Run all default models at once (probably around minimum of 25 processors)

```
./example.bash
```
Ensure this file is executable. On linux, you can make a file executable with the following command:
```
chmod +x path/to/example.bash
```

## Output

Results are stored in the results folder and model object in the models folder. Confusion matrices and ROC curves can be saved wherever user prefers but have dedicated folders in results. The results.csv must already be created with the following headers exactly: **Name,Description,AUC,F1score,Accuracy,crossfoldScore,Path,Importances**

## Parameters

As explained above in the options, each model takes an input of its respective .json paramaters file. Viewing the .json file in the 'params' folder will show the parameters we have set for our data. If you wish to change the parameters, you can do so by editing the model's .json file or by creating an entirely new one and calling that file with the '-j' flag when running a model script.

# Future Work
- Combine ROC curves for multiple models
- Create results.csv without having user first create it
- Add processor adjustment flag to deep learning model
- Allow for even greater model flexibility (i.e. different sized layers) in deep learning
- Recurrent Neural Networks for Timeseries Prediction
- Store weights for deep learning instead of pickle object
