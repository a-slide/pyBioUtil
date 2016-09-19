# -*- coding: utf-8 -*-

"""
@package pyBioUtil
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2016
* <aleg@ebi.ac.uk>
* [Github](https://github.com/a-slide)
"""

# Strandard library imports
from gzip import open as gopen
from glob import glob

# Third party packages


#~~~~~~~ Fasta tools ~~~~~~~#

def gencode_fasta_info (gencode_fasta, gencode_info="gencode_info.tsv"):
    """
    Parse the header lines of a gencode transcript fasta file and generate a simple tabulated file containing
    the information for each transcript
    @param gencode_fasta Input file gziped or not, containing the gencode fasta sequences of transcripts 
    @param gencode_info Output file, simple tsv file with the transcript information [DEFAULT:gencode_info.tsv]
    """
    # init empty handlers
    infile = outfile = None
    try:
        # open the fasta file gziped or not
        infile = gopen(gencode_fasta, "rt") if gencode_fasta[-2:].lower()=="gz" else open(gencode_fasta, "r")
        outfile = open (gencode_info, "wt")
        
        # write header in the info file
        h = "GENCODE_transcript_id\tGENCODE_gene_id\tHAVANA_gene_id\tHAVANA_transcript_id\ttranscript_name\tgene_name\tlength\tRNA_type\n"
        outfile.write (h)
        
        nseq=0
        # iterate over the fasta file line and get only the headers
        for line in infile:
            if line[0] == ">":
                # Decompose the line and recompose it in tabulated format
                outfile.write("\t".join([i.strip() for i in line[1:].split("|")])+"\n")
                nseq+=1
                
    # Close the files properly
    finally:
        if infile:
            infile.close()
        if outfile:
            outfile.close()
    
    print ("Found {} Sequences".format(nseq))


def gencode_fasta_clean (gencode_fasta, clean_fasta="gencode_clean.fa.gz"):
    """
    Parse the header lines of a gencode transcript fasta file and rewrite a file with clean header fields containing only the sequence ID
    If you want to keep a track of the additional information you can use gencode_fasta_info beforehand
    @param gencode_fasta Input file gziped or not, containing the gencode fasta sequences of transcripts 
    @param gencode_info Output file, fasta file whith clean sequence header field [DEFAULT:gencode_clean.fa.gz]
    """
    # init empty handlers
    infile = outfile = None
    try:
        # open the fasta file gziped or not
        infile = gopen(gencode_fasta, "rt") if gencode_fasta[-2:].lower()=="gz" else open(gencode_fasta, "r")
        outfile = gopen(clean_fasta, "wt") if clean_fasta[-2:].lower()=="gz" else open(clean_fasta, "r")
        
        nseq=0
        # iterate over the fasta file line and get only the headers
        for line in infile:
            # Decompose the line get only the sequence ID
            if line[0] == ">":
                ls = [i.strip() for i in line[1:].split("|")]
                outfile.write (">{}\n".format(ls[0]))
                nseq+=1
            
            else:
                outfile.write(line)
                
    # Close the files properly
    finally:
        if infile:
            infile.close()
        if outfile:
            outfile.close()
    
    print ("Found {} Sequences".format(nseq))


#~~~~~~~ Fastq tools ~~~~~~~#

def fastqc_summary(fastqc_res_dir, table_if=["pass", "warn", "fail"] , plot_if=["pass", "warn", "fail"], max_table_row=10):
    """
    Summarize and display fastqc results in IPYTHON/JUPYTER. Don't try to use in another interface or in terminal
    Requires PANDAS third party python package to work. Works directly from the zipped results obtained by FastQC v0.11.5+
    @param fastqc_res_dir Directory containing the zipped folders obtained with fastQC
    @param table_if Output the tables for each individual section if in the following list ["pass", "warn", "fail"] [DEFAULT ALL]
    @param plot_if Output the plot for each individual section if in the following list ["pass", "warn", "fail"] [DEFAULT ALL]
    @param max_table_row  Maximal number of row for the tables [DEFAULT 10]
    """
    from IPython.core.display import Image, display, Markdown, HTML
    import pandas as pd
    import zipfile
    import tempfile
    
    # List modules used by FastQC:
    table_modules = ["Basic Statistics", "Overrepresented sequences", "Kmer Content"]
    
    plot_modules = {
        'Per base sequence quality': 'per_base_quality.png',
        'Per tile sequence quality': 'per_tile_quality.png',
        'Per sequence quality scores': 'per_sequence_quality.png',
        'Per base sequence content':'per_base_sequence_content.png',
        'Per sequence GC content':'per_sequence_gc_content.png',
        'Per base N content':'per_base_n_content.png',
        'Sequence Length Distribution':'sequence_length_distribution.png',
        'Sequence Duplication Levels':'duplication_levels.png',
        'Kmer Content':"kmer_profiles.png",
        'Adapter Content':'adapter_content.png'}

    for fp in glob(fastqc_res_dir+"/*_fastqc.zip"):
        display(Markdown("---\n## {}".format(fp.rpartition("/")[2])))
        
        # Extract the zip file in a temporary folder 
        with tempfile.TemporaryDirectory() as tmpdirname:
            zfile = zipfile.ZipFile(fp)
            zfile.extractall(tmpdirname)    
            base_dir = glob(tmpdirname+"/*/")[0]
            
            # Iterate over the fastqdata file
            with open(base_dir+"fastqc_data.txt", "r") as fastqc_data:
                # Set bool flag to false
                found_module = found_header = False
                
                # Parse the file
                for line in fastqc_data:
                    line = line.strip()

                    # If End of section a table section can be rendered
                    if line.startswith(">>END_MODULE"):
                        if found_module and found_header:
                            df_html = df.head(max_table_row).to_html()
                            display(HTML('<font size=2>'+df_html+'</font>'))
                            found_module = found_header = False             
                    
                    # Detect the start of a mmodule section
                    elif line.startswith(">>"):
                        ls = line[2:].split("\t")
                        module_name = ls[0]
                        module_status = ls[1]
                        display(Markdown("**{} : {}**".format(ls[0], ls[1])))
                        
                        # Detect the start of a table section
                        if module_status in table_if and module_name in table_modules:
                            found_module=True
                        
                        # Detect the header of a plot section
                        if module_status in plot_if and module_name in plot_modules:                
                            image = "{}Images/{}".format(base_dir, plot_modules[module_name])
                            display(Image(image, width=500))
                    
                    # Fill a dataframe if a table section was found
                    elif found_module:
                        if not found_header:
                            df = pd.DataFrame(columns=line.split("\t")[1:])
                            found_header = True
                        else:
                            ls = line.split("\t")
                            df.loc[ls[0]] = ls[1:]

############################################################################################################################################
