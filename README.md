# read-pipeline
Short read pipeline of CIWARS for taxonomy classification and ARG normalization based on rpoB marker genes
# Requirements
<ol>
  <li>Linux operating system</li>
  <li>conda</li>
</ol>

# Installation
<pre>
git clone https://github.com/muhit-emon/read-pipeline.git
cd read-pipeline
bash install.sh
conda env create -f environment.yml
</pre>
# conda environment activation
After installation, a conda environment named <b>read_pipeline</b> will be created.<br>
To activate the environment, run the following command <br>
<pre>
conda activate read_pipeline
</pre>

# Download the compressed Blast Database file from Zenodo (25 GB) to run MetaCompare and uncompress it
Go inside <b>Assembly_Pipeline</b> directory
<pre>
wget https://zenodo.org/records/10471551/files/BlastDB.tar.gz
tar -zxvf BlastDB.tar.gz
</pre>

# Download the compressed DeepARG-DB and mobileOG database (DB.tar.gz) from one drive, put it inside the "Assembly_Pipeline" directory and uncompress it
Go to <a href="https://virginiatech-my.sharepoint.com/:u:/g/personal/muhitemon_vt_edu/EQjIpLhOmMVFotPRDiK5Id0BGQFVSVVbyEWOmRWiz-rYUA">One Drive</a> and download <b>DB.tar.gz</b>. Put it inside the <b>Assembly_Pipeline</b> directory and uncompress it.
<pre>
tar -zxvf DB.tar.gz
</pre>

# Usage on metagenomic paired-end short read data
Go inside <b>Assembly_Pipeline</b> directory. <br> <br>
<b>To run the assembly pipeline on metagenomic paired-end short read data (<span> &#42; </span>.fastq/<span> &#42; </span>.fq/<span> &#42; </span>.fastq.gz/<span> &#42; </span>.fq.gz), use the following command</b> <br>
<pre>
nextflow run assembly_pipeline.nf --R1 &ltabsolute/path/to/forward/read/file&gt --R2 &ltabsolute/path/to/reverse/read/file&gt --out_fname &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>assembly_pipeline.nf</b>) are: <br><br>
<b>--R1</b>: The absolute path of the fastq file containing forward read sequences <br>
<b>--R2</b>: The absolute path of the fastq file containing reverse read sequences <br>
<b>--out_fname</b>: The prefix of the output file name <br><br>

With <b>--out_fname S1</b>, three output files named <b>S1_resistome_risk.txt</b>, <b>S1_ARGs.faa</b>, and <b>S1_ARGs_and_mobility.tsv</b> will be generated inside <b>Assembly_Pipeline</b> directory. <br><br>
