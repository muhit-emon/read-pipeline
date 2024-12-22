# read-pipeline
Short read pipeline of CIWARS for taxonomy classification of bacterial pathogens and computation of ARG abundance based on rpoB marker gene normalization
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

# Download the compressed standard 16 GB kraken2 DB and uncompress it
Go inside <b>read-pipeline</b> directory and create a folder named k2_DB
<pre>
mkdir k2_DB
cd k2_DB
</pre>

Now download the kraken2 DB inside this folder and uncompress it
<pre>
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240112.tar.gz
tar -zxvf k2_standard_16gb_20240112.tar.gz
rm k2_standard_16gb_20240112.tar.gz
</pre>

# Usage on metagenomic paired-end short read data
Go inside <b>read-pipeline</b> directory. <br> <br>
<b>To run the short read pipeline on metagenomic paired-end short read data (<span> &#42; </span>.fastq/<span> &#42; </span>.fq/<span> &#42; </span>.fastq.gz/<span> &#42; </span>.fq.gz), use the following command</b> <br>
<pre>
nextflow run short-read-pipeline.nf --R1 &ltabsolute/path/to/forward/read/file&gt --R2 &ltabsolute/path/to/reverse/read/file&gt --out_fname &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>short-read-pipeline.nf</b>) are: <br><br>
<b>--R1</b>: The absolute path of the fastq file containing forward read sequences <br>
<b>--R2</b>: The absolute path of the fastq file containing reverse read sequences <br>
<b>--out_fname</b>: The prefix of the output file name <br><br>

With <b>--out_fname S1</b>, output files named <b>S1.k2report</b>, <b>S1_rpoB_ARG_norm.tsv</b>, and <b>S1_drug_wise_rpoB_norm.tsv</b> will be generated inside <b>read-pipeline</b> directory. <br><br>
