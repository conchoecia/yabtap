"""
Author: Darrin Schultz @conchoecia
File: Transcriptome assembly, annotation, and db creation

Instructions:
  - To run this script and all of its analyses: make sure that python 3,
    snakemake, and biopython are installed on your Unix computer.
  - Execute the following command: `snakemake --cores 45`, replacing `45`
    with the number of threads available on your machine.
"""
import os
import sys
from Bio import SeqIO
import shutil

configfile: "config.yaml"
trimmomatic = "/usr/local/bin/Trimmomatic-0.35"
maxthreads = 90
dockeruser = "dschultz"
curr_dir = os.getcwd()
#print(curr_dir)

config["rna_f"] = {}
config["rna_r"] = {}
config["input_reads"] = []


# Now assert that the sample IDs are max 12 characters
offenders = []
for sample in config["samples"]:
   thisid = config["samples"][sample]["id"]
   if len(thisid) > 12:
      offenders.append("{} ID:{}".format(sample, thisid))
if len(offenders) > 0:
   print( "READ BELOW for instructions on how to fix your config file.")
   print("The following sample IDs are too long. Max 12 chars")
   for element in offenders:
      print("  - {}".format(element))
   print("""\nTip: When shortening the sample IDs to 12 chars max, make sure that the 
         ID doesn't contain a \"_\" char. Also, ensure that the sample ID and parameters ID
         are exactly the same.

         ALLOWED - Sample ID and params ID match:
            POL_Poeo_meseres_H678-98:
              id: "H678-98"
            CTE_Horm_cali_U789089
              id: "U789089"
         NOT ALLOWED:
            POL_Poeo_meseres_H678_98: (Underscore in sampleid)
              id: "H678_98"
            CTE_Horm_cali_U789089: (Mismatch of IDs)
              id: "U7_89089"
         """)

   sys.exit()

# Now assert that the sample IDs don't have any underscores.
#  We care because we use this information to parse sample names later on
#  and delimit using underscores.
offenders = []
for sample in config["samples"]:
   thisid = config["samples"][sample]["id"]
   if "_" in thisid:
      offenders.append("{} ID:{}".format(sample, thisid))
if len(offenders) > 0:
   print("The following sample IDs contain the underscore character, \"_\". ")
   for element in offenders:
      print("  - {}".format(element))
   print("Please use a dash, \"-\", or another character instead.")
   print("""\nTip: When removing \"_\" chars from the sample IDs, make sure that the 
         that the sample ID and parameters ID are exactly the same.

         ALLOWED - Sample ID and params ID match:
            POL_Poeo_meseres_H678-98:
              id: "H678-98"
            CTE_Horm_cali_U789089
              id: "U789089"
         NOT ALLOWED:
            POL_Poeo_meseres_H678_98: (Underscore in sampleid)
              id: "H678_98"
            CTE_Horm_cali_U789089: (Mismatch of IDs)
              id: "U7_89089"
         """)
   sys.exit()

# Parse the samplename in the configfile and split the species and id information
offenders = {}
for sample in config["samples"]:
   species = "_".join(sample.split("_")[:-1])
   sample_id = sample.split("_")[-1]
   if sample_id != config["samples"][sample]["id"]:
      offenders[sample] = config["samples"][sample]["id"]
if len(offenders) > 0:
   print("The following sample IDs in the samplename don't match the id in the sample parameters.")
   for key in offenders:
      print("  - sample: {} id: {}".format(key, offenders[key]))
   print("""Please ensure that the sample ID and parameters ID are exactly the same.

         ALLOWED - Sample ID and params ID match:
            POL_Poeo_meseres_H678-98:
              id: "H678-98"
            CTE_Horm_cali_U789089
              id: "U789089"
         NOT ALLOWED:
            POL_Poeo_meseres_H678_98: (Underscore in sampleid)
              id: "H678_98"
            CTE_Horm_cali_U789089: (Mismatch of IDs)
              id: "U7_89089"
         """)
   sys.exit()

# Now that we have parsed all of the sequence IDs we are confident that we can
#  correctly split the sample name and sample ID
# First, assert that the sample names are all max 12 characters + 4 for 3-letter code and _
offenders = {}
for sample in config["samples"]:
   realsample = "_".join(sample.split("_")[:-1])
   if len(realsample) > 16:
      offenders[sample] = realsample
if len(offenders) > 0:
   print("The following sample names are too long. Max 16 chars")
   for key in offenders:
      print("  - String:{} Sample:{}".format(key, offenders[key]))
   sys.exit()



# Now we need to change how the SS_lib_type_param is handled
for sample in config["samples"]:
   thisss_lib = config["samples"][sample]["SS_lib_type"]
   if thisss_lib in ["None", "NA", "none", "no", "nein"]:
      config["samples"][sample]["SS_lib_type"] = ""
   elif thisss_lib in ["rf", "RF", "Rf", "rF"]:
      config["samples"][sample]["SS_lib_type"] = "--SS_lib_type RF"
   elif thisss_lib in ["fr", "FR", "Fr", "fR"]:
      config["samples"][sample]["SS_lib_type"] = "--SS_lib_type FR"

#make config_lib pairs for the purpose of making reads
config["sample_lib"] = {}
for sample in config["samples"]:
   for lib in config["samples"][sample]["libs"]["short"]:
      thisd = {"read1": config["samples"][sample]["libs"]["short"][lib]["read1"],
               "read2": config["samples"][sample]["libs"]["short"][lib]["read2"]}
      config["sample_lib"]["{}_{}".format(sample, lib)] = thisd


def read_number_from_file(filename):
    with open(filename, "r") as f:
       for line in f:
           if line.strip():
               return line.strip()

def singleline_to_multiline(infile, outfile):
   """
   This function takes a fasta file and wraps it to N characters
   """
   out_handle = open(outfile, 'w') #open outfile for writing

   with open(infile, "rU") as handle:
      for record in SeqIO.parse(handle, "fasta"):
         SeqIO.write(record, out_handle, "fasta")
   out_handle.close()

#sample to
rule all:
    input:
        #make the symlinks
        expand("reads/{sample_lib}_f.fastq.gz", sample_lib = config["sample_lib"]),
        expand("reads/{sample_lib}_r.fastq.gz", sample_lib = config["sample_lib"]),
        ##trim the reads
        expand("trimmed/{sample_lib}_{readtype}.trim.fastq.gz", sample_lib = config["sample_lib"], readtype = ["f", "r"]),
        ##now assemble the transcriptomes
        expand("txomes/raw/{sample}_TRI_raw.fasta", sample = config["samples"]),
        ## make it multi-line and rename
        expand("txomes/final/{sample}_TRI.fasta", sample = config["samples"]),
#        ## transdecoder and translate, rename the pep files
#        #expand("pepfiles/final/{sample}.pep", sample = config["samples"]),
#        ## softlinks for /data/ncbi/db
#        #expand("/data/ncbi/db/{sample}.fasta",  sample = config["samples"]),
#        #expand("/data/ncbi/db/{sample}.pep",  sample = config["samples"]),
#        # make the blast databases
#        #expand("/data/ncbi/db/{sample}.fasta.nhr",  sample = config["samples"]),
#        #expand("/data/ncbi/db/{sample}.pep.phr",  sample = config["samples"]),
#        #expand("/data/ncbi/db/{sample}.dmnd",  sample = config["samples"]),
#        ## make the report
#        #expand("info/counts/raw/{sample}_raw_count.txt", sample = config["samples"]),
#        #expand("info/counts/trimmed/{sample}_trimmed_count.txt", sample = config["samples"]),
#        #expand("info/trimming/txome/{sample}_qual_rejects.txt", sample = config["samples"]),
#        #expand("info/trimming/txome/{sample}_adapter_rejects.txt", sample = config["samples"]),
#        #expand("info/trimming/txome/{sample}_gc_rejects.txt", sample = config["samples"]),
#        #expand("info/counts/fasta/{sample}_fasta_count.txt", sample = config["samples"]),
#        #expand("info/counts/peps/{sample}_pep_count.txt", sample = config["samples"]),
#        #expand("info/counts/fasta/{sample}_N50.txt", sample = config["samples"]),
#        #expand("info/counts/peps/{sample}_N50.txt", sample = config["samples"]),
#        #expand("info/counts/fasta/{sample}_meanlen.txt", sample = config["samples"]),
#        #expand("info/counts/fasta/{sample}_N90.txt", sample = config["samples"]),
#        #expand("info/counts/fasta/{sample}_mb.txt", sample = config["samples"]),
#        #"report/final_report.txt"
#
###############################################################
#  __   __   ___  __   __   __   __   ___  __   __  _       __
# |__) |__) |__  |__) |__) /  \ /  ` |__  /__` /__` | |\ | / _`
# |    |  \ |___ |    |  \ \__/ \__, |___ .__/ .__/ | | \| \__>
#
#  These are the rules for preprocessing the inputs in order
#   to prepare for the rest of the pipeline.

# preprocessing rule 1
rule make_symlinks:
   output:
       readsf = "reads/{sample_lib}_f.fastq.gz",
       readsr = "reads/{sample_lib}_r.fastq.gz"
   threads:
       1
   run:
      if not os.path.exists("reads"):
         print("Making a directory called 'reads'.", file = sys.stderr)
         os.makedirs("reads")
      else:
         print("The 'reads' directory exists already.", file = sys.stderr)

      #make symlinks for the illumina reads
      print("Making read symlinks.", file =sys.stderr)
      print("  - Checking files for sample {}".format(wildcards.sample_lib), file=sys.stderr)
      for direction in ["read1", "read2"]:
          rna_f=""
          if direction == "read1":
              rna_f = "reads/{}_f.fastq.gz".format(wildcards.sample_lib)
          elif direction == "read2":
              rna_f = "reads/{}_r.fastq.gz".format(wildcards.sample_lib)
          if not os.path.exists(rna_f):
             print("    - Making a symlink for {}".format(rna_f), file=sys.stderr)
             reads_path = config["sample_lib"][wildcards.sample_lib][direction]
             if not os.path.exists(reads_path):
                raise Exception("{} does not exist.".format(reads_path))
             os.symlink(reads_path, rna_f)
          else:
             print("    - A symlink already exists for {}".format(rna_f), file=sys.stderr)

rule rename_f:
   input:
       illumina_f = "reads/{sample_lib}_f.fastq.gz",
   output:
       illumina_f = "reads/renamed/{sample_lib}_f.renamed.fastq.gz",
   threads:
       1
   shell:
       """
       zcat {input.illumina_f} | awk '{{print (NR%4 == 1) ? "@1_" ++i "/1": $0}}' | gzip -c > {output.illumina_f}
       """

rule rename_r:
   input:
       illumina_r = "reads/{sample_lib}_r.fastq.gz",
   output:
       illumina_r = "reads/renamed/{sample_lib}_r.renamed.fastq.gz"
   threads:
       1
   shell:
       """
       zcat {input.illumina_r} | awk '{{print (NR%4 == 1) ? "@1_" ++i "/2": $0}}' | gzip -c > {output.illumina_r}
       """

#
# end of rules for preprocessing
#
#############################################################


#############################################################
#       _ ____                _
#     (_) / /_  ______ ___  (_)___  ____ _
#    / / / / / / / __ `__ \/ / __ \/ __ `/
#   / / / / /_/ / / / / / / / / / / /_/ /
#  /_/_/_/\__,_/_/ /_/ /_/_/_/ /_/\__,_/
#
#  These are the rules for working with the Illumina reads.

# illumina rule 1
rule trim_pairs:
    input:
        f = "reads/renamed/{sample_lib}_f.renamed.fastq.gz",
        r = "reads/renamed/{sample_lib}_r.renamed.fastq.gz",
        trim_jar = os.path.join(trimmomatic, "trimmomatic-0.35.jar"),
        adapter_path = "adapters/TruNexx.fa"
    output:
        f = "trimmed/{sample_lib}_f.trim.fastq.gz",
        r = "trimmed/{sample_lib}_r.trim.fastq.gz",
        u1= temp("trimmed/{sample_lib}_f.trim.unpaired.fastq.gz"),
        u2= temp("trimmed/{sample_lib}_r.trim.unpaired.fastq.gz")
    threads:
        maxthreads
    shell:
        """java -jar {input.trim_jar} PE \
        -phred33 -threads {threads} \
        {input.f} {input.r} \
        {output.f} \
        {output.u1} \
        {output.r} \
        {output.u2} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10:8:TRUE \
        SLIDINGWINDOW:10:30 MINLEN:50"""

# illumina rule 2
rule assemble_txome:
    input:
        f_paired = lambda w: ["trimmed/{}_{}_f.trim.fastq.gz".format(w.sample, x) for x in config["samples"][w.sample]["libs"]["short"]] ,
        r_paired = lambda w: ["trimmed/{}_{}_r.trim.fastq.gz".format(w.sample, x) for x in config["samples"][w.sample]["libs"]["short"]]
    output:
        assemblypath = "txomes/raw/{sample}_TRI_raw.fasta"
    params:
        outpath = "txomes/trinity_{sample}",
        outfasta = "txomes/trinity_{sample}/Trinity.fasta",
        dockeruser = dockeruser,
        sslibtype = lambda wildcards: config["samples"][wildcards.sample]["SS_lib_type"],
        freads = lambda w: ",".join(["`pwd`/trimmed/{}_{}_f.trim.fastq.gz".format(w.sample, x) for x in config["samples"][w.sample]["libs"]["short"]]),
        rreads = lambda w: ",".join(["`pwd`/trimmed/{}_{}_r.trim.fastq.gz".format(w.sample, x) for x in config["samples"][w.sample]["libs"]["short"]])
    threads:
        maxthreads
    shell:
        """
        docker run \
        -u $(id -u):$(id -g) --rm  \
        -v `pwd`:`pwd` \
        -v /etc/passwd:/etc/passwd \
        trinity_bowtie2.3.4.1 Trinity \
        --seqType fq \
        --left  {params.freads} \
        --right {params.rreads} \
        --max_memory 100G \
        --CPU {threads} \
        --trimmomatic \
        {params.sslibtype} \
        --output `pwd`/{params.outpath};
        mv {params.outfasta} {output.assemblypath}
        rm -rf {params.outpath}
        """

# illumina rule 3
rule correct_trinity_names:
    """
    The goal of this rule is to rename the fasta headers for the transcriptomes
     into a standardized and more useful string, while ensuring that the headers
     are still unique identifiers

    The headers start by looking like this:
      >TRINITY_DN19333_c0_g4_i1 len=202 path=[0:0-201]

    But they actually need to look like this:

                                 11111111112222222222333333333344444444445
                        12345678901234567890123456789012345678901234567890

    Trinity Headers    >CTE_Genus_specie|201905R_RNA4|19335c11g1i17.23
    RNAspades Headers  >CTE_Genus_specie|201905R_RNA4|592239g10470i10.23
    Space for fields        123456789111 123456789111 123456789111111111
                                     012          012          012345678
    12 chars for genus_species. Can be anything you want without | char.
    12 chars for unique_sample identifier
    """
    input:
        assem = "txomes/raw/{sample}_TRI_raw.fasta"
    output:
        fassem = "txomes/final/{sample}_TRI.fasta"
    params:
        samplename = "{sample}"
    run:
        out_handle = open(output.fassem, 'w') #open outfile for writing

        with open(infile, "rU") as handle:
           for record in SeqIO.parse(handle, "fasta"):
              record.name = ""
              record.description = ""
              recid = record.id.replace("TRINITY_DN", "")
              recid = recid.replace("_", "")
              record.id = "{}|{}|{}".format(wildcards.sample,
                        config["samples"][wildcards.sample]["id"],
                        recid)
              SeqIO.write(record, out_handle, "fasta")
              out_handle.close()


## illumina rule 4
#rule translate_txome:
#    input:
#        txome = "txomes/final/{sample}.fasta"
#    output:
#        outpeps = temp("{sample}.fasta.transdecoder_dir/longest_orfs.pep")
#    shell:
#        """
#        TransDecoder.LongOrfs -t {input.txome}
#        """
## illumina rule 5
#rule trim_transdecoder_names:
#    """
#    This trims transdecoder names from an extended string like:
#      >CHA_Caecosagitta_macrocephala_DS244_DN57546_c0_g1_i1.p1 type:internal len:122 gc:universal CHA_Caecosagitta_macrocephala_DS244_DN57546_c0_g1_i1:1-363(+)
#    to:
#      >CHA_Caecosagitta_macrocephala_DS244_DN57546_c0_g1_i1.p1
#    """
#    input:
#        pepfile = "{sample}.fasta.transdecoder_dir/longest_orfs.pep"
#    output:
#        renamed = temp("pepfiles/temp/{sample}_longest_orfs_renamed.pep")
#    shell:
#        """
#        bioawk -cfastx '{{printf(">%s\\n%s\\n", $name, $seq)}}' {input.pepfile} > {output.renamed}
#        """
#
## illumina rule 6
#rule move_and_multiline_peps:
#    """
#    This moves the translated txome to its final resting place and
#    turns it from singleline to multiline.
#    """
#    input:
#        pepfile = "pepfiles/temp/{sample}_longest_orfs_renamed.pep"
#    output:
#        pepfinal = "pepfiles/final/{sample}.pep"
#    params:
#        rmdir = lambda wildcards: "{}.fasta.transdecoder_dir".format(wildcards.sample),
#        checkpoints = lambda wildcards: "{}.fasta.transdecoder_dir.__checkpoints_longorfs".format(wildcards.sample)
#    run:
#        singleline_to_multiline(input.pepfile, output.pepfinal)
#        if os.path.exists(params.rmdir):
#           shutil.rmtree(params.rmdir)
#        if os.path.exists(params.checkpoints):
#           shutil.rmtree(params.checkpoints)
#
## illumina rule 7 fasta softlink to db
#rule softlink_fasta_data:
#    input:
#        assem = "txomes/final/{sample}.fasta"
#    output:
#        outlink = "/data/ncbi/db/{sample}.fasta"
#    params:
#        absln = lambda wildcards: "{}/txomes/final/{}.fasta".format(curr_dir, wildcards.sample)
#    shell:
#        """
#        ln -s {params.absln} {output.outlink}
#        """
#
## illumina rule 8 pep softlink to db
#rule softlink_pep_data:
#    input:
#        pep = "pepfiles/final/{sample}.pep"
#    output:
#        outlink = "/data/ncbi/db/{sample}.pep"
#    params:
#        absln = lambda wildcards: "{}/pepfiles/final/{}.pep".format(curr_dir, wildcards.sample)
#    shell:
#        """
#        ln -s {params.absln} {output.outlink}
#        """
#
## illumina rule 9
#rule nucl_db_data:
#    input:
#        inp = "/data/ncbi/db/{sample}.fasta"
#    output:
#        out = "/data/ncbi/db/{sample}.fasta.nhr"
#    shell:
#        """
#        makeblastdb -in {input.inp} -input_type fasta -dbtype nucl -parse_seqids
#        """
#
## illumina rule 10
#rule prot_db_data:
#    input:
#        inp = "/data/ncbi/db/{sample}.pep"
#    output:
#        out = "/data/ncbi/db/{sample}.pep.phr"
#    shell:
#        """
#        makeblastdb -in {input.inp} -input_type fasta -dbtype prot -parse_seqids
#        """
#
## ILLUMINA rule 11
#rule diamond_data:
#    input:
#        pepfile="/data/ncbi/db/{sample}.pep"
#    output:
#        ddb = "/data/ncbi/db/{sample}.dmnd"
#    shell:
#        """
#        diamond makedb --in {input.pepfile} -d {output.ddb}
#        """
## report
## this section handles writing a report on all of the samples.
#
## Fields to acquire:
##  - number of reads in untrimmed
##  - number of reads in trimmed
##  - number of transcripts
##  - number of peps
#
#rule number_reads_untrimmed:
#    input:
#        raw = "reads/{sample}_f.fastq.gz"
#    output:
#        raw_count = "info/counts/raw/{sample}_raw_count.txt"
#    shell:
#        """
#        echo $(bioawk -cfastx 'END{{print NR}}' {input.raw}) > {output.raw_count}
#        """
#
#rule number_reads_trimmed:
#    input:
#        trimmed = "trimmed/{sample}_f.trim.fastq.gz"
#    output:
#        trimmed_count = "info/counts/trimmed/{sample}_trimmed_count.txt"
#    shell:
#        """
#        echo $(bioawk -cfastx 'END{{print NR}}' {input.trimmed}) > {output.trimmed_count}
#        """
#
#rule number_transcripts:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        fasta_counts = "info/counts/fasta/{sample}_fasta_count.txt"
#    shell:
#        """
#        echo $(bioawk -cfastx 'END{{print NR}}' {input.fasta}) > {output.fasta_counts}
#        """
#
#rule number_peps:
#    input:
#        pep = "pepfiles/final/{sample}.pep"
#    output:
#        pep_count = "info/counts/peps/{sample}_pep_count.txt"
#    shell:
#        """
#        echo $(bioawk -cfastx 'END{{print NR}}' {input.pep}) > {output.pep_count}
#        """
#
#rule calcN50_tx:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        fasta_N50 = "info/counts/fasta/{sample}_N50.txt"
#    shell:
#        """
#        bioawk -cfastx '{{print(length($seq))}}' {input.fasta} | \
#         sort -n | \
#         awk '{{len[i++]=$1;sum+=$1}} \
#              END {{for (j=0;j<i+1;j++) {{csum+=len[j]; \
#              if (csum>=sum/2) {{print len[j];break}} }} }}' > {output.fasta_N50}
#        """
#
#rule calcN90_tx:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        fasta_N90 = "info/counts/fasta/{sample}_N90.txt"
#    shell:
#        """
#        bioawk -cfastx '{{print(length($seq))}}' {input.fasta} | \
#         sort -n | \
#         awk '{{len[i++]=$1;sum+=$1}} \
#              END {{for (j=0;j<i+1;j++) {{csum+=len[j]; \
#              if (csum>=(sum*0.9)) {{print len[j];break}} }} }}' > {output.fasta_N90}
#        """
#
#rule calcN50_pep:
#    input:
#        pep = "pepfiles/final/{sample}.pep"
#    output:
#        pep_N50 = "info/counts/peps/{sample}_N50.txt"
#    shell:
#        """
#        bioawk -cfastx '{{print(length($seq))}}' {input.pep} | \
#         sort -n | \
#         awk '{{len[i++]=$1;sum+=$1}} \
#              END {{for (j=0;j<i+1;j++) {{csum+=len[j]; \
#              if (csum>=sum/2) {{print len[j];break}} }} }}' > {output.pep_N50}
#        """
#
#rule quality_rejects:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        q_rejects = "info/trimming/txome/{sample}_qual_rejects.txt"
#    shell:
#        """
#        echo "0" > {output.q_rejects}
#        """
#
#rule adapter_rejects:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        a_rejects = "info/trimming/txome/{sample}_adapter_rejects.txt"
#    shell:
#        """
#        echo "0" > {output.a_rejects}
#        """
#
#rule GC_rejects:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        g_rejects = "info/trimming/txome/{sample}_gc_rejects.txt"
#    shell:
#        """
#        echo "0" > {output.g_rejects}
#        """
#
#rule calcmean:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        fasta_meanlen = "info/counts/fasta/{sample}_meanlen.txt"
#    shell:
#        """
#        bioawk -cfastx '{{sum += length($seq)}} \
#          END{{printf("%.2f", sum/NR)}}' {input.fasta} > {output.fasta_meanlen}
#        """
#
#rule calcMB:
#    input:
#        fasta = "txomes/final/{sample}.fasta"
#    output:
#        fasta_numMB = "info/counts/fasta/{sample}_mb.txt"
#    shell:
#        """
#        bioawk -cfastx '{{sum += length($seq)}} \
#          END{{printf("%.2f", sum/1000000)}}' {input.fasta} > {output.fasta_numMB}
#        """
#
#
#
#rule collate_report:
#    """ this method prints out the final qc report of all the libraries."""
#    input:
#        raw_num     = expand("info/counts/raw/{sample}_raw_count.txt", sample = config["samples"]),
#        qual_reject = expand("info/trimming/txome/{sample}_qual_rejects.txt", sample = config["samples"]),
#        adap_reject = expand("info/trimming/txome/{sample}_adapter_rejects.txt", sample = config["samples"]),
#        gc_reject   = expand("info/trimming/txome/{sample}_gc_rejects.txt", sample = config["samples"]),
#        trim_num    = expand("info/counts/trimmed/{sample}_trimmed_count.txt", sample = config["samples"]),
#        txome_count = expand("info/counts/fasta/{sample}_fasta_count.txt", sample = config["samples"]),
#        pep_count   = expand("info/counts/peps/{sample}_pep_count.txt", sample = config["samples"]),
#        txome_N50   = expand("info/counts/fasta/{sample}_N50.txt", sample = config["samples"]),
#        pep_N50     = expand("info/counts/peps/{sample}_N50.txt", sample = config["samples"]),
#        tx_meanlen  = expand("info/counts/fasta/{sample}_meanlen.txt", sample = config["samples"]),
#        tx_N90      = expand("info/counts/fasta/{sample}_N90.txt", sample = config["samples"]),
#        tx_MB       = expand("info/counts/fasta/{sample}_mb.txt", sample = config["samples"]),
#
#    output:
#        "report/final_report.txt"
#    run:
#        finalout_handle = open(output[0], "w")
#        #print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
#        print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
#                        "sample",
#                        "trim_efficiency",
#                        "num_reads",
#                        "num_reads_trimmed",
#                        "quality_rejects",
#                        "adapter_rejects",
#                        "gc_rejects",
#                        "transcripts_count",
#                        "transcripts_meanlen",
#                        "transcripts_N50",
#                        "transcripts_N90",
#                        "transcripts_MB",
#                        "peps_count",
#                        "peps_N50"),
#              file = finalout_handle)
#        for thiss in sorted(config["samples"]):
#            raw_num     = read_number_from_file("info/counts/raw/{}_raw_count.txt".format(thiss))
#            trim_num    = read_number_from_file("info/counts/trimmed/{}_trimmed_count.txt".format(thiss))
#            txome_count = read_number_from_file("info/counts/fasta/{}_fasta_count.txt".format(thiss))
#            pep_count   = read_number_from_file("info/counts/peps/{}_pep_count.txt".format(thiss))
#            txome_N50   = read_number_from_file("info/counts/fasta/{}_N50.txt".format(thiss))
#            pep_N50     = read_number_from_file("info/counts/peps/{}_N50.txt".format(thiss))
#            qual_reject = read_number_from_file("info/trimming/txome/{}_qual_rejects.txt".format(thiss))
#            adap_reject = read_number_from_file("info/trimming/txome/{}_adapter_rejects.txt".format(thiss))
#            gc_reject   = read_number_from_file("info/trimming/txome/{}_gc_rejects.txt".format(thiss))
#            tx_meanlen  = read_number_from_file("info/counts/fasta/{}_meanlen.txt".format(thiss))
#            txome_N90   = read_number_from_file("info/counts/fasta/{}_N90.txt".format(thiss))
#            txome_MB    = read_number_from_file("info/counts/fasta/{}_mb.txt".format(thiss))
#            print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
#               thiss,
#               int(trim_num)/int(raw_num),
#               raw_num,
#               trim_num,
#               qual_reject,
#               adap_reject,
#               gc_reject,
#               txome_count,
#               tx_meanlen,
#               txome_N50,
#               txome_N90,
#               txome_MB,
#               pep_count,
#               pep_N50
#            ), file = finalout_handle)
#        finalout_handle.close()
