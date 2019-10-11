# Renaming sequences

## TL/DR; The header format

```
                             11111111112222222222333333333344444444445 
                    12345678901234567890123456789012345678901234567890
Trinity Headers    >CTE_Genus_specie|201905R_RNA4|19335c11g1i17.23
RNAspades Headers  >CTE_Genus_specie|201905R_RNA4|592239g10470i10.23
Space for fields        123456789111 123456789111 123456789111111111
                                 012          012          012345678

12 chars for genus_species
12 chars for unique_sample identifier
```

## Working through the data to come to a file format

BATs renames the default Trinity fasts headers to make them more useful when analyzing individual sequences.

The header limit for some software is 50 characters. Trinity outputs sequence IDs like this: `DN10001_c0_g1_i1`. This is 16 characters, however a quick analysis of an entire transcriptome shows that the longest ID that we will likely encounter is 18 characters long.

```
bioawk -cfastx '{print(length(substr($name, 4)))}' CTE_hcal_onlynames_short.fa | sort | uniq -c | column -t | sort -k2 -n
3      12
60     13
556    14
7438   15
97578  16
8927   17
128    18
```

The information in the Trinity names can be compressed. `DN10001_c0_g1_i1` can become `10001c0g1i1` without losing information. If the Trinity headers are renamed, then the longest sequence that one is likely to find becomes 13 characters. This pattern held true even when looking at a transcriptome with over 400,000 entries.

When translating the transcriptome, `TransDecoder` adds a `.p1`, `.p2`, et cetera to the header to indicate the nth translation product from any particular transcript. A look at a transcriptome shows that there are occassionally double-digt numbers of proteins per transcript. So, the information `.p23` can be compressed to `.23` without losing information, and at most we add 3 characters when translating the transcripts in a peptide file.

## First look at fasta header

The fasta header format is starting to take shape. We now know that headers will look something like this, leaving 5 characters for buffer space for unexpected corner cases in Trinity, RNAspades, et cetera:

```
            11111111112222222222333333333344444444445 
   12345678901234567890123456789012345678901234567890
  >                             19335c11g1i17.23
```

It looks like there are now 29 characters to use to make the sequence attributable to a specific species and sample. 

## Sample provenance information.

It is important to have some information to track the sample provenance, such as an accession number, or a collection date plus individual number. Below are some possible provenance tracking numbers that we may use internally.

```
char count
         11111111112
12345678901234567890
----------------------
P914-T1
SRR7533167
R785-G11
GLO064
X8273-UY28
20020612-BW2
10158X12
2006BW29SEPT
26MAR19-BW3-19
M8734-827
201905RNA_RNA4
UM27UP-5
IOI-TY7
20190417-4
```

From this naive analysis, it looks like the longest practical sample identifier name is 14 characters long. Adding this information, along with `|` pipe delimiters, and a 3-letter clade code (CTE, CNI, VER) to our headers will look something like this: 

```
            11111111112222222222333333333344444444445 
   12345678901234567890123456789012345678901234567890
  >CTE_         |201905RNA_RNA4|19335c11g1i17.23
```

This only leaves 9 characters to define the genus and species. Perhaps we can afford to use some of the 5 buffer characters by determining the compressability of the RNAspades transcript output?

```
bioawk -cfastx '{print($name)}' RNAspades.fasta | head
NODE_1_length_25573_cov_14.379691_g0_i0
NODE_2_length_15939_cov_16.138120_g0_i1
NODE_3_length_14863_cov_16.250135_g0_i2
NODE_4_length_13985_cov_16.931339_g0_i3
NODE_5_length_13727_cov_50.846199_g1_i0
NODE_6_length_13727_cov_49.417178_g1_i1
NODE_7_length_13727_cov_49.283406_g1_i2
NODE_8_length_13712_cov_77.985071_g2_i0
```

Is everything after the period a unique identifier?

```
bioawk -cfastx '{print($name)}' RNAspades.fasta | awk -F '.' '{print($2)}' | sort | uniq -c | column -t | awk '{print($1)}' | sort | uniq -c
 357629 1
```

Yes, everything after the period in the RNAspades output headers is a unique identifier. Now, what is the length range of these unique identifiers?

```
bioawk -cfastx '{print($name)}' RNAspades.fasta | awk -F '.' '{print(length($2))}' | sort | uniq -c
     28 12
    213 13
   1932 14
  16201 15
 118899 16
 220356 17
```

The longest header length is 17 characters long. Some examples are:

```
bioawk -cfastx '{print($name)}' RNAspades.fasta | awk -F '.' '{if (length($2) == 17 ) { print($2)}} ' | head
592239_g10470_i10
251118_g17053_i10
219345_g10954_i10
497678_g17911_i10
361366_g17911_i11
```

This information can be compressed from `592239_g10470_i10` to `592239g10470i10`. This is 15 characters long. Given that we now know the length of both the RNAspades and Trinity unique sequence IDs, we can predict the length of the longest headers for each.

```
                             11111111112222222222333333333344444444445 
                    12345678901234567890123456789012345678901234567890
Trinity Headers    >CTE_          |201905RNA_RNA4|19335c11g1i17.23
RNAspades Headers  >CTE_          |201905RNA_RNA4|592239g10470i10.23
Space for Gen.sp.       1234567891
                                 0
```

In the above example more space was added to the genus/species field, and now even the longest of transcripts found in RNAspades have two additional buffer characters. Currently, there are ten spaces to put genus and species information. Some possible names for species like _Beroe forskalii_ could be:

```
         1
1234567890
B_forskali
Be_forskal
Ber_forska
Bero_forsk
Beroe_fors
```

Is this useful for species that have similar names?

```
Paraphyllina
Periphyllopsis
Periphylla
         1
1234567890
Periphylli
Periphyllo
Periphylla
```

It could still be confusing. More characters for the Genus-species delimiter would be useful. Changing the 14 chars of sample ID to 12 chars would be useful to make more space for the genus-species field.

```
                             11111111112222222222333333333344444444445 
                    12345678901234567890123456789012345678901234567890
Trinity Headers    >CTE_            |201905R_RNA4|19335c11g1i17.23
RNAspades Headers  >CTE_            |201905R_RNA4|592239g10470i10.23
Space for Gen.sp.       123456789111
                                 012
```

This is the final header format that will be useful for both Trinity assemblies and RNAspades assemblies.

