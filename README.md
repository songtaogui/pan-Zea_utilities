# pan-Zea-utilities

- [pan-Zea-utilities](#pan-zea-utilities)
  - [Prerequisites](#prerequisites)
    - [Runtime environment:](#runtime-environment)
    - [Dependencies](#dependencies)
  - [Detailed usage](#detailed-usage)
    - [Genetic](#genetic)
      - [PANZ_SVflankSNP_LD.sh and SV_LD_type_draw.r](#panz_svflanksnp_ldsh-and-sv_ld_type_drawr)
      - [PANZ_QTL_FineMap.sh](#panz_qtl_finemapsh)
      - [PANZ_part_h2.sh](#panz_part_h2sh)
      - [PANZ_MAGMA.sh](#panz_magmash)
    - [Genomic](#genomic)
      - [PANZ_determine_core_dispensable.sh](#panz_determine_core_dispensablesh)
      - [PANZ_SubG_PAV.sh](#panz_subg_pavsh)
      - [PANZ_gene_age.sh](#panz_gene_agesh)
      - [PANZ_SV_Annotation.sh](#panz_sv_annotationsh)
    - [Statistic](#statistic)
      - [PANZ_freq_enrich.sh](#panz_freq_enrichsh)
      - [PANZ_rankINT.sh](#panz_rankintsh)
      - [PANZ_regional_enrich.sh](#panz_regional_enrichsh)
  - [Citations](#citations)

This repository collected the miscellaneous analysis scripts used in the research of pan-Zea genome and genetics.

Class | Script_name | Descriptions
------- | ------- | -------
Genetic | `PANZ_SVflankSNP_LD.sh` | Investigate the LD level of given query variants to nearby SNPs
Genetic | `SV_LD_type_draw.r` | Companion of `PANZ_SVflankSNP_LD.sh`. Can also be used for visualization.
Genetic | `PANZ_QTL_FineMap.sh` | Identify QTLs from GWAS results and perform statistical fine mapping
Genetic | `PANZ_part_h2.sh` | Partition genetic variants and calculate h2 for each part.
Genetic | `PANZ_MAGMA.sh` | Perform regional association analysis of genic regions
Genomic | `PANZ_determine_core_dispensable.sh` | Determine if a gene is core or dispensable based on gene presence and absence matrix
Genomic | `PANZ_SubG_PAV.sh` | Identify subgroup unblanced PAV genes
Genomic | `PANZ_gene_age.sh` | Calculate gene age based on sequence similarity within the species tree
Genomic | `PANZ_SV_Annotation.sh` | Annotate SVs relative to nearby genes and TE classes
Statistic | `PANZ_freq_enrich.sh` | Enrichment analysis given frequency stat query and ref
Statistic | `PANZ_rankINT.sh` | Normalizing input with RankINT
Statistic | `PANZ_regional_enrich.sh` | Genome regional enrichment analysis


## Prerequisites

### Runtime environment:
>
> - Linux, tested with version 3.10.0-862.el7.x86_64 (Red Hat 4.8.5-28)
> - bash, tested with version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
> - perl 5, tested with v5.30.1

Most of the pipelines were wrote with bash, R, and Perl, and tested on Linux platform. However, other unix-like platform could also work. The only thing to note is that some of the pipelines used the `bash` redirection features that not supported in basic `sh`:

```sh
# in bash or zsh it works
bash -c 'cat <(echo "hello")'
# hello

# in sh it failed
sh -c 'cat <(echo "hello")'
# sh: 1: Syntax error: "(" unexpected
```

Thus, try not to run the pipelines through `sh PANZ_some_pipe.sh`.

### Dependencies

Most of the pipelines would require the presence of other third-party softwares, such as [GNU parallel](https://www.gnu.org/software/parallel/). The detailed dependencies would be listed in the usage of certain pipeline.

## Detailed usage

### Genetic

#### PANZ_SVflankSNP_LD.sh and SV_LD_type_draw.r

These script could do the similar analysis on the representation of the query genetic variant with nearby SNPs, as described in:

>Stuart T, Eichten S R, Cahn J, et al. [Population scale mapping of transposable element diversity reveals links to gene regulation and epigenomic variation[J]](https://elifesciences.org/articles/20777). elife, 2016, 5: e20777.

- PANZ_SVflankSNP_LD.sh

```
------------------------------------------------------------
Calculate SV LDs with nearby SNP/InDels use method in
(Stuart et al. eLife 2016;5:e20777. DOI: 10.7554/eLife.20777)
------------------------------------------------------------
Dependency: bcftools vcftools perl Rscript csvtk GNU-parallel
------------------------------------------------------------
USAGE:
    bash PANZ_SVflankSNP_LD.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -t, --threads    <num>    [O]    set threads (default 2)
    -q, --query      <str>    [R]    Input query vcf file, bcf or gz and indexed (SV)
    -r, --ref        <str>    [R]    Input ref vcf file, bcf or gz and indexed (SNP/Indel)
    -s, --size       <num>    [O]    Number of flanking ref records for each query (Default: 150)
    -p, --prefix     <str>    [O]    Prefix for outputs (Default: Calc_LD_OUT)
    -R, --rscript    <str>    [O]    PATH to the 'SV_LD_type_draw.r' script (Default: ./SV_LD_type_draw.r)
    -m, --mode       <str>    [O]
                    Analysis mode:{ 0; 1; 2 } (Default: 0)
                        0 ---> only out put LD type
                        1 ---> only draw heatmap and linesfig
                        2 ---> out LD type and draw figs
                    use 1,2 on all LD files is SLOW, you can run the stand-alone rscript on each LD file later.
    --clean                   [O]    Clean each variant LD result, just keep the main summary output, useful when there
                                    are too many variants (if disk file number limit is an issue), and you donot need
                                    each result file to draw the ld heatmap figs.
Note:
    needs these software in PATH to run the pipeline:
    GNU parallel; bcftools; Rscript; vcftools; csvtk;
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

- SV_LD_type_draw.r

```sh
--------------------------------------------------
Calculate sv_snp-snp_mid num and draw.
--------------------------------------------------
Usage:
        Rscript <Program> <mode> <vcftools.LD> <size> <outdir>
--------------------------------------------------
OPTIONs:
        mode:
                0 ---> only out put LD type
                1 ---> only draw heatmap and linesfig
                2 ---> out LD type and draw figs
        haploview.LD:
                LD file generated from the main program
        size:
                flanking ref size used in the main program
    outdir:
        path of the outputs
--------------------------------------------------
```

Outputs:

If you ran `PANZ_SVflankSNP_LD.sh` with `--mode 2`, you would get 4 outputs:

1. output LD file records the pairwise LD of the variants flanking the query. The position of the variants were modified as the order of the query and flanking variants:
    ```sh
        # e.g. if you set the flanking of 150 SNPs, that means the resulting LD files would be the pairwise LD values of 301 items (left 150 + query + right 150):
        CHR   POS1   POS2   N_INDV   R^2
        1     1      2      662      0.000105541
        1     1      3      661      5.96288e-05
        1     1      4      663      6.86625e-06
        1     1      5      662      6.88707e-06
        1     1      6      664      6.84552e-06
        ...
        ...
    ```

2. the LD rank summary output of the query:
    
```sh
    # format <query_ID>   <# query-SNP ranks over nearby SNP-SNP>   <# Total nearby SNPs>  <LD-leve tag>
    Zm00001d027244_PAV  293     300     high
```

3. and 4. visualisations in PDF format of the detailed LD heatmap and the LD rank comparison line plots

NOTE:
if you have ran the pipeline with `--mode 0` or `--mode 1`, and you would like to draw plots for specific query, you could just take the output LD file of that query to `SV_LD_type_draw.r` :

```sh
SV_LD_type_draw.r 2 myQuery_flank150_hplv.geno.ld 150 .
# Outputs:
# myQuery  293     300     high
# heatmap Fig: ./myQuery_heatmap.pdf
# lines Fig: ./myQuery_lines.pdf
```

#### PANZ_QTL_FineMap.sh

```sh
------------------------------------------------------------
Identify QTLs and perform statistical fine mapping for PANZ 
------------------------------------------------------------
Dependence:csvtk perl HARVESTER DAP-G
------------------------------------------------------------
USAGE:
    

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -t, --threads       <num>   [O]     set threads (default: 2)

# I/O:
    -o, --out           <str>   [O]     Output prefix, will create a dir with this string in the current path,
                                    so make sure the dir does not exist (default: PANZ_FineMap_Out)
    -v, --gvcf          <str>   [R]     input gvcf file (should be indexed).
    --trait             <str>   [R]     input trait file name, make sure the file was named as "trait_name.prefix",
                                    File Format (TSV format, the key word "<trait>" is needed):
                                        <trait>    trait_name
    --covar             <str>   [R]     Input quantitative covariates from a plain text file in format (TSV):
                                        <covar>    covar_ID
                                        sample1        1.0
                                        sample2        2.5
    --gwas              <str>   [R]     gwas file.tsv[.gz] with header( header could be any string),
                                    should contain columns of:
                                        variant_ID, chr, start, end, p-value, effect size
    --info_cols         <str>   [R]     a set of numbers seperated by comma to indicate the column of
                                    [variant_ID, chr, start, end, p-value, effect size] in the gwas file, for example,
                                    if your gwas file looks like this:
                                        Chr    start0    End    variant_ID    p-value    beta
                                        1       1234    1235    SNP01         1E-7      0.22
                                    you should set this option as "4,1,2,3,5,6"

# gwas peak calling
    --min_peak_logP     <int>   [O]     the minimum -log(P-value) of the peak to keep. (default: 5)
    --inlimit           <int>   [O]     the minimum p-value of the variant to include in a peak region. (default: 0.001)
    --flank             <int>   [O]     flanking length (in bp) of peak to include as final peak region. (default: 10000)
    --min_count         <int>   [O]     the minimum No. of variants that passed the --inlimit cutoff within a peak. (default: 3)

# fine mapping (Bayes)
    --ld_control        <0-1>   [O]     r^2 threshold within a signal cluster (default: 0.25)
    --credible_prob     <0-1>   [O]     Prob of credible set (default: 0.95)

# Miscellaneous
    --tissue            <str>   [O]     tissue id of the trait that used in the final annotation vcf file.(default: kernel)
    --tar                       [O]     tar intermediate dirs if set. (default: do not tar intermediate dirs)

------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

#### PANZ_part_h2.sh

```sh
------------------------------------------------------------
Identify QTLs and perform statistical fine mapping for PANZ
------------------------------------------------------------
Dependence:csvtk perl HARVESTER DAP-G
------------------------------------------------------------
USAGE:


OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -t, --threads       <num>   [O]     set threads (default: 2)

# I/O:
    -o, --out           <str>   [O]     Output prefix, will create a dir with this string in the current path,
                                    so make sure the dir does not exist (default: PANZ_FineMap_Out)
    -v, --gvcf          <str>   [R]     input gvcf file (should be indexed).
    --trait             <str>   [R]     input trait file name, make sure the file was named as "trait_name.prefix",
                                    File Format (TSV format, the key word "<trait>" is needed):
                                        <trait>    trait_name
    --covar             <str>   [R]     Input quantitative covariates from a plain text file in format (TSV):
                                        <covar>    covar_ID
                                        sample1        1.0
                                        sample2        2.5
    --gwas              <str>   [R]     gwas file.tsv[.gz] with header( header could be any string),
                                    should contain columns of:
                                        variant_ID, chr, start, end, p-value, effect size
    --info_cols         <str>   [R]     a set of numbers seperated by comma to indicate the column of
                                    [variant_ID, chr, start, end, p-value, effect size] in the gwas file, for example,
                                    if your gwas file looks like this:
                                        Chr    start0    End    variant_ID    p-value    beta
                                        1       1234    1235    SNP01         1E-7      0.22
                                    you should set this option as "4,1,2,3,5,6"

# gwas peak calling
    --min_peak_logP     <int>   [O]     the minimum -log(P-value) of the peak to keep. (default: 5)
    --inlimit           <int>   [O]     the minimum p-value of the variant to include in a peak region. (default: 0.001)
    --flank             <int>   [O]     flanking length (in bp) of peak to include as final peak region. (default: 10000)
    --min_count         <int>   [O]     the minimum No. of variants that passed the --inlimit cutoff within a peak. (default: 3)

# fine mapping (Bayes)
    --ld_control        <0-1>   [O]     r^2 threshold within a signal cluster (default: 0.25)
    --credible_prob     <0-1>   [O]     Prob of credible set (default: 0.95)

# Miscellaneous
    --tissue            <str>   [O]     tissue id of the trait that used in the final annotation vcf file.(default: kernel)
    --tar                       [O]     tar intermediate dirs if set. (default: do not tar intermediate dirs)

------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

#### PANZ_MAGMA.sh

```sh
------------------------------------------------------------
Perform regional association analysis of genic regions for PANZ.
A wrapper of MAGMA
------------------------------------------------------------
Dependence: MAGMA plink
------------------------------------------------------------
USAGE:
    bash PANZ_MAGMA.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -r, --ref           <str>   [R]     Plink bed prefix of all variants
    --annote            <str>   [R]     Gene-Variant Annotation file in format:
                                        <Gene_ID>    <Variants_sep_by_space>
                                        Gene1        SNP1 SV2 INDEL3
                                        ...          ...
    --gene_sets         <str>   [R]     Gene set file in format:
                                        <Gene_SET_ID>    <Genes_sep_by_space>
                                        Gene_family_1    Gene1 Gene2 Gene3
                                        ...              ...
    --pheno             <str>   [R]     Phenotype in plink format
                                    *** NOTE:   Phenotype should not contain duplicate samples,
                                                and the first Phenotype should not be "NA"
    --covar             <str>   [R]     covariant in plink format
    --region            <str>   [O]     Physical region included for the analysis,
                                        in format "CHR:START-END", use all variants
                                        if not provided.
    -o, --out           <str>   [O]     Output prefix, will create a dir accordingly for the outputs (default: PANZ_Gene_GWAS_out)

------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

### Genomic

#### PANZ_determine_core_dispensable.sh

```sh
------------------------------------------------------------
Determine core and dispensable genes/gene families
using pvalue of binomial test of gene loss rate.
------------------------------------------------------------
Dependency: parallel perl Rscript
------------------------------------------------------------
USAGE:
    bash PANZ_determine_core_dispensable.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -m, --matrix     <str>    [R]    path of PAV matrix file.
                            format (tab-separated):
                            ID     Sample1    Sample2    Sample3
                            Gene1  1          0          1
                            Gene2  1          1          1
                            Gene3  0          0          1
    -r, --lossrate   <0-1>    [O]    min loss rate to determine a gene dispensable (Default: 0.01)
    -p,--pvalue      <0-1>    [O]    max pvalue cutoff to determine a gene dispensable (Default: 0.05)
    -o,--output      <str>    [O]    output file (default: core_disp_out.tsv)
    --prefix         <str>    [O]    prefix tag of output gene types (default: PANZ)
    -t, --threads    <num>    [O]    set threads (default: 2)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

#### PANZ_SubG_PAV.sh

```sh
------------------------------------------------------------
Determine subgroup unbalanced genes/gene families
using adjusted-pvalue (q-value, Storey's Method) of two-sided
Fisher's exact test of gene PAV matrix among different subgroups.
------------------------------------------------------------
Dependencies: GNU-parallel perl Rscript csvtk (All should be in PATH)
------------------------------------------------------------
USAGE:
    bash PANZ_SubG_PAV.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -m, --matrix     <str>    [R]    path of PAV matrix file.
                            format (tab-separated):
                            ID     Sample1    Sample2    Sample3
                            Gene1  1          0          1
                            Gene2  1          1          1
                            Gene3  0          0          1
    -s,--subgroup    <file>   [R]    a list of samples in subgroup, one per line.
                            Note: make sure that all samples were in the provided matrix. The comparasion
                            will be performed between the whole matrix provided and the subgroup matrix
                            extracted based on this subgroup name list.
    -p,--prefix      <str>    [O]    prefix tag of outputs, you may use subgroup name (default: PANZ)
    -f,--FDR         <0-1>    [O]    FDR cutoff to determine a gene unbalance (Default: 0.05)
    --localfdr                [O]    Use local FDR values (default: Global FDR, aka qvalue)
    -o,--output      <str>    [O]    output file (default: PAV_subgroup_unbalanced_out.tsv)
    -t, --threads    <num>    [O]    set threads (default: 2)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com

```

#### PANZ_gene_age.sh

```sh
------------------------------------------------------------
Infering Gene ages using methods discribed in:
    Zhang, Y. E. et al. Chromosomal redistribution of male-biased
    genes in mammalian evolution with two bursts of gene gain on
    the X chromosome. PLoS Biol. 8, e1000494 (2010).
------------------------------------------------------------
Dependency in PATH:
    diamond, csvtk, taxonkit, perl, GNU-parallel, blastdbcmd
------------------------------------------------------------
USAGE:
    bash PANZ_gene_age.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -d, --database   <str>    [O]    Path to NCBI-NR protein database (Default: /mnt/d/Works/NnRAD/mirror_1111/ref/NR/nr/nr)
    -i, --input      <str>    [R]    Input protein query fasta sequences
    -o,--output      <str>    [O]    Output prefix (Default: PANZ_Gene_ages)
    --taxnomy        <str>    [R]    Tab-seperated Hierarchy of NCBI Taxnomy IDs. Format:
                            <Label> <Taxnomy ID> <Hierarchy Order>
                            Example:
                                For taxnomy like:
                                        |- cellular organisms (131567)
                                            |--- Eukaryota (2759)
                                                |--- Alveolata (33630)
                                The taxnomy order file would be:
                                        P1_Cellular    131567    1
                                        P2_Eukaryota   2759      2
                                        P3_Alveolata   33630     3
                            ** 1. Taxnomy with same Label would be merged in the output.
                            ** 2. Please make sure all the taxnomy IDs are in SAME taxnomy tree branch
                            ** 3.  Hierarchy order should be increased by taxnomy age from old to new
    --evalue         <num>    [O]    Evalue cut-off of the blastp output (Default 1e-5)
    --identity       <0-1>    [O]    Identity cut-off of the blastp output (Default 0.3)
    -t, --threads    <num>    [O]    set threads (default: 2)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

#### PANZ_SV_Annotation.sh

```sh
------------------------------------------------------------
PANZ annotate SV with TE and Gene gff
------------------------------------------------------------
Dependence in PATH: parallel, RepeatMasker, bedtools, bcftools
------------------------------------------------------------
USAGE:
    bash PANZ_SV_Annotation.sh [OPTIIONS]
OPTIONS:
    1. SV.vcf
    2. gene.gff
    3. TE.bed: chr start end ID class strand
    4. TE.fa: NR TE sequences
    5. Flanking_length(bp)
    6. Threads
    7. output_name
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

### Statistic

#### PANZ_freq_enrich.sh

```sh
------------------------------------------------------------
PANZ freq enrichment analysis: input query and ref frequency file, and do enrichment analysis.

frequency file format:(no header)
    <Item_ID>  <Freq>
    ITEM1      123
    ITME2      888
    ...        ...

------------------------------------------------------------
USAGE:
    bash PANZ_freq_enrich.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -q, --query      <str>    [R]    Query frequency file
    -r, --ref        <str>    [R]    Referance frequency filen
    -c, --cutoff     <0-1>    [O]    Qvaule cutoff (default: 0.05)
    -o, --output     <str>    [O]    Output file name (default: PANZ_KEGG_Enrich.tsv)
    -t, --threads    <num>    [O]    set threads (default: 2)

Output format:
    <Item>  <Type>  <ratio_in_query>  <ratio_in_ref>  <raw_Pvalue>  <Qvalue>  <Lfdr>
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

#### PANZ_rankINT.sh

```sh
------------------------------------------------------------
PANZ Rank-Based Inverse Normal Transformation
------------------------------------------------------------
Dependence: Rscript
------------------------------------------------------------
USAGE:
    bash PANZ_rankINT.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -t, --trait         <str>   [R]     trait file in tsv format
    -c, --col           <int>   [O]     trait value column (default: 2)
    --na                <str>   [O]     NA represent string (default: NA)
    -o, --out           <str>   [O]     Output prefix (default: header of trait column)
    --SW_test                   [O]     Out put a Shapiro-Wilk norm test result to <out>.sw_test
    --adj               <0-1>   [O]     Adjust parameter for Rank-Based Inverse Normal Transformation (default: 0.5)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com

```

#### PANZ_regional_enrich.sh

```sh
------------------------------------------------------------
PANZ_Regional_enrich: a wrapper of R/regioneR function, with
some DIY options.
------------------------------------------------------------
Dependence:
    R/regioneR package (http://bioconductor.org/packages/release/bioc/html/regioneR.html)
------------------------------------------------------------
USAGE:
    bash PANZ_regional_enrich.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -o, --out <str>            [R]      Output prefix.
    -g, --genome <str>         [R]      Genome range in bed format. Used as the bondaries of region manipulating. Eg:
                                            chr1    0   1234567
                                            chr2    0   1345678

    -a, --query <str>          [R]      Query region file in bed format, will do permutation based on this file.

    -b, --feature <str>        [R]      Feature region file in bed format, will count overlaps of each permutation of query region with this file to evaluate the enrichment.

    -r, --ref <str>            [O]      Set the range of query region permutation. Three type of parameters are supported:
                1. [bed:your_bed_file.bed] --> Use a region file as the permutation boundary, so make sure the query region file is subset of the input region file.
                2. [flank:<int>] --> use a flanking of <int> bp length of query region (both side) as the permutation boundary, eg: "flank:10000" will flank left 10kb and right 10kb of query region.
                3. [time:<int>] --> use a flanking of <int> * <length of each query region> bp of query region as the permutation boundary, eg:
                    if your query region is:
                    #    chr1    1000    1100
                    #    chr2    5010    5020
                    and you set ref as "--ref time:5", the permutation boundary would be:
                    #    chr1    500     1600    (query length=100,flanking=100*5)
                    #    chr2    4960    5070    (query length=10,flanking=10*5)
                The default behavior of --ref (if unset) is to use "--genome" file as boundary, that is permutation along the whole genome.

    -n, --ntimes <int>         [O]      Number of permutation times.A large number of permutations will produce more accurate results and a nicer-looking plot but a permutation test can be computationally expensive.(default: 100)

    --cutoff <0-1>             [O]      P-value cutoff to call the result as significantly non-random. (default: 0.05)

    --seed <int>               [O]      Set random seeds to create reproducible results. (default: 1234)

    --plot                     [O]      Plot the permutation results.

    --force_save               [O]      Force save the permutation result rdata. The default behavior is to save only the result that passed the P-value cutoff.
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

## Citations

If you use these pipelines in your work, or you would like to know more details about them, please refer to:

<!-- > Gui, S. (2021). TITLE HERE.
> *Journal HERE*, **34**:3094-3100.- [doi:DOIhere][doi] -->
