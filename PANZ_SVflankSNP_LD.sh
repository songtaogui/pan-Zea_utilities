#!/usr/bin/env bash

# set -euo pipefail
# >>>>>>>>>>>>>>>>>>>>>>>> common functions >>>>>>>>>>>>>>>>>>>>>>>>
function gst_log () {
    local info=$1
    echo -e "\033[36m[$(date +'%y-%m-%d %H:%M')]\033[0m $info" >&2
}
function gst_err () {
    local info=$1
    echo -e "\033[31m\033[7m[ERROR]\033[0m --> $info" >&2
}
function gst_warn () {
    local info=$1
    echo -e "\033[35m[WARNING]\033[0m --> $info" >&2
}

# usage check_files_exists
function check_files_exists(){
    local num_related_file=1
    local related_file=""
    for related_file in  "$@"
    do
        if [[ ! -s "$related_file" ]]; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> No file: $related_file $" >&2
            let num_related_file++
        fi
    done
    [ "$num_related_file" -ne 1 ] && exit 1
}

function check_var_empty () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        case ${var} in
            '')
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name is empty: '$var' " >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}
function check_sftw_path(){
    local num_tp_program=1
    local tp_program=""
    for tp_program in "$@"
    do
        if ! which $tp_program >/dev/null 2>&1 ; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> Program not in PATH: $tp_program $" >&2
            let num_tp_program++
        fi
    done
    [ "$num_tp_program" -ne 1 ] && exit 1
}

export -f check_sftw_path check_files_exists gst_log gst_err gst_warn
# <<<<<<<<<<<<<<<<<<<<<<<< common functions <<<<<<<<<<<<<<<<<<<<<<<<

usage="
------------------------------------------------------------
Calculate SV LDs with nearby SNP/InDels use method in
(Stuart et al. eLife 2016;5:e20777. DOI: 10.7554/eLife.20777)
------------------------------------------------------------
Dependency: bcftools vcftools perl Rscript csvtk GNU-parallel
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -t, --threads    <num>    [O]    set threads (default 2)
    -q, --query      <str>    [R]    Input query vcf file, bcf or gz and indexed (SV)
    -r, --ref        <str>    [R]    Input ref vcf file, bcf or gz and indexed (SNP/Indel)
    -s, --size       <num>    [O]    Number of flanking ref records for each query (Default: 150)
    -p, --prefix     <str>    [O]    Prefix for outputs (Default: Calc_LD_OUT)
    -R, --rscript    <str>    [O]    PATH to the 'SV_LD_type_draw.r' script (Default: $(dirname $0)/SV_LD_type_draw.r)
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
"
if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export threads=2
export query=
export ref=
export size=150
export pre=Calc_LD_OUT
export rscript="$(dirname $0)/SV_LD_type_draw.r"
export mode=0
export clean="FALSE"
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -t|--threads)
            #echo "set argument \"$1\" with value: $2" >&2
            threads=$2
            shift 2
        ;;
        -q|--query)
            query=$2
            shift 2
        ;;
        -r|--ref)
            ref=$2
            shift 2
        ;;
        -s|--size)
            size=$2
            shift 2
        ;;
        -v|--vcftools)
            vcftools_path=$2
            shift 2
        ;;
        -p|--prefix)
            pre=$2
            shift 2
        ;;
        -R|--rscript)
            rscript=$2
            shift 2
        ;;
        -m|--mode)
            mode=$2
            shift 2
        ;;
        --clean)
            clean="TRUE"
            shift 1
        ;;
        *) # unknown flag/switch
            UNKOWN_ARGS+=("$1")
            shift
        ;;
    esac
done
if [ "${#UNKOWN_ARGS[@]}" -gt 0 ];then
    echo "[ERROR] --> Wrong options: \"${UNKOWN_ARGS[@]}\"" >&2
    exit 1
fi
# check related args:
echo "$size" | [ -n "`sed -n '/^[0-9][0-9]*$/p'`" ] || gst_log "-s,--size: $size is not numeric !"
echo "$mode" | [ -n "`sed -n '/^[0-9][0-9]*$/p'`" ] || gst_log "-m,--mode: $mode is not numeric !"
check_files_exists $query $ref $rscript
check_var_empty query ref
unset UNKOWN_ARGS # restore UNKOWN_ARGS params
check_sftw_path bcftools vcftools perl Rscript csvtk parallel
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

gst_log "Start running ...
------------------------------------------------------------
            Program start with Parameters:
------------------------------------------------------------
    query: $query
    ref  : $ref
    size : $size
    CPUs : $threads
    Mode : $mode
------------------------------------------------------------
Will generate these outputs:

Total_records   : ${pre}_Total_records.tsv
Each_record_dir :
    ${pre}/LDs   : dir of each query LDs
    ${pre}/Rcalc : dir of each query types (and figs if mode not 0)
------------------------------------------------------------
"
# >>>>>>>>>>>>>>>>>>>>>>>> Checking >>>>>>>>>>>>>>>>>>>>>>>>
md5_query_sample=$( bcftools view -h $query | tail -n 1 | md5sum | cut -d " " -f 1 )
md5_ref_sample=$( bcftools view -h $ref | tail -n 1 | md5sum | cut -d " " -f 1)

if [ "$md5_query_sample" != "$md5_ref_sample" ]; then
    gst_err "Query and Ref file should have same sample orders."
    exit 1
fi

# ? mk output dirs
mkdir -p $pre &&\
mkdir -p ${pre}/tmp &&\
mkdir -p ${pre}/LDs &&\
mkdir -p ${pre}/Rcalc
if [[ $? != 0 ]]; then
    gst_err "Cannot creat output dir: $pre"
    exit 1
fi
bcftools view -h $ref -o ${pre}/vcf_header.txt
if [[ $? != 0 ]]; then
    gst_err "non zero exit for bcftools"
    exit 1
fi
# <<<<<<<<<<<<<<<<<<<<<<<< Checking <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> Indexing >>>>>>>>>>>>>>>>>>>>>>>>
gst_log "Indexing ..."

if [[ -s "${pre}/merge_idx.tmp" && -s "${pre}/merge_idx-Query.tmp" && -s "${pre}/merge_idx-Ref.tmp" ]];then
    gst_warn "Found pre-exits index files, skip running."
else
    #? Create index for query and ref files
    # get index: 1chr 2posi 3id 4line_number 5q/r
    # bcftools query -f "%CHROM\t%POS\t%ID\n" $query | perl -pe 's/$/\t$.\tQuery/' > ${pre}/query_idx.tmp
    # bcftools query -f "%CHROM\t%POS\t%ID\n" $ref | perl -pe 's/$/\t$.\tRef/' > ${pre}/ref_idx.tmp
    csvtk cut -tTH -j $threads -f 1,2,3 $query | perl -lane '
        BEGIN{$type="Query";$,="\t";}
        $ord_in_chr=++$h{$F[0]};
        print "$_\t$ord_in_chr\t$type";
    ' > ${pre}/query_idx.tmp &&\
    csvtk cut -tTH -j $threads -f 1,2,3 $ref | perl -lane '
        BEGIN{$type="Ref";$,="\t";}
        $ord_in_chr=++$h{$F[0]};
        print "$_\t$ord_in_chr\t$type";
    ' > ${pre}/ref_idx.tmp &&\
    # merge query ref idx and sort and re-number
    # 1chr 2posi 3id 4raw_no 5q/r 6new_no
    cat ${pre}/query_idx.tmp ${pre}/ref_idx.tmp | csvtk sort -tTH -j $threads -k 1 -k 2:n |perl -lane '
        BEGIN{$,="\t";}
        $ord_in_chr=++$h{$F[0]};
        print "$_\t$ord_in_chr";
    ' > ${pre}/merge_idx.tmp &&\
    #? ${pre}/merge_idx-Query.tmp ${pre}/merge_idx-Ref.tmp
    csvtk split -tTH -j $threads -f 5 ${pre}/merge_idx.tmp -o ${pre}
    if [[ $? != 0 ]]; then
        gst_err "Non zero exit when generate the index files"
        exit 1
    fi
fi
# <<<<<<<<<<<<<<<<<<<<<<<< Indexing <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>
function Calc_LD () {
    query_id=$1
    #? Get ref vcf flanking query record and calc the LD
    query_rcd=$(grep -m 1 -w "$query_id" ${pre}/merge_idx-Query.tmp)
    query_ord=$(echo "$query_rcd" | cut -f 6)
    query_chr=$(echo "$query_rcd" | cut -f 1)
    query_region=$(echo "$query_rcd" | perl -F"\t" -lane 'print "$F[0]:$F[1]";')
    ref_just_before_query=$(
        cat ${pre}/merge_idx-Ref.tmp |\
        csvtk grep -tTH -p "$query_chr" -f 1 |\
        csvtk filter -tTH -f "6<$query_ord" |\
        tail -n 1
    )
    ref_just_after_query=$(
        cat ${pre}/merge_idx-Ref.tmp |\
        csvtk grep -tTH -p "$query_chr" -f 1 |\
        csvtk filter -tTH -f "6>$query_ord" |\
        head -n 1
    )
    ref_jb_ord=$( echo "$ref_just_before_query" | cut -f 4 )
    ref_ja_ord=$( echo "$ref_just_after_query" | cut -f 4 )
    #? set ord of the bondary to 0 or max
    [ -z "$ref_jb_ord" ] && ref_jb_ord=0;
    [ -z "$ref_ja_ord" ] && ref_ja_ord=9999999999;
    #? get ref records before/after query
    ref_after_records=$(
        cat ${pre}/merge_idx-Ref.tmp |\
        csvtk grep -tTH -p "$query_chr" -f 1 |\
        csvtk filter -tTH -f "4>=$ref_ja_ord" |\
        head -n $size
    )
    ref_after_num=$(echo "$ref_after_records" | wc -l)
    ref_before_num=$ref_jb_ord
    [ -z "$ref_after_records" ] && ref_after_num=0
    #? skip querys if there are not enough refs flanking
    if [[ "$ref_before_num" -lt "$size" || "$ref_after_num" -lt "$size" ]];then
        echo -e "${query_id}\tNaN\tNaN\tSkipped"
        gst_warn "SKIP ${query_id}: Not enough ref-records (Need:${size}/${size}, have:${ref_before_num}/${ref_after_num})"
    else
        ref_before_records=$(
            cat ${pre}/merge_idx-Ref.tmp |\
            csvtk grep -tTH -p "$query_chr" -f 1 |\
            csvtk filter -tTH -f "4<=$ref_jb_ord" |\
            tail -n $size
        )
        #? get ref records and query records and merge to one vcf
        ref_records=$( echo -e "${ref_before_records}\n${ref_after_records}")
        ref_id=$( echo "$ref_records" | cut -f 3 )
        ref_region=$( echo "$ref_records" | cut -f 1,2 | sed -n '1p;$p' | sed 'N;s/\n/\t/' | perl -F"\t" -lane 'print "$F[0]:$F[1]-$F[3]";')
        #? get ref vcf and format them: change format to "."
        ref_vcf=$(
            bcftools view -H -r $ref_region $ref 2>/dev/null |\
            csvtk grep -tTH -f 3 -P <(echo "$ref_id") |\
            perl -F"\t" -lane '
                BEGIN{$,="\t";}
                $F[7]=".";
                print @F;
            ')
        #? get query vcf and format it: change format col to ".", change ref/alt geno to T/A
        # gst_log "33"
        query_vcf=$(
            bcftools view -H -r $query_region $query 2>/dev/null |\
            csvtk grep -tTH -f 3 -p "$query_id" |\
            perl -F"\t" -lane '
                BEGIN{$,="\t";}
                $F[3]="T";
                $F[4]="A";
                $F[7]=".";
                print @F;
            ')
        #? merge ref query vcf, sort, header and generate tmp vcf file
        # gst_log "44"
        merge_vcf=$(echo -e "${ref_vcf}\n${query_vcf}" | sort -k1,1 -k2,2n | perl -lane '$F[1]=$.;$_=join("\t",@F);print;' )
        cat ${pre}/vcf_header.txt <(echo "${merge_vcf}") > ${pre}/tmp/${query_id}_flank${size}.vcf
        #? using vcftools: vcftools --vcf in.vcf --min-r2 0 --geno-r2 --out vcftools.r2
        # will generate out.geno.ld, redirect error to prevent out.log
        vcftools --vcf ${pre}/tmp/${query_id}_flank${size}.vcf --ld-window-bp 100000000 --min-r2 0 --geno-r2 --out ${pre}/LDs/${query_id}_flank${size}_hplv 2>/dev/null
        #? run Rscript to parse LD result
        # Rscript <Program> <mode> <haploview.LD> <size> <outdir>
        if [ ! -s "${pre}/Rcalc/${query_id}_ldtype.txt" ];then
            Rscript $rscript $mode ${pre}/LDs/${query_id}_flank${size}_hplv.geno.ld $size ${pre}/Rcalc
        else
            cat ${pre}/Rcalc/${query_id}_ldtype.txt
        fi
        rm -f ${pre}/tmp/${query_id}_flank${size}.vcf
    fi
    if [[ $? != 0 ]]; then
        gst_err "Non zero exit for $query_id"
        exit 1
    fi
    if [[ "$clean" == "TRUE" ]];then
        rm -f ${pre}/Rcalc/${query_id}_ldtype.txt ${pre}/LDs/${query_id}_flank${size}_hplv.geno.ld
    fi
}
export -f Calc_LD
# <<<<<<<<<<<<<<<<<<<<<<<< Main function <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> MAIN >>>>>>>>>>>>>>>>>>>>>>>>

parallel -j $threads --bar -k Calc_LD :::: <(cut -f 3 ${pre}/merge_idx-Query.tmp) 1> ${pre}_Total_records.tsv &&\
rm -rf ${pre}/*.tmp ${pre}/tmp ${pre}/vcf_header.txt &&\
gst_log "All Done."
if [[ $? != 0 ]]; then
    gst_err "Non zero exit"
    exit 1
fi

# <<<<<<<<<<<<<<<<<<<<<<<< MAIN <<<<<<<<<<<<<<<<<<<<<<<<
