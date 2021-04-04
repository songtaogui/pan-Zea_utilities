#!/usr/bin/env bash

# set -o xtrace
# set -o errexit
set -o nounset
set -o pipefail

# >>>>>>>>>>>>>>>>>>>>>>>> Common functions >>>>>>>>>>>>>>>>>>>>>>>>
gst_log () {
    local info=$1
    echo -e "\033[36m[$(date +'%y-%m-%d %H:%M')]\033[0m $info" >&2
}

gst_rcd () {
    local info=$1
    echo -e "\033[32m>>>------------>\033[0m $info" >&2
}

gst_err () {
    local info=$1
    echo -e "\033[31m\033[7m[ERROR]\033[0m --> $info" >&2
}

gst_warn () {
    local info=$1
    echo -e "\033[35m[WARNING]\033[0m --> $info" >&2
}

check_files_exists(){
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

check_abs_path() {
    local var_cc=1
    local check_file=""
    for check_file in "$@";do
        if [[ "${check_file:0:1}" != "/" ]]; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> $check_file was not an ABSOLUTE path." >&2
            let var_cc++
        fi
    done
    [ "$var_cc" -ne 1 ] && exit 1
}

check_sftw_path(){
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

check_var_empty () {
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

check_var_numeric () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        # add ${var#prefix} substitution to trim sign
        case ${var#[-+]} in
            '')
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name is empty: '$var' " >&2
                let var_cc++ ;;
            *.*.*)
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name has more than one decimal point: '$var' " >&2
                let var_cc++ ;;
            *[!0-9.]*)
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name has a non-digit somewhere in it: '$var' " >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}

check_suffix () {
    check_suffix_file=$( basename $1 )
    check_suffix=$2
    # add x incase file has no suffix
    if [[ "${check_suffix_file##*.}"x != "$check_suffix"x ]];then
        echo "[ERROR] --> $check_suffix_file should have suffix: '$check_suffix'." >&2
        exit 1
    fi
}

export -f gst_log gst_rcd gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix check_files_exists check_abs_path
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<
usage=$(
cat <<EOF
------------------------------------------------------------
Perform regional association analysis of genic regions for PANZ.
A wrapper of MAGMA
------------------------------------------------------------
Dependence: MAGMA plink
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

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

EOF
)
if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export annote=""
export ref=""
export gene_sets=""
export pheno=""
export covar=""
export region=""
export out="PANZ_Gene_GWAS_out"
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -r|--ref)
            #echo "set argument \"$1\" with value: $2" >&2
            ref=$2
            shift 2
        ;;
        --annote)
            #echo "set argument \"$1\" with value: $2" >&2
            annote=$2
            shift 2
        ;;
        --gene_sets)
            #echo "set argument \"$1\" with value: $2" >&2
            gene_sets=$2
            shift 2
        ;;
        --pheno)
            pheno=$2
            shift 2
        ;;
        --covar)
            #echo "set argument \"$1\" with value: $2" >&2
            covar=$2
            shift 2
        ;;
        --region)
            region=$2
            shift 2
        ;;
        -o|--out)
            out=$2
            shift 2
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
unset UNKOWN_ARGS # restore UNKOWN_ARGS params
# ! Check if required vars are legal
check_var_empty ref annote gene_sets pheno covar out
check_files_exists $ref.bed $ref.bim $ref.fam $annote $gene_sets $pheno $covar
if [ -n "$region" ]; then
    echo "$region" | grep -P "^\S+\:\d+\-\d+$" > /dev/null
    if [ $? -ne 0 ];then gst_err "Wrong region format: $region."; exit 1;fi
fi

# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>> Main >>>>>>>>>>>>>>>>>>>>>>>>
gst_log "All start. Will output to $out ..."
# ? make output dir
mkdir -p $out


gst_rcd "Get ref ..."
if [ ! -s "$out/get_ref_done" ];then
    # ? get ref bed
    if [ -n "$region" ];then
        gst_rcd "Subset with $region ..."
        echo "$region" | sed 's/[\:\-]/\t/g' | while read chr start end
        do
            export chr start end
            gst_rcd "Get ref bed ..."
            plink --make-bed --bfile $ref --chr $chr --from-bp $start --to-bp $end --out $out/ref 1>$out/plink.log 2>&1
            if [ $? -ne 0 ];then gst_err "subset plink bed failed: check $out/plink.log for details";rm -f $out/ref.{bed,bim,fam}; exit 1;fi
            gst_rcd "get ref annote .."
            cat $annote | perl -F"\t" -lane '
                BEGIN{use List::Util qw/max min/;$,="\t";}
                if(/^#/){print;next;}
                ($c,$s,$e)=split(":",$F[1]);
                #max(A.start,B.start)<=min(A.end,B.end)
                if( $c eq $ENV{chr} && max($s,$ENV{start})<=min($e,$ENV{end}) ){
                    print @F;
                }
            ' > $out/ref.v_annote
            if [ $? -ne 0 ];then gst_err "get subset annote failed: Non-zero exit";rm -f $out/ref.v_annote; exit 1;fi
        done
    else
        ln -f -s $ref.bed $out/ref.bed &&\
        ln -f -s $ref.bim $out/ref.bim &&\
        ln -f -s $ref.fam $out/ref.fam &&\
        ln -f -s $annote $out/ref.v_annote
        if [ $? -ne 0 ];then gst_err "link ref failed: Non-zero exit";rm -f $out/ref.*; exit 1;fi
    fi
    if [ $? -ne 0 ];then gst_err "Get ref failed: Non-zero exit"; exit 1;fi
    echo "done" > $out/get_ref_done
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
fi

gst_rcd "Gene analysis ..."
magma="magma"
if [ ! -s "$out/gene_ana_done" ];then
    # magma --bfile $out/ref --covar file=$covar --pheno file=$pheno --gene-annot $out/ref.v_annote --gene-settings snp-max-miss=0.25 adap-permp --seed 1234 --out $out/${magma} 1> $out/gene_ana.log 2>&1
    magma --bfile $out/ref --covar file=$covar --pheno file=$pheno --gene-annot $out/ref.v_annote --gene-settings snp-max-miss=0.25 --out $out/${magma} 1> $out/gene_ana.log 2>&1
    if [ $? -ne 0 ];then gst_err "gene_analysis failed: check $out/gene_ana.log for details"; exit 1;fi
    echo "done" > $out/gene_ana_done
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
fi

gst_rcd "Gene set analysis ..."
if [ ! -s "$out/gene_set_done" ];then
    magma --gene-results $out/${magma}.genes.raw --set-annot $gene_sets --out $out/$magma 1> $out/gene_set.log 2>&1
    if [ $? -ne 0 ];then gst_err "gene set analysis failed: check $out/gene_set.log for details"; exit 1;fi
    echo "done" > $out/gene_set_done
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
fi


gst_log "All Done!"

# <<<<<<<<<<<<<<<<<<<<<<<< Main <<<<<<<<<<<<<<<<<<<<<<<<

