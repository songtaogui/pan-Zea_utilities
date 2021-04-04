#!/usr/bin/env bash

# >>>>>>>>>>>>>>>>>>>>>>>> Common functions >>>>>>>>>>>>>>>>>>>>>>>>
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
function check_files_exists(){
    local num_related_file=1
    local related_file=""
    for related_file in  "$@"
    do
        if [[ ! -s "$related_file" ]]; then
            echo "\033[31m\033[7m[ERROR]\033[0m --> No file: $related_file $" >&2
            let num_related_file++
        fi
    done
    [ "$num_related_file" -ne 1 ] && exit 1
}
function check_sftw_path(){
    local num_tp_program=1
    local tp_program=""
    for tp_program in "$@"
    do
        if ! which $tp_program >/dev/null 2>&1 ; then
            echo "\033[31m\033[7m[ERROR]\033[0m --> Program not in PATH: $tp_program $" >&2
            let num_tp_program++
        fi
    done
    [ "$num_tp_program" -ne 1 ] && exit 1
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
function check_var_numeric () {
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
function check_suffix () {
    check_suffix_file=$( basename $1 )
    check_suffix=$2
    # add x incase file has no suffix
    if [[ "${check_suffix_file##*.}"x != "$check_suffix"x ]];then
        echo "[ERROR] --> $check_suffix_file should have suffix: '$check_suffix'." >&2
        exit 1
    fi
}
export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<
usage="
------------------------------------------------------------
Determine subgroup unbalanced genes/gene families
using adjusted-pvalue (q-value, Storey's Method) of two-sided
Fisher's exact test of gene PAV matrix among different subgroups.
------------------------------------------------------------
Dependencies: GNU-parallel perl Rscript csvtk (All should be in PATH)
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

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
    -p,--prefix      <str>    [O]    prefix tag of outputs, you may use subgroup name (default: "PANZ")
    -f,--FDR         <0-1>    [O]    FDR cutoff to determine a gene unbalance (Default: 0.05)
    --localfdr                [O]    Use local FDR values (default: Global FDR, aka qvalue)
    -o,--output      <str>    [O]    output file (default: PAV_subgroup_unbalanced_out.tsv)
    -t, --threads    <num>    [O]    set threads (default: 2)
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
matrix=
output=PAV_subgroup_unbalanced_out.tsv
subgroup=
export prefix="PANZ"
export fdr=0.05
localfdr=FALSE
threads=2
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -m|--matrix)
            #echo "set argument \"$1\" with value: $2" >&2
            matrix=$2
            shift 2
        ;;
        -s|--subgroup)
            #echo "set argument \"$1\" with value: $2" >&2
            subgroup=$2
            shift 2
        ;;
        -f|--FDR)
            #echo "set argument \"$1\" with value: $2" >&2
            fdr=$2
            shift 2
        ;;
        -o|--output)
            #echo "set argument \"$1\" with value: $2" >&2
            output=$2
            shift 2
        ;;
        --localfdr)
            #echo "set argument \"$1\" with value: $2" >&2
            localfdr=TURE
            shift 1
        ;;
        -p|--prefix)
            #echo "set argument \"$1\" with value: $2" >&2
            prefix=$2
            shift 2
        ;;
        -t|--threads)
            #echo "set argument \"$1\" with value: $2" >&2
            threads=$2
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
# Check if required vars are legal
check_var_empty matrix output prefix subgroup localfdr
check_var_numeric fdr threads
check_files_exists $matrix $subgroup
check_sftw_path parallel perl Rscript csvtk
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> core function >>>>>>>>>>>>>>>>>>>>>>>>
# ? input: ID present1 presen2 lost1 lost2
function fisher_pvalue () {
    local line=$1
    local id=""
    local present1=""
    local present2=""
    local lost1=""
    local lost2=""
    local cur_pv=""
    # ? gene_id No_total No_present No_loss
    echo "$line" | while read id present1 present2 lost1 lost2
        do
            check_var_empty id present1 present2 lost1 lost2
            # run fisher test in R
            cur_pv=$(Rscript -e 'argv=as.numeric(commandArgs(TRUE));
                p1=argv[1];p2=argv[2];l1=argv[3];l2=argv[4];
                x=c(p1,p2,l1,l2);dim(x)=c(2,2);
                pv=fisher.test(x,alternative="two.sided")$p.value;
                cat(pv);' $present1 $present2 $lost1 $lost2)
            check_var_empty $cur_pv
            echo -e "$id\t$present1\t$present2\t$lost1\t$lost2\t$cur_pv"
        done
}
function stat_pav () {
    perl -lane '$n_p=grep { $_ == 1 } @F[1..$#F];
        $n_l=grep { $_ == 0 } @F[1..$#F];
        $,="\t";
        print $F[0],$n_p,$n_l;
    '
}
function qvalue () {
    # ? usage: qvalue input output
    # ? input: file with header, and "pvalue" in colum
    # ? ouput: file name
    local pv_file=$1
    local out_file=$2
    check_files_exists $pv_file
    check_var_empty out_file
    Rcmd=$(cat<<EOF
library(qvalue)
argv=commandArgs(TRUE)
tt=read.table(as.character(argv[1]),stringsAsFactors=F,header=T,check.names=F)
tt_qvalue=qvalue(tt\$Pvalue)
tt\$Qvalue=tt_qvalue\$qvalues
tt\$Lfdr=tt_qvalue\$lfdr
write.table(tt,file=as.character(argv[2]),quote=F,sep="\t",row.names=F)
EOF
)
    Rscript -e "$Rcmd" $pv_file $out_file
}

function determine_type () {
    perl -lane '$cur_fdr=$F[$ENV{fdr_idx}];$type="Unknown";$type="$ENV{prefix}_Unbalanced" if $cur_fdr <= $ENV{fdr};print "$_\t$type";'
}
export -f fisher_pvalue stat_pav qvalue determine_type
# <<<<<<<<<<<<<<<<<<<<<<<< core function <<<<<<<<<<<<<<<<<<<<<<<<
gst_log "Program start ...
--------------------------
input:      $matrix
subgroup:   $subgroup
prefix:     $prefix
output:     $output
fdr:        $fdr
local:      $localfdr
--------------------------
"
if [ -s "$output" ];then
    gst_warn "Skip running: $output exists."
else
    # ? calc pop1 and pop2 present and lost numbers
    subgroup_csv=$(cat $subgroup | perl -pane 's/\n/,/ unless eof')
    gst_log "Get main group PAV stats ..."
    sed '1d' $matrix | parallel --pipe --jobs $threads -k stat_pav > ${output}_maingroup_stat.tmp
    if [ $? -ne 0 ];then gst_err "get main group stat failed: Non-zero exit"; rm -f ${output}_maingroup_stat.tmp; exit 1;fi
    gst_log "Get subgroup PAV stats ..."
    first_col_id=$(head -1 $matrix | cut -f 1)
    check_var_empty first_col_id
    csvtk cut -tT -f "$first_col_id,$subgroup_csv" -j $threads $matrix | sed '1d' | parallel --pipe --jobs $threads -k stat_pav > ${output}_subgroup_stat.tmp
    if [ $? -ne 0 ];then gst_err "get subgroup stat failed: Non-zero exit"; rm -f ${output}_subgroup_stat.tmp; exit 1;fi
    # ? merge main and sub group stat
    csvtk join -tTH -j $threads -f 1 -k ${output}_maingroup_stat.tmp ${output}_subgroup_stat.tmp | csvtk cut -tTH -f 1,2,4,3,5 -j $threads -o ${output}_ppll.tmp
    if [ $? -ne 0 ];then gst_err "get merged stats failed: Non-zero exit"; rm -f ${output}_ppll.tmp exit 1;fi
    # ? calc fisher pv
    gst_log "Calc fisher test pvalues ..."
    echo -e "ID\tPresent1\tPresent2\tLost1\tLost2\tPvalue" > ${output}_fisher_pv.tmp
    parallel --bar -j $threads -k fisher_pvalue :::: ${output}_ppll.tmp >> ${output}_fisher_pv.tmp
    if [ $? -ne 0 ];then gst_err "parallel run fisher_pvalue failed: Non-zero exit";rm -f ${output}_fisher_pv.tmp; exit 1;fi
    # ? correct with qvalue
    qvalue ${output}_fisher_pv.tmp ${ouput}_qvalue.tmp
    if [ $? -ne 0 ];then gst_err "corrrect qvalue failed: Non-zero exit"; rm -f ${ouput}_qvalue.tmp; exit 1;fi
    # ? decide if unbalance based on qvalue
    gst_log "Determine balance type based on fdr ..."
    export fdr_idx=6
    [ "$localfdr" == "TRUE" ] && fdr_idx=7
    echo -e "ID\tPresent1\tPresent2\tLost1\tLost2\tPvalue\tQvalue\tLfdr\t${prefix}" > ${output}
    sed '1d' ${ouput}_qvalue.tmp | parallel --pipe -j $threads -k determine_type >> ${output}
    if [ $? -ne 0 ];then gst_err "Determine type failed: Non-zero exit"; rm -f ${output}; exit 1;fi
fi
csvtk cut -tTH -j $threads -f 1,9 -o ${output}.summary ${output}
if [ $? -ne 0 ];then gst_err "get summary failed: Non-zero exit"; rm -f ${output}.summary; exit 1;fi
gst_log "All Done !
--------------------------------------
Detailed output:    ${output}
Summary output:     ${output}.summary
--------------------------------------
"



