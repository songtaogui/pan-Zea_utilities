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
            gst_err "No file: $related_file"
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
            gst_log "Program not in PATH: $tp_program"
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
                gst_err "$var_name is empty: '$var'"
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
                gst_err "$var_name is empty: '$var'"
                let var_cc++ ;;
            *.*.*)
                gst_err "$var_name has more than one decimal point: '$var'"
                let var_cc++ ;;
            *[!0-9.]*)
                gst_err "$var_name has a non-digit somewhere in it: '$var'"
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
        gst_err "$check_suffix_file should have suffix: '$check_suffix'."
        exit 1
    fi
}
export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<

usage="
------------------------------------------------------------
Determine core and dispensable genes/gene families
using pvalue of binomial test of gene loss rate.
------------------------------------------------------------
Dependency: parallel perl Rscript
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
    -r, --lossrate   <0-1>    [O]    min loss rate to determine a gene dispensable (Default: 0.01)
    -p,--pvalue      <0-1>    [O]    max pvalue cutoff to determine a gene dispensable (Default: 0.05)
    -o,--output      <str>    [O]    output file (default: core_disp_out.tsv)
    --prefix         <str>    [O]    prefix tag of output gene types (default: "PANZ")
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
output=core_disp_out.tsv
export lossrate=0.01
export pvalue=0.05
export prefix="PANZ"
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
        -r|--lossrate)
            #echo "set argument \"$1\" with value: $2" >&2
            lossrate=$2
            shift 2
        ;;
        -p|--pvalue)
            #echo "set argument \"$1\" with value: $2" >&2
            pvalue=$2
            shift 2
        ;;
        -o|--output)
            #echo "set argument \"$1\" with value: $2" >&2
            output=$2
            shift 2
        ;;
        --prefix)
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
check_var_empty matrix output prefix
check_var_numeric lossrate pvalue threads
check_files_exists $matrix
check_sftw_path parallel perl Rscript
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> core function >>>>>>>>>>>>>>>>>>>>>>>>
function core_disp () {
    local line=$1
    local stat
    local id
    local total
    local present
    local lost
    local cur_pv
    local cur_rcd
    # ? gene_id No_total No_present No_loss
    cur_rcd=$(echo "$line" | perl -lane '
        $n_p=grep { $_ == 1 } @F[1..$#F];
        $n_l=grep { $_ == 0 } @F[1..$#F];
        $n_t=$#F;
        $,="\t";
        print $F[0],$n_t,$n_p,$n_l;
    ' | while read id total present lost
        do
            check_var_empty id total present lost
            # run binorm test in R
            cur_pv=$(Rscript -e 'argv=as.numeric(commandArgs(TRUE));
                pv=binom.test(argv[1],argv[2],argv[3],alternative="greater")$p.value;
                cat(pv);' $lost $total $lossrate )
            echo -e "$id\t$total\t$present\t$lost\t$cur_pv" | perl -F"\t" -lane '
                ($id,$total,$present,$lost,$cpv)=@F;
                $type="NA";
                if($present == $total){
                    $type="Core";
                }elsif($cpv <= $ENV{pvalue}){
                    $type="Dispensable";
                }else{
                    $type="Candidate-Core";
                }
                print "$_\t$ENV{prefix}_$type";'
        done
    )
    echo "$cur_rcd"
}
export -f core_disp gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix
# <<<<<<<<<<<<<<<<<<<<<<<< core function <<<<<<<<<<<<<<<<<<<<<<<<
gst_log "Program start ...
--------------------------
input:      $matrix
prefix:     $prefix
output:     $output
loss rate:  $lossrate
p-value:    $pvalue
--------------------------
"
if [ -s "$output" ];then
    gst_warn "Skip running: $output exists."
else
    echo -e "ID\tTotal_num\tPresent_num\tLost_num\tPvalue\tType" > $output
    parallel --bar -j $threads -k core_disp :::: <(sed '1d' $matrix) >> $output
    if [ $? -ne 0 ];then gst_err "parallel run core_disp failed: Non-zero exit";rm -f $output; exit 1;fi
fi
gst_log "All Done !"
