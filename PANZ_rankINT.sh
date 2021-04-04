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
PANZ Rank-Based Inverse Normal Transformation
------------------------------------------------------------
Dependence: Rscript
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

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

EOF
)
if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export trait=""
export col=2
export na_str="NA"
export out=""
export SW_test="FALSE"
export adj=0.5
export outdir=$PWD
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -t|--trait)
            #echo "set argument \"$1\" with value: $2" >&2
            trait=$2
            shift 2
        ;;
        -c|--col)
            #echo "set argument \"$1\" with value: $2" >&2
            col=$2
            if [[ "$col" -lt 1 ]];then
                gst_warn "Trait column <= 1, which is unusual, please check."
                exit 1
            fi
            shift 2
        ;;
        --na)
            #echo "set argument \"$1\" with value: $2" >&2
            na_str=$2
            shift 2
        ;;
        -o|--out)
            #echo "set argument \"$1\" with value: $2" >&2
            out=$2
            outdir=$(dirname $out)
            if [ -z "$outdir" ]; then
                outdir=$PWD
            fi
            shift 2
        ;;
        --adj)
            #echo "set argument \"$1\" with value: $2" >&2
            adj=$2
            shift 2
        ;;
        --SW_test)
            #echo "set argument \"$1\" with value: $2" >&2
            SW_test="TRUE"
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
unset UNKOWN_ARGS # restore UNKOWN_ARGS params
# ! Check if required vars are legal
check_var_empty trait col na_str outdir adj
check_var_numeric col adj
check_files_exists $trait
if [ ! -w "$outdir" ]; then
    gst_err "The output dir $outdir is not writable".
    exit 1
fi
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# USAGE:$0 non_NA_trait.tsv
rank_int_r () {
    local trait_name=$(head -n 1 $trait | cut -f "$col" )
    local outf=$trait_name
    if [ -n "$out" ]; then
        outf=$out;
    fi
    check_var_empty outf
    cut -f "1,$col" $trait > $outdir/$outf.tmp
    if [ $? -ne 0 ];then gst_err "cut trait matrix failed: Non-zero exit"; rm -f $outdir/$out.tmp; exit 1;fi
    local in=$outdir/$outf.tmp
    local oo=$outdir/$outf
    check_files_exists $in
    gst_log "Dealing with $trait_name ..."
    head -n 1 $in > $oo.RINT.tsv
    if [ $? -ne 0 ];then gst_err "Write to $oo.RINT.tsv failed: Non-zero exit"; rm -f $oo.RINT.tsv; exit 1;fi
    Rscript -e '
        argv=as.character(commandArgs(TRUE));
        rankTransPheno <- function(pheno,para_c)
            {
                pheno<-qnorm((rank(pheno)-para_c)/(length(pheno)-2*para_c+1))
                return(pheno)
            }
        # ? input
        tt=read.table(argv[1],sep="\t",header=T,na.string=argv[2],stringsAsFactors=F)
        # tt.value=tt[[as.numeric(argv[3])]]
        tt.value=tt[[2]]
        # ? SW test
        if(argv[6] == "TRUE"){
            sw=shapiro.test(tt.value[!is.na(tt.value)])
            cat(
                paste0("Method = ",sw$method,"\n",names(sw$statistic)," = ",sw$statistic,"\nP_value = ",sw$p.value,"\n"),
                file=paste0(argv[5],".sw_test")
            )
        }
        # ? Rank trans
        # tt[[as.numeric(argv[3])]][!is.na(tt.value)]=rankTransPheno(tt.value[!is.na(tt.value)],as.numeric(argv[4]))
        tt[[2]][!is.na(tt.value)]=rankTransPheno(tt.value[!is.na(tt.value)],as.numeric(argv[4]))
        write.table(tt,file=paste0(argv[5],".RINT.tsv"),sep="\t",quote=F,row.names=F,col.names=F,append=T)
    ' $in $na_str $col $adj $oo $SW_test
    if [ $? -ne 0 ];then gst_err "rankTrans in R failed: Non-zero exit"; exit 1;fi
    rm -f $in
    if [ $? -ne 0 ];then gst_err "rm $in failed: Non-zero exit"; exit 1;fi
    gst_log "Done! \nOutput: $oo.RINT.tsv"
}
export -f rank_int_r

# >>>>>>>>>>>>>>>>>>>>>>>> Main >>>>>>>>>>>>>>>>>>>>>>>>
rank_int_r
# <<<<<<<<<<<<<<<<<<<<<<<< Main <<<<<<<<<<<<<<<<<<<<<<<<
