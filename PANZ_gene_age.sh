#!/usr/bin/env bash
set -uo pipefail

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
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> No file: $related_file $" >&2
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
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> Program not in PATH: $tp_program $" >&2
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
Infering Gene ages using methods discribed in:
    Zhang, Y. E. et al. Chromosomal redistribution of male-biased
    genes in mammalian evolution with two bursts of gene gain on
    the X chromosome. PLoS Biol. 8, e1000494 (2010).
------------------------------------------------------------
Dependency in PATH:
    diamond, csvtk, taxonkit, perl, GNU-parallel, blastdbcmd
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                       show help and exit.
    -d, --database   <str>    [O]    Path to NCBI-NR protein database (Default: $WORK/ref/NR/nr/nr)
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
"

if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
database=$WORK/ref/NR/nr/nr
input=
output=PANZ_Gene_ages
taxnomy=
evalue=1e-5
identity=0.3
threads=2
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -d|--database)
            #echo "set argument \"$1\" with value: $2" >&2
            database=$2
            shift 2
        ;;
        -i|--input)
            #echo "set argument \"$1\" with value: $2" >&2
            input=$2
            shift 2
        ;;
        -o|--output)
            #echo "set argument \"$1\" with value: $2" >&2
            output=$2
            shift 2
        ;;
        --taxnomy)
            #echo "set argument \"$1\" with value: $2" >&2
            taxnomy=$2
            shift 2
        ;;
        --evalue)
            #echo "set argument \"$1\" with value: $2" >&2
            evalue=$2
            shift 2
        ;;
        --identity)
            #echo "set argument \"$1\" with value: $2" >&2
            identity=$2
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
# ! Check if required vars are legal
check_var_empty database input taxnomy evalue output
check_var_numeric threads identity
check_files_exists ${database}.pal $input $taxnomy
check_sftw_path diamond csvtk taxonkit perl parallel blastdbcmd

# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> work functions >>>>>>>>>>>>>>>>>>>>>>>>
function get_tid_list () {
    # usage: get_tid_list $label $tid $order $threads
    local label=$1
    local tid=$2
    local order=$3
    local threads=$4
    check_var_empty label
    check_var_numeric tid order threads
    local fmt_order=$(printf "${output}_%07d" $order)
    # get tid list using taxonkit and append label to the end, then output to fmt_order.tid_list
    gst_log "Get taxnomy ID for $label ..."
    if [ -s "${fmt_order}.tid_list" ];then
        gst_warn "Output file exists, skip running"
    else
        taxonkit list -j $threads --ids $tid --indent "" | parallel --pipe -q -k -j $threads sed "/^[[:space:]]*$/d;s/\$/\t$label/" >${fmt_order}.tid_list
        if [ $? -ne 0 ];then gst_err "get tid list failed: Non-zero exit"; rm -f ${fmt_order}.tid_list; exit 1;fi
    fi
    gst_log "Resulted file: ${fmt_order}.tid_list"
}
function get_tid_dict () {
    # usage: get_tid_dict $output
    local output=$1
    # use tac to put newest order first
    gst_log "Get unique taxnomy file ..."
    if [ -s "$output" ];then
        gst_warn "Output file exists, skip running"
    else
        cat $(ls *.tid_list | tac) | perl -F"\t" -lane 'print if ++$h{$F[0]}==1' > $output
        if [ $? -ne 0 ];then gst_err "get_tid_dict failed: Non-zero exit"; rm -f $output; exit 1;fi
    fi
    gst_log "Resulted unique taxnomy file: $output"
}

function get_ref_fa () {
    # usage: get_ref_fa $database $tid_dict $output
    local db=$1
    local tid_dict=$2
    local output=$3
    export cur_all_taxnomy_dictionary=$tid_dict
    gst_log "Getting reference fasta ..."
    if [ -s "$output" ];then
        gst_warn "Output file exists, skip running"
    else
        blastdbcmd -db $db -entry all -outfmt "%g,%T,%s" | perl -F"," -lane '
        BEGIN{
            open(IN,"$ENV{cur_all_taxnomy_dictionary}") or die "Cannot open file: $ENV{cur_all_taxnomy_dictionary}";
            while(<IN>){
                chomp;
                ($tid,$label)=split(/\t/,$_);
                $tid_hash{$tid}=$label;
            }
        }
        $uniq_tag="$tid_hash{$F[1]},$F[2]";
        # ? keep only one record for records with same group type and same sequences
        # ? output fa: group_type,gid,tid
        print ">$tid_hash{$F[1]},$F[0],$F[1]\n$F[2]" if $tid_hash{$F[1]} and ++$h{$uniq_tag}==1;
        ' > $output
        if [ $? -ne 0 ];then gst_err "get_ref_fa failed: Non-zero exit"; rm -f $output; exit 1;fi
    fi
    gst_log "Resulted reference fasta file: $output
    Name Format: Group_type,gi,taxid"
    unset cur_all_taxnomy_dictionary
}

function run_diamond_blastp () {
    # usage: run_diamond_blastp ref.fa query.fa evalue identity threads outputs
    local ref=$1
    local query=$2
    local ev=$3
    local pid=$4
    local cpu=$5
    local output=$6
    gst_log "Make diamond index ..."
    if [ ! -s "${ref}.dmnd" ];then
        diamond makedb --in $ref -d $ref --threads $cpu
        if [ $? -ne 0 ];then gst_err "make diamond db failed: Non-zero exit"; rm -f ${ref}.dmnd; exit 1;fi
    fi
    gst_log "Run diamond blastp alignment ... "
    if [ -s "$output" ];then
        gst_warn "Output file exists, skip running"
    else
        diamond blastp --threads $cpu --query $query --db $ref -o $output --evalue $ev --id $pid --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp
        if [ $? -ne 0 ];then gst_err "diamond blastp failed: Non-zero exit";rm -f $output; exit 1;fi
    fi
    gst_log "Resulted blastp tabular file: $output
    format:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp
    "
}
function fill_sort_age () {
    # usage: echo -e "gene\tage1;age2;age3" | fill_sort_age
    while read id ages
    do
        cur_age_line=$(echo "$ages" | sed 's/;/\n/g' | perl -lane 'print "$_\t1";')
        cur_age_code=$( csvtk join -tTH -k <(echo "$all_age_line") <(echo "$cur_age_line") |\
        perl -F"\t" -lane '$a="0";$a="1" if $F[1];print "$F[0]:$a";'|perl -pne 's/\n/\t/ unless eof')
        echo -e "$id\t$cur_age_code"
    done
}

function get_final_age_all () {
    # usage: get_final_age_all $diamond_out $threads $output
    # output: Gene1 age4:0 age3:1 age2:1 age1:1
    local dmdout=$1
    local threads=$2
    local output=$3
    gst_log "Get final age details ..."
    if [ -s "$output" ];then
        gst_warn "Output file exists, skip running"
    else
        cut -f 1,2 $dmdout |\
        perl -F"\t" -lane '$F[1]=~s/,.*//; print "$F[0]\t$F[1]" if ++$h{"$F[0],$F[1]"}==1;' |\
        csvtk collapse -tTH -f 1 -v 2 -j $threads | parallel --pipe -k -j $threads fill_sort_age >$output
        if [ $? -ne 0 ];then gst_err "get final age failed: Non-zero exit";rm -f $output; exit 1;fi
    fi
    gst_log "age detail file: $output"
}
function get_final_age_rep () {
    # usage: get_final_age_rep $age_detail_file $threads $output
    local input=$1
    local threads=$2
    local output=$3
    gst_log "Get final age represents ..."
    gst_warn "This step will remove '[-_]T\d+$' in the protein IDs, if this influence your real gene ID, please comment out Line 14 of 'get_final_age_rep' function"
    if [ -s "$output" ];then
        gst_warn "Output file exists, skip running"
    else
        cat $input |\
        parallel --pipe -j $threads -k -q perl -lane '
            $F[0]=~s/[\-_]T\d+$//;
            foreach $a (@F[1..$#F]){
                $age=$a;
                $age=~s/:.*//;
                if($a=~/1$/){
                    print "$F[0]\t$age";break;
                }
            }
        ' | csvtk sort -tTH -j $threads -k 1:r -k 2 | perl -lane 'print if ++$h{$F[0]}==1' >$output
        if [ $? -ne 0 ];then gst_err "get final age represents failed: Non-zero exit"; rm -f $output; exit 1;fi
    fi
    gst_log "Final age represent file: $output"
}
export -f get_tid_list get_tid_dict get_ref_fa run_diamond_blastp fill_sort_age get_final_age_all get_final_age_rep
# <<<<<<<<<<<<<<<<<<<<<<<< work functions <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> Main >>>>>>>>>>>>>>>>>>>>>>>>
cat $taxnomy | while read label tid order
do
    get_tid_list $label $tid $order $threads
done
get_tid_dict ${output}.taxnomyID_dict
get_ref_fa $database ${output}.taxnomyID_dict ${output}_NR_ref_seq.fa
run_diamond_blastp ${output}_NR_ref_seq.fa $input $evalue $identity $threads ${output}_diamond_blastp_out
export all_age_line=$(sort -k 3,3n $taxnomy | cut -f 1)
get_final_age_all ${output}_diamond_blastp_out $threads ${output}_all_age.tsv
get_final_age_rep ${output}_all_age.tsv $threads ${output}_age_represent.tsv
# <<<<<<<<<<<<<<<<<<<<<<<< Main <<<<<<<<<<<<<<<<<<<<<<<<

