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
export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix check_files_exists
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<
usage="
------------------------------------------------------------
PANZ freq enrichment analysis: input query and ref frequency file, and do enrichment analysis.

frequency file format:(no header)
    <Item_ID>  <Freq>
    ITEM1      123
    ITME2      888
    ...        ...

------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

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
"

if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export threads=2
export query=""
export ref=""
export database=""
export cutoff=0.05
export output="PANZ_KEGG_Enrich.tsv"
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -q|--query)
            #echo "set argument \"$1\" with value: $2" >&2
            query=$2
            shift 2
        ;;
        -r|--ref)
            #echo "set argument \"$1\" with value: $2" >&2
            ref=$2
            shift 2
        ;;
        # -d|--database)
        #     #echo "set argument \"$1\" with value: $2" >&2
        #     database=$2
        #     shift 2
        # ;;
        -c|--cutoff)
            #echo "set argument \"$1\" with value: $2" >&2
            cutoff=$2
            shift 2
        ;;
        -o|--output)
            #echo "set argument \"$1\" with value: $2" >&2
            output=$2
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
unset UNKOWN_ARGS
check_sftw_path csvtk parallel Rscript
check_var_numeric cutoff threads
check_files_exists $query $ref
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<
# USAGE:0 input1 input2
function fmt_freq () {
    # usage: fmt_freq output
    # use query ref database threads
    local output=$1
    local query_ko=""
    local all_ko=""
    local merge_ko_q=""
    local merge_ko_qr=""
    local ref_sum=""
    local query_sum=""
    all_ko=$(cat $ref $query | cut -f 1| sort -u)
    merge_ko_q=$(csvtk join -tTH -f 1 -k <(echo "$all_ko") $query | perl -F"\t" -lane '$F[1]=0 unless $F[1];$,="\t";print @F;' )
    merge_ko_qr=$(csvtk join -tTH -f 1 -k <(echo "$merge_ko_q") $ref | perl -F"\t" -lane '$F[2]=0 unless $F[2];$,="\t";print @F;')
    # ? add rest number to merge_ko_2
    query_sum=$(echo "$merge_ko_qr" | cut -f 2 | perl -lane '$sum+=$_;END{print $sum;}')
    ref_sum=$(echo "$merge_ko_qr" | cut -f 3 | perl -lane '$sum+=$_;END{print $sum;}')
    # ? ID pq pr lq lr
    echo "$merge_ko_qr" | csvtk mutate2 -L 0 -tTH -e "$query_sum - \$2" |\
    csvtk mutate2 -L 0 -tTH -e "$ref_sum - \$3" |\
    csvtk grep -tTH -f 1 -P <(cat $query | cut -f 1) -o $output
    if [ $? -ne 0 ];then gst_err "Get stat kegg failed: Non-zero exit"; exit 1;fi
}



function fisher_pvalue () {
    #usage: input line "KO p1 p2 rest1 rest2"
    local line=$1
    local id=""
    local present1=""
    local present2=""
    local rest1=""
    local rest2=""
    local cur_pv=""
    # ? gene_id No_total No_present No_loss
    echo "$line" | while read id present1 present2 rest1 rest2
        do
            check_var_empty id present1 present2 rest1 rest2
            # run fisher test in R
            cur_pv=$(Rscript -e 'argv=as.numeric(commandArgs(TRUE));
                p1=argv[1];p2=argv[2];l1=argv[3];l2=argv[4];
                x=c(p1,p2,l1,l2);dim(x)=c(2,2);
                pv=fisher.test(x,alternative="two.sided")$p.value;
                cat(pv);' $present1 $present2 $rest1 $rest2)
            check_var_empty $cur_pv
            echo -e "$id\t$present1\t$present2\t$rest1\t$rest2\t$cur_pv"
        done
}

function qvalue () {
    # ? usage: qvalue input output
    # ? input: file with header, and "Pvalue" in colum
    # ? ouput: file name
    local pv_file=$1
    local out_file=$2
    check_files_exists $pv_file
    check_var_empty out_file
    Rcmd=$(cat<<EOF
library(qvalue)
argv=commandArgs(TRUE)
tt=read.table(as.character(argv[1]),stringsAsFactors=F,header=T,check.names=F)
tt_qvalue=qvalue(as.numeric(tt\$Pvalue),pi0=1)
tt\$Qvalue=tt_qvalue\$qvalues
tt\$Lfdr=tt_qvalue\$lfdr
write.table(tt,file=as.character(argv[2]),quote=F,sep="\t",row.names=F)
EOF
)
    Rscript -e "$Rcmd" $pv_file $out_file
}

function filter_results () {
    # usage: $0 input output
    # ID P1 P2 L1 L2 Pv Qv Lfdr
    # <KO>  <Type>  <ratio_in_query>  <ratio_in_ref>  <raw_Pvalue>  <Qvalue>  <Lfdr>  <Description>
    local input=$1
    local output=$2
    check_files_exists $input
    sed '1d' $input | perl -F"\t" -lane '
        BEGIN{
            $,="\t";
            # $inputfile="$ENV{database}";
            # open(IN,"$inputfile") or die("Cannot open file: $inputfile");
            # while(<IN>){
            #     chomp;
            #     ($gene,$ko,$des)=split(/\t/,$_);
            #     @Ades=split(/\|/,$des);
            #     $fdes=join("\|",$Ades[-2,-1]);
            #     $hash{$ko}=$fdes;
            # }
        }
        $sig="Not-significant";
        $sum_q=$F[1]+$F[3];
        $sum_r=$F[2]+$F[4];
        $ratio_q=$F[1]/$sum_q;
        $ratio_r=$F[2]/$sum_r;
        $qr="$F[1]/$sum_q";
        $rr="$F[2]/$sum_r";
        $type="Enriched";
        $type="Diminished" if $ratio_q < $ratio_r;
        if($F[6] <= $ENV{cutoff}){
            # print $F[0],$type,$qr,$rr,@F[5..$#F],$hash{$F[0]} if $hash{$F[0]};
            # print $F[0],$type,$qr,$rr,@F[5..$#F],"Significant";
            $sig="Significant";
        }
        print $F[0],$type,$qr,$rr,$ratio_q/$ratio_r,@F[5..$#F],$sig' > $output
}
# function get_rep () {
#     # usage: $0 in_line
#     local line=$1
#     local cur_ko=$(echo "$line" | cut -f 1)
#     local cur_qgene=""
#     check_files_exists $map
#     check_var_empty cur_ko
#     cur_qgene=$(csvtk grep -tTH -f 2 -p "$cur_ko" ${output}.query_map_tmp | cut -f 1 | csvtk transpose)
#     check_var_empty cur_qgene
#     echo -e "${line}\t${cur_qgene}"
# }

export -f fisher_pvalue qvalue filter_results fmt_freq
# >>>>>>>>>>>>>>>>>>>>>>>> Main >>>>>>>>>>>>>>>>>>>>>>>>
# gst_log "Calc query and ref KO stats ..."

fmt_freq ${output}.fmt_freq_tmp

gst_log "Calc fisher's exact P value ..."
echo -e "ID\tP1\tP2\tL1\tL2\tPvalue" > ${output}.pvalue_tmp
parallel --bar -k fisher_pvalue :::: ${output}.fmt_freq_tmp >> ${output}.pvalue_tmp
if [ $? -ne 0 ];then gst_err "calc pvalue failed: Non-zero exit"; exit 1;fi

gst_log "Adjust Pvalues to Qvalues ..."
qvalue ${output}.pvalue_tmp ${output}.Qvalue_tmp

gst_log "Filter out Qvalue > $cutoff ..."
filter_results ${output}.Qvalue_tmp ${output}

if [ ! -s "${output}" ];then
    gst_warn "No data passed the filter cutoff: $cutoff"
    exit 1
fi

perl -pi -e '$_="Item\tType\tratio_in_query\tratio_in_ref\tfold_enrichment\traw_Pavalue\tQvalue\tLfdr\tSignificant_Tag\n".$_ if $.==1' ${output}
gst_log "All done ! Result in $output"
rm -f ${output}.*_tmp
# <<<<<<<<<<<<<<<<<<<<<<<< Main <<<<<<<<<<<<<<<<<<<<<<<<



