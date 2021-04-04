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
Partition genetic variants and calculate variant heritability
using GWAS Summary Statistic.
This is a wrapper of LDAK pipelines.
------------------------------------------------------------
Dependence: LDAK plink csvtk parallel bedtools bgzip perl5
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -t, --threads       <num>   [O]     set threads (default: 2)
# I/O
    -v, --gvcf          <str>   [R]     bgzipped and index vcf file. SHOULD BE ABSOLUTE PATH.
    -p, --pheno         <str>   [R]     plink format pheno file. SHOULD BE ABSOLUTE PATH.
    --pheno_tag         <str>   [R]     A string tagging the pheno file (e.g. "Cob_color"), used to generate different reml outputs, in order to multiple-run of several phenotypes based on same gvcf files.
    -c, --covar         <str>   [R]     covar file in plink format. SHOULD BE ABSOLUTE PATH.
    -o, --out           <str>   [O]     output prefix (default: PANZ)
# mode
    --part_tsv          <str>   [O]     ABSOLUTE path of the files included Variant IDs and patition rules, will create a dir based on the prefix of this file.
                                    format (With header, header of Col1 must be VID, header of other cols would be used as feature tag in output. Mark excluded data with "NA"):
                                        VID        Annotations     MAF_Bins
                                        SNP1          CDS_SNPs     MAF0-0d05
                                        SNP2   INTERGENIC_SNPS     MAF0d05-0d1
                                        SNP3          CDS_SNPS     NA
    --part_bed          <str>   [O]     ABSOLUTE path of the bed files if you would like to estimatie the patition h2 by regions, will create a dir based on the prefix of this file.
                                    format (UCSC bed format, no header):
                                        <chr>   <start0>    <end>   <Feature_ID>
                                        1       12345       13456   TE_Rich_Region1
                                        2       22345       33456   Gene_Rich_Region1
                                        3       32345       43456   TE_Rich_Region1
# prune & keep
    --keep                      [O]     Keep only unrelated individuals for h2 estimation. (default: off, that is use all individuals provided)
    --window_kb         <int>   [O]     Window size to compute allelic correlations to thin variants.(default:1000)
    --window_prune      <0-1>   [O]     Max correlation squared cut-off to thin variants. (default: 0.2)
# PCA
    --PCA                       [O]     Use PCs as additional covar, that is, if --covar was provided, will combine --covar with PCs as a final covar file.
    --axes              <int>   [O]    Along with --PCA, set number of PCs to keep in the PCA of kept individuals by pruned variants, and add PCs as covar. (default: 20)
# OTHERs
    --constrain                 [O]     By default, LDAK will not restrict heritability estimates to be within [0,1]. This is generally our preference, as to obtain unbiased estimates of heritabilities near zero, it is necessary to accept the possibility of negative estimates. Set this option will force all heritability estimates within [0,1].
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
export threads=2
export gvcf=""
export pheno=""
export pheno_tag=""
export covar=""
export out="PANZ"
export mode="ALL"
export part_tsv=""
export part_bed=""
export ref_pre="00_ref"
export PCA="FALSE"
export keep="FALSE"
export axes=20
export window_prune=0.2
export window_kb=1000
export constrain="FALSE"
# export summary=""
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
        -v|--gvcf)
            #echo "set argument \"$1\" with value: $2" >&2
            gvcf=$2
            shift 2
        ;;
        -p|--pheno)
            #echo "set argument \"$1\" with value: $2" >&2
            pheno=$2
            shift 2
        ;;
        --pheno_tag)
            #echo "set argument \"$1\" with value: $2" >&2
            pheno_tag=$2
            shift 2
        ;;
        -c|--covar)
            #echo "set argument \"$1\" with value: $2" >&2
            covar=$2
            shift 2
        ;;
        -o|--out)
            #echo "set argument \"$1\" with value: $2" >&2
            out=$2
            shift 2
        ;;
        # --mode)
        #     #echo "set argument \"$1\" with value: $2" >&2
        #     mode=$2
        #     echo "$mode" | grep -P "ALL|PART" >/dev/null
        #     if [ $? -ne 0 ];then
        #         gst_err "Wrong mode value: $mode, should be one of ALL/PART";
        #         exit 1;
        #     fi
        #     shift 2
        # ;;
        --part_tsv)
            #echo "set argument \"$1\" with value: $2" >&2
            part_tsv=$2
            shift 2
        ;;
        --part_bed)
            #echo "set argument \"$1\" with value: $2" >&2
            part_bed=$2
            shift 2
        ;;
        --PCA)
            #echo "set argument \"$1\" with value: $2" >&2
            PCA="TRUE"
            shift 1
        ;;
        --axes)
            #echo "set argument \"$1\" with value: $2" >&2
            axes=$2
            shift 2
        ;;
        --keep)
            #echo "set argument \"$1\" with value: $2" >&2
            keep="TRUE"
            shift 1
        ;;
        --window_kb)
            #echo "set argument \"$1\" with value: $2" >&2
            window_kb=$2
            shift 2
        ;;
        --window_prune)
            #echo "set argument \"$1\" with value: $2" >&2
            window_prune=$2
            shift 2
        ;;
        --constrain)
            #echo "set argument \"$1\" with value: $2" >&2
            constrain="TRUE"
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
check_var_empty gvcf pheno covar out mode pheno_tag
check_var_numeric threads axes
check_abs_path $gvcf $pheno $covar
check_files_exists $gvcf $pheno $covar
gst_log "ALL Start ..."
if [ -n "$part_tsv" ];then
    check_files_exists $part_tsv
    check_abs_path $part_tsv
    gst_rcd "Found partition table file, will also calculate partitional h2 accordingly ..."
    mode="PART"
elif [ -n "$part_bed" ];then
    check_files_exists $part_bed
    check_abs_path $part_bed
    gst_rcd "Found partition region file, will also calculate partitional h2 accordingly ..."
    mode="PART"
fi

# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# USAGE:0 input1 input2
chr_cutw () {
    # cut weight and calc kin by chr
    local chr=$1
    local data=$ref_pre
    check_files_exists $data.bed $data.bim $data.fam
    mkdir -p  ${data}_chr
    gst_rcd "cut weights for $chr (check $out/${data}_chr/chr_$chr.log for details)..."
    # ? cut weights
    {
    if [ ! -s "${data}_chr/chr_$chr.cutw_done" ];then
        ldak --cut-weights ${data}_chr/chr_$chr --bfile ${data} --chr $chr &&\
        ldak --calc-weights-all ${data}_chr/chr_$chr --bfile ${data} --chr $chr
        if [ $? -ne 0 ];then gst_err "chr_kin for $chr failed: Non-zero exit"; exit 1;fi
        check_files_exists ${data}_chr/chr_$chr/weights.short
        echo "Done chr_cutw for $chr" >${data}_chr/chr_$chr.cutw_done;
    fi
    } 1> ${data}_chr/chr_$chr.log 2>&1
}
export -f chr_cutw

# USAGE:0 input1 input2
chr_kins () {
    # compute kinships for each chr
    local chr=$1
    local data=$ref_pre
    check_files_exists ${data}.all_weights $data.bed $data.bim $data.fam ${data}_kin/partition.list
    gst_rcd "calculate kinships for $chr ..."
    {
    if [ ! -s "${data}_kin/chr$chr.kins_done" ];then
        ldak --calc-kins ${data}_kin --bfile ${data} --partition $chr --weights ${data}.all_weights --power -0.25
        if [ $? -ne 0 ];then gst_err "chr_kins for $chr failed: Non-zero exit"; exit 1;fi
        echo "Done chr_kins for $chr" > ${data}_kin/chr$chr.kins_done
    fi
    } 1>>${data}_chr/chr_$chr.log 2>&1
}
export -f chr_kins

# USAGE:0 input1 input2
vcf2bed_kin () {
    #input gvcffile,output plink bed file, prune and keep and PCA if set, cut weights and calc all kinship
    # ? make bed
    local ingvcf=$1
    local data=$ref_pre
    gst_rcd "get plink bed ..."
    if [ -s "${data}.bed" -a -s "${data}.bim" -a -s "${data}.fam" ];then
        gst_warn "Found pre-exists bed files of ${data}. Skip running."
    else
        plink --vcf $ingvcf --make-bed --threads $threads --out ${data} 1>${data}.v2bed_log 2>&1
        if [ $? -ne 0 ];then gst_err "get ${data} bed failed: check ${data}.v2bed_log for details"; exit 1;fi
    fi

    gst_rcd "prune variants within $window_kb for correlation square < $window_prune ..."
    # ? prune
    if [ ! -s "${data}.prune_done" ];then
        ldak --thin ${data}_prune --bfile ${data} --window-kb $window_kb --window-prune $window_prune 1>${data}.prune_log 2>&1
        if [ $? -ne 0 ];then gst_err "prune failed: check ${data}.prune_log for details"; exit 1;fi
        echo "done" > ${data}.prune_done
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
    fi

    gst_rcd "Calculate direct kinships from pruned variants ..."
    # ? calc direct kins from prune.in variants
    if [ ! -s "${data}.driect_kins_done" ];then
        ldak --calc-kins-direct ${data}_prune --bfile ${data} --extract ${data}_prune.in --ignore-weights YES --power -0.25 1> ${data}.driect_kins_log 2>&1
        if [ $? -ne 0 ];then gst_err "keep failed: check ${data}.driect_kins_log for details"; exit 1;fi
        echo "done" > ${data}.driect_kins_done
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
    fi

    gst_rcd "Get covar ..."
    # ? if --PCA, merge provided covar with PCA
    if [ ! -s "${data}.covar_done" ];then
        if [ "$PCA" == "TRUE" ];then
            # ? run pca and join with covar
            ldak --pca ${data}_pca --grm 00_ref_prune --axes $axes 1>${data}.pca_log 2>&1
            if [ $? -ne 0 ];then gst_err "run pca failed: check ${data}.pca_log for details"; exit 1;fi
            # ? join pca with covar
            cat $covar | perl -lane '$,="\t";print $F[0],"$F[1],$F[2]",@F[3..$#F];' > ${data}.covar_tmp1 &&\
            cat ${data}_pca.vect | perl -lane '$,="\t";print $F[0],"$F[1],$F[2]",@F[3..$#F];' > ${data}.covar_tmp2 &&\
            csvtk join -TH -f 1 ${data}.covar_tmp1 ${data}.covar_tmp2 | sed 's/"//g' > ${data}.covar &&\
            rm -f ${data}.covar_tmp1 ${data}.covar_tmp2
            if [ $? -ne 0 ];then gst_err "Join covar and PCs failed: Non-zero exit"; exit 1;fi
        else
            # ? just cp covar
            cat $covar > ${data}.covar
        fi
        if [ $? -ne 0 ];then gst_err "Get covar failed: Non-zero exit"; exit 1;fi
        check_files_exists ${data}.covar
        echo "done" > ${data}.covar_done
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
    fi

    # get chr ids from ${data}.bim
    gst_rcd "get chrs ..."
    if [ ! -s "${data}.chr_list" ];then
        csvtk cut -tTH -f 1 -j $threads ${data}.bim | csvtk uniq -tTH -f 1 -j $threads -o ${data}.chr_list
        if [ $? -ne 0 ];then gst_err "get chr list failed: Non-zero exit";rm -f ${data}.chr_list; exit 1;fi
    fi

    # get weights
    gst_log "Get weights ..."
    if [ ! -s "${data}.weights_done" ];then
        # ? split cut weights by chr
        gst_rcd "Get each chr weights ..."
        parallel -j $threads chr_cutw :::: ${data}.chr_list
        if [ $? -ne 0 ];then gst_err "split run chr_cutw failed: Non-zero exit"; exit 1;fi
        # ? merge all weights
        gst_rcd "Merge weights ..."
        cat ${data}_chr/chr_*/weights.short > ${data}.all_weights
        if [ $? -ne 0 ];then gst_err "merge weights failed: Non-zero exit"; rm -f ${data}.all_weights; exit 1;fi
        check_files_exists ${data}.all_weights
        echo "Done" > ${data}.weights_done
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
    fi

    # get kinships
    gst_log "Get kinships ..."
    if [ ! -s "${data}.kinship_done" ];then
        # ? cut kins
        if [ ! -s "${data}_kin/partition.list" ]; then
            gst_rcd "Cut kins ..."
            ldak --cut-kins ${data}_kin --bfile ${data} --by-chr YES 1>${data}_cut_kin.log 2>&1
            if [ $? -ne 0 ];then gst_err "Cut kins by chr failed: Non-zero exit. Check ${data}_cut_kin.log for details"; rm -rf ${data}_kin;exit 1;fi
            rm -f ${data}_cut_kin.log
        fi
        # ? each chr kin
        parallel -j $threads chr_kins :::: ${data}.chr_list
        if [ $? -ne 0 ];then gst_err "each chr kin failed: Non-zero exit"; exit 1;fi
        # ? merge kins
        ldak --join-kins ${data}_kin 1>${data}_join_kin.log 2>&1
        if [ $? -ne 0 ];then gst_err "merge kinship failed: Non-zero exit. Check ${data}_join_kin.log for details.";rm -f ${data}_kin/kinships.all.grm*; exit 1;fi
        rm -f ${data}_join_kin.log
        echo "Done" > ${data}.kinship_done
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running." >&2
    fi
}
export -f vcf2bed_kin

# USAGE:0 input1 input2
sub_grm () {
    # sub set grm based on all bed and grm file, and extract VID file
    local extract=$1
    local data=$ref_pre
    local cur_outpre=${extract%.*}
    local kin=${data}_kin/kinships.all
    check_files_exists $extract ${data}.bed ${kin}.grm.bin
    if [[ ! -s "$cur_outpre.sub_kin.done" ]];then
        gst_rcd "subset grm for $(basename $cur_outpre) ..."
        ldak --bfile ${data} --sub-grm $cur_outpre.sub_kin --grm $kin --extract $extract 1>$cur_outpre.sub_kin.log 2>&1
        if [ $? -ne 0 ];then gst_err "get sub kin for $cur_outpre failed: Non-zero exit. Check $cur_outpre.sub_kin.log for details"; exit 1;fi
        echo "Done for sub kin of $cur_outpre" > $cur_outpre.sub_kin.done
    else
        gst_warn "Found pre-exist sub grm file for $(basename $cur_outpre). Skip."
    fi
}
export -f sub_grm

# USAGE:0 input1 input2
h2snp () {
    # get h2snp using LDAK
    local data=$ref_pre
    local grm=$1
    local pheno_pre=$pheno_tag
    local outpre=$(basename ${grm})
    local cur_pheno=$pheno
    local cur_covar=${data}.covar
    check_files_exists $grm.grm.bin $cur_pheno $cur_covar $data.bed
    mkdir -p 01_REML_$pheno_pre
    if [[ ! -s "01_REML_$pheno_pre/$outpre.h2_done" ]];then
        # ? keep only unrelated samples
        if [ "$keep" == "TRUE" ];then
            gst_rcd "Keep only unrelated individuals ..."
            if [ ! -s "01_REML_$pheno_pre/keep.keep_done" ];then
                ldak --filter 01_REML_$pheno_pre/keep --grm ${data}_prune --pheno $cur_pheno 1>01_REML_$pheno_pre/keep.keep_log 2>&1
                if [ $? -ne 0 ];then gst_err " CMD failed: check 01_REML_$pheno_pre/keep.keep_log for details"; exit 1;fi
                check_files_exists 01_REML_$pheno_pre/keep.keep
                echo "done" > 01_REML_$pheno_pre/keep.keep_done
            fi

            gst_rcd "Calculate h2 for $outpre ..."
            if [ "$constrain" == "TRUE" ];then
                # ? constrain h2 within [0,1]
                ldak --reml 01_REML_$pheno_pre/$outpre --grm $grm --pheno $cur_pheno --covar $cur_covar --bfile $data --keep 01_REML_$pheno_pre/keep.keep --constrain YES 1>01_REML_$pheno_pre/$outpre.log 2>&1
            else
                # ? do not constrain
                ldak --reml 01_REML_$pheno_pre/$outpre --grm $grm --pheno $cur_pheno --covar $cur_covar --bfile $data --keep 01_REML_$pheno_pre/keep.keep 1>01_REML_$pheno_pre/$outpre.log 2>&1
            fi
            if [ $? -ne 0 ];then gst_err "Calc h2 for $grm failed: Non-zero exit. Check 01_REML_$pheno_pre/$outpre.log for details";rm -f 01_REML_$pheno_pre/${outpre}.reml; exit 1;fi
            echo "Done for h2 of $grm" > 01_REML_$pheno_pre/$outpre.h2_done
        else
            gst_rcd "Use all provided individuals ..."
            gst_rcd "Calculate h2 for $outpre ..."
            if [ "$constrain" == "TRUE" ];then
                # ? constrain h2 within [0,1]
                gst_rcd "Force h2 within [0,1] ..."
                ldak --reml 01_REML_$pheno_pre/$outpre --grm $grm --pheno $cur_pheno --covar $cur_covar --bfile $data --constrain YES 1>01_REML_$pheno_pre/$outpre.log 2>&1
            else
                # ? do not constrain
                gst_rcd "h2 may out range of [0,1] ..."
                ldak --reml 01_REML_$pheno_pre/$outpre --grm $grm --pheno $cur_pheno --covar $cur_covar --bfile $data 1>01_REML_$pheno_pre/$outpre.log 2>&1
            fi
            if [ $? -ne 0 ];then gst_err "Calc h2 for $grm failed: Non-zero exit. Check 01_REML_$pheno_pre/$outpre.log for details";rm -f 01_REML_$pheno_pre/${outpre}.reml; exit 1;fi
            echo "Done for h2 of $grm" > 01_REML_$pheno_pre/$outpre.h2_done
        fi
    else
        gst_warn "Found pre-exist h2 results for $outpre . Skip."
    fi
}
export -f h2snp

# USAGE:0 feature_pre
calc_kinship_each_feature () {
    # ? get each feature grm
    local cur_wd=$PWD
    local data=$ref_pre
    local feature_pre=$1
    gst_log "Get grm for each feature of $feature_pre ..."
    # cd ${data}_Feature_${feature_pre}
    if [ $? -ne 0 ];then gst_err "set work dir failed: Non-zero exit"; exit 1;fi
    # ? get each feature kinship
    parallel -j $threads -q -k sub_grm :::: <(ls ${data}_Feature_${feature_pre}/*.vid)
    if [ $? -ne 0 ];then gst_err "Paraellel run sub_grm failed: Non-zero exit"; exit 1;fi
}
export -f calc_kinship_each_feature

# get each feature grm
# USAGE:0 input1 input2
parse_feature_tsv () {
    local cur_outpre=""
    local data=$ref_pre
    local cur_tag=""
    local cur_clean_file=$(basename ${part_tsv}.clean)
    check_var_empty cur_clean_file
    gst_log "Dealing with feature table $(basename $part_tsv) ..."
    # ? check if feature table was correctly headed
    head -n 1 $part_tsv | grep -P "^VID\t" >/dev/null
    if [ $? -ne 0 ];then gst_err "Wrong feature table format of $part_tsv!"; exit 1;fi
    # ? make sure only VIDs in ${data}.bim was included
    # if [ ! -s "${data}.Check_feature_tsv_done" ];then
    #     gst_rcd "Check VIDs of feaure table ..."
    #     csvtk join -tT -j $threads -f 1 <(csvtk cut -tTH -f 2 -j $threads ${data}.bim | sed '1i\VID') $part_tsv -o $cur_clean_file
    #     if [ $? -ne 0 ];then gst_err "Clean feature table failed: Non-zero exit";rm -f $cur_clean_file; exit 1;fi
    #     check_files_exists $cur_clean_file
    #     echo "Done" > ${data}.Check_feature_tsv_done
    # fi
    ln -s $part_tsv $cur_clean_file
    csvtk head -j $threads -tTH -n 1 $cur_clean_file | csvtk transpose -tTH | sed '1d' | while read cur_tag
    do
        gst_log "Getting features of $cur_tag ..."
        cur_outpre=$cur_tag
        if [[ ! -s "${data}_Feature_${cur_outpre}/parse_$cur_outpre.done" ]];then
            mkdir -p ${data}_Feature_$cur_outpre
            csvtk cut -tT -j $threads -f "VID,$cur_tag" $cur_clean_file| sed '1d' > ${cur_tag}.feature_tmp
            if [ $? -ne 0 ];then gst_err "get ${cur_tag}.feature_tmp failed: Non-zero exit";rm -f ${cur_tag}.feature_tmp; exit 1;fi
            csvtk cut -j $threads -tTH -f 2 ${cur_tag}.feature_tmp |\
            csvtk uniq -tTH -f 1 -j $threads | grep -v -P "^NA$" |\
            while read cur_feature
            do
                gst_rcd "$cur_feature ..."
                csvtk grep -tTH -f 2 -p "$cur_feature" -j $threads ${cur_tag}.feature_tmp | csvtk cut -tTH -f 1 -j $threads |csvtk uniq -tTH -f 1 -j $threads -o ${data}_Feature_${cur_outpre}/${cur_tag}_${cur_feature}.vid
                if [ $? -ne 0 ];then gst_err "get vid for $cur_feature failed: Non-zero exit";rm -f ${data}_Feature_${cur_outpre}/${cur_tag}_${cur_feature}.vid; exit 1;fi
                local cur_feature_num=0
                cur_feature_num=$(cat ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid | wc -l)
                if [[ "$cur_feature_num" -lt "50" ]];then
                    gst_warn "Feature ${cur_outpre}_${cur_feature} has only $cur_feature_num variants (Ignore this if it's a normal situation). "
                fi
                if [[ "$cur_feature_num" -eq "0" ]];then
                    gst_warn "Feature ${cur_outpre}_${cur_feature} has NO variants, skipping this feature. "
                    rm -f ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid
                fi
            done
            # get all done mark file
            echo "Done parse feature of $cur_outpre" > ${data}_Feature_${cur_outpre}/parse_$cur_outpre.done
        else
            gst_warn "Found pre-exist feature files for $cur_outpre. Skip generating new. This may CAUSE PROBLEMS if you have changed the content of the input feature file but kept the name as is. Remove \"$out/${data}_Feature_${cur_outpre}/parse_$cur_outpre.done\" to force generating new files."
        fi
        # ? calc kinship for each feature
        calc_kinship_each_feature $cur_outpre
        # ? clean tmp files
        rm -f ${cur_tag}.feature_tmp
    done
}
export -f parse_feature_tsv

# USAGE:0 input1 input2
parse_feature_bed () {
    local cur_outpre=$(basename ${part_bed%.*})
    local data=$ref_pre
    check_files_exists $data.bim $part_bed
    gst_log "Dealing with feature region $cur_outpre ... "
    mkdir -p ${data}_Feature_$cur_outpre
    # ? convert bim to bed format
    if [ ! -s "${data}_Feature_$cur_outpre/${data}.bed" ];then
        cat $data.bim | parallel --pipe -j $threads -q -k perl -lane 'BEGIN{$,="\t";} print $F[0],$F[3]-1,$F[3],$F[1];' > ${data}_Feature_$cur_outpre/${data}.bed
        if [ $? -ne 0 ];then gst_err "bim2bed failed: Non-zero exit";rm -f ${data}_Feature_$cur_outpre/${data}.bed; exit 1;fi
    fi
    # ? get each feature vid
    if [[ ! -s "${data}_Feature_${cur_outpre}/parse_$cur_outpre.done" ]];then
        csvtk cut -j $threads -tTH -f 2 $part_bed| csvtk uniq -tTH -f 1 -j $threads | while read cur_feature
        do
            gst_rcd "$cur_feature ..."
            bedtools intersect -a ${data}_Feature_$cur_outpre/${data}.bed -b <(csvtk grep -tTH -f 4 -p "$cur_feature" -j $threads $part_bed) -wa | csvtk cut -tTH -f 4 -j $threads |csvtk uniq -tTH -j $threads -f 1 -o ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid
            if [ $? -ne 0 ];then gst_err "get vid for $cur_feature failed: Non-zero exit";rm -f ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid; exit 1;fi
            local cur_feature_num=0
            cur_feature_num=$(cat ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid | wc -l)
            if [[ "$cur_feature_num" -lt "50" ]];then
                gst_warn "Feature ${cur_outpre}_${cur_feature} has only $cur_feature_num variants (Ignore this if it's a normal situation). "
            fi
            if [[ "$cur_feature_num" -eq "0" ]];then
                gst_warn "Feature ${cur_outpre}_${cur_feature} has NO variants, skipping this feature. "
                rm -f ${data}_Feature_${cur_outpre}/${cur_outpre}_${cur_feature}.vid
            fi
        done
        # get all done mark file
        echo "Done parse feature of $cur_outpre" > ${data}_Feature_${cur_outpre}/parse_$cur_outpre.done
    else
        gst_warn "Found pre-exist feature files for $cur_outpre. Skip generating new. This may CAUSE PROBLEMS if you have changed the content of the input feature file but kept the name as is. Remove \"$out/${data}_Feature_${cur_outpre}/parse_$cur_outpre.done\" to force generating new files."
    fi
    # ? calc kinship for each feature
    calc_kinship_each_feature $cur_outpre
}
export -f parse_feature_bed

# USAGE:0 ldak.reml.file
fmt_reml_out () {
    # format ldak reml output into one line
    local reml_file=$1
    local cur_reml_prefix=$(basename ${reml_file%.reml} | perl -pe 's/kinships\.all/ALL/g;s/\.sub_kin//g;')
    check_var_empty cur_reml_prefix
    export cur_reml_prefix
    # ? fmt to one line with perl
    cat $reml_file | perl -lane '
        BEGIN{$,="\t";}
        if($#F<=1){
            $hash{$F[0]}=$F[1];
        }
        if(/^Her_ALL/){
            $info=join("\t",$ENV{cur_reml_prefix},@F[1..$#F]);
        }
        END{
            # ? header
            # print "Trait","Feature","Heritability","Her_SD","Size","Mega_Intensity","Int_SD","Covar_Heritability","Total_Samples","With_Phenotypes","Null_Likelihood","Alt_Likelihood","LRT_Stat","LRT_P";
            # ? main
            print $ENV{pheno_tag},$info,$hash{"Covar_Heritability"},$hash{"Total_Samples"},$hash{"With_Phenotypes"},$hash{"Null_Likelihood"},$hash{"Alt_Likelihood"},$hash{"LRT_Stat"},$hash{"LRT_P"};
        }'
}
export -f fmt_reml_out


# >>>>>>>>>>>>>>>>>>>>>>>> MAIN >>>>>>>>>>>>>>>>>>>>>>>>
data=$ref_pre
mkdir -p $out
# ? clean input vcf to make sure uniq VIDs and chrs are all numbers
gst_log "Cleaning INPUT variants ..."
if [ ! -s "$out/clean.done" ];then
    bcftools view -h --threads $threads $gvcf -o $out/${data}.vcf
    if [ $? -ne 0 ];then gst_err "get header failed: Non-zero exit"; exit 1;fi
    csvtk uniq -tTH -f 3 -j $threads $gvcf | parallel --pipe -j $threads -q -k perl -F"\t" -lane '
    BEGIN{$,="\t"}
    print if $F[0] =~ /^\d+$/;' >> $out/${data}.vcf
    if [ $? -ne 0 ];then gst_err "uniq vid failed: Non-zero exit"; exit 1;fi
    bgzip -@ $threads $out/${data}.vcf
    if [ $? -ne 0 ];then gst_err "bgzip $out/${data}.vcf failed: Non-zero exit"; exit 1;fi
    check_files_exists $out/${data}.vcf.gz
    echo "Clean done" > $out/clean.done
fi

# ! change work dir to $out
cd $out
if [ $? -ne 0 ];then gst_err "link $gvcf and cd to $out failed: Non-zero exit"; exit 1;fi

gst_log "Get kinships ..."
# ? get kinships
if [ ! -s "kinship.done" ];then
    gst_rcd "Current work dir: $PWD"
    gst_log "Format variants ..."
    vcf2bed_kin ${data}.vcf.gz
    check_files_exists $data.bed $data.bim $data.fam ${data}_kin/kinships.all.grm.bin

    if [ "$mode" == "PART" ];then
        gst_log "Calculate each feature kinships ..."
        # ? get each feature grm
        if [ -s "$part_bed" ];then
            parse_feature_bed
        fi
        if [ -s "$part_tsv" ];then
            parse_feature_tsv
        fi
        gst_log "Calculate each feature variant heritability ..."
        # ? get each feature kinship list
        ls ${data}_Feature_*/*.grm.bin | sed 's/.grm.bin//' > ${data}_All_Feature.list
        if [ $? -ne 0 ];then gst_err "get All Feature list failed: Non-zero exit";rm -f ${data}_All_Feature.list; exit 1;fi
    fi
    echo "Done " > kinship.done
else
    gst_warn "Skip kinship step. Use pre-exist results ..."
fi

# ? calc h2
if [ ! -s "02_h2_${pheno_tag}.done" ];then
    # ? calc h2 for all
    gst_log "Calculate All variant heritability ($PWD)..."
    h2snp ${data}_kin/kinships.all
    if [ $? -ne 0 ];then gst_err "Calc All h2 failed: Non-zero exit"; exit 1;fi
    # ? Calc h2 for each feature
    if [ "$mode" == "PART" ];then
        gst_log "Calculate Each feature variant heritablity ..."
        check_files_exists ${data}_All_Feature.list
        parallel -j $threads -q -k h2snp :::: ${data}_All_Feature.list
        if [ $? -ne 0 ];then gst_err "Parallel run h2snp failed: Non-zero exit"; exit 1;fi
    fi
    # ? merge results
    gst_log "Get merged h2 results to 02_h2_${pheno_tag}.tsv"
    ls 01_REML_${pheno_tag}/*.reml > 01_REML_${pheno_tag}.all_her_list
    check_files_exists 01_REML_${pheno_tag}.all_her_list
    perl -e '$,="\t";print "Trait","Feature","Heritability","Her_SD","Size","Mega_Intensity","Int_SD","Covar_Heritability","Total_Samples","With_Phenotypes","Null_Likelihood","Alt_Likelihood","LRT_Stat","LRT_P\n";' > 02_h2_${pheno_tag}.tsv &&\
    parallel -j $threads -q -k fmt_reml_out :::: 01_REML_${pheno_tag}.all_her_list 1>>02_h2_${pheno_tag}.tsv
    if [ $? -ne 0 ];then gst_err "Merge h2 output for $pheno_tag failed: Non-zero exit"; rm -f 02_h2_${pheno_tag}.tsv; exit 1;fi
    echo "Done !" > 02_h2_${pheno_tag}.done

    # ? clean temp result
    rm -rf 01_REML_${pheno_tag} 01_REML_${pheno_tag}.all_her_list
    if [ $? -ne 0 ];then gst_err "Remove 01_REML_${pheno_tag} failed: Non-zero exit"; exit 1;fi

    gst_log "All Done ! Final heritability result in
    $(ls $PWD/02_h2_${pheno_tag}.tsv)"
else
    gst_warn "Found pre-exist output for ${pheno_tag}. Skip running."
fi
# <<<<<<<<<<<<<<<<<<<<<<<< MAIN <<<<<<<<<<<<<<<<<<<<<<<<



