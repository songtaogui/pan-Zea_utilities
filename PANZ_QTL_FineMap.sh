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

export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix check_files_exists check_abs_path
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<
usage=$(
cat <<EOF
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
EOF
)

if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export threads=2
export out=PANZ_FineMap_Out
export gvcf=""
export trait_file=""
export covar_file=""
export trait=""
export control=""
export gwas=""
export info_cols=""
export min_peak_logP=5
export inlimit=0.001
export flank=10000
export min_count=3
export ld_control=0.25
export all="False"
export credible_prob=0.95
export tissue="kernel"
export tar_value="FALSE"
export heritability="FALSE"
# export dump="NAN"
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
        -o|--out)
            #echo "set argument \"$1\" with value: $2" >&2
            out=$2
            shift 2
        ;;
        -v|--gvcf)
            #echo "set argument \"$1\" with value: $2" >&2
            gvcf=$2
            shift 2
        ;;
        --trait)
            #echo "set argument \"$1\" with value: $2" >&2
            trait_file=$2
            shift 2
        ;;
        --covar)
            #echo "set argument \"$1\" with value: $2" >&2
            covar_file=$2
            shift 2
        ;;
        --gwas)
            #echo "set argument \"$1\" with value: $2" >&2
            gwas=$2
            shift 2
        ;;
        --info_cols)
            #echo "set argument \"$1\" with value: $2" >&2
            info_cols=$2
            shift 2
        ;;
        --min_peak_logP)
            #echo "set argument \"$1\" with value: $2" >&2
            min_peak_logP=$2
            shift 2
        ;;
        --flank)
            #echo "set argument \"$1\" with value: $2" >&2
            flank=$2
            shift 2
        ;;
        --inlimit)
            #echo "set argument \"$1\" with value: $2" >&2
            inlimit=$2
            shift 2
        ;;
        --ld_control)
            #echo "set argument \"$1\" with value: $2" >&2
            ld_control=$2
            shift 2
        ;;
        --credible_prob)
            #echo "set argument \"$1\" with value: $2" >&2
            credible_prob=$2
            shift 2
        ;;
        --tissue)
            #echo "set argument \"$1\" with value: $2" >&2
            tissue=$2
            shift 2
        ;;
        --tar)
            #echo "set argument \"$1\" with value: $2" >&2
            tar_value="TRUE"
            shift 1
        ;;
        --heritability)
            #echo "set argument \"$1\" with value: $2" >&2
            heritability="TRUE"
            shift 1
        ;;
        # --dump)
        #     #echo "set argument \"$1\" with value: $2" >&2
        #     dump=$2
        #     shift 2
        # ;;
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
check_files_exists $gvcf $gwas $covar_file
if [ ! -f "${trait_file}" ];then
    gst_err "No file: ${trait_file}"
    exit 1
fi
trait=$(basename ${trait_file%.*})
check_var_empty gvcf trait gwas info_cols covar_file flank tissue
check_var_numeric threads flank inlimit min_peak_logP min_count ld_control credible_prob
check_sftw_path summarize_dap2enloc.pl
# check output dir exists
if [ -d "$out" ];then
    gst_warn "output dir $out already exists !"
    # exit 1
else
    gst_log "Creat output dir $out ..."
fi
mkdir -p $out/99_Misc
if [ $? -ne 0 ];then gst_err "mkdir $out/99_Misc failed: Non-zero exit"; exit 1;fi

# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> Functions >>>>>>>>>>>>>>>>>>>>>>>>

# USAGE:0 
calc_zscore () {
    # calc_zscore from input gwas info file in format:
    # <variant_ID> <Peak_ID> <chr> <start> <end> <P-value> <effect_size>
    local in=$1
    local out=$2
    check_files_exists $in
    Rscript -e '
        argv=as.character(commandArgs(TRUE));
        tt=read.table(argv[1],sep="\t",header=F,comment.char="$")
        tt$V8=sign(tt$V7)*abs(qnorm(tt$V6/2))
        write.table(tt,file=argv[2],sep="\t",quote=F,row.names=F,col.names=F)
    ' $in $out
    if [ $? -ne 0 ];then gst_err "calc_zscore failed: Non-zero exit";
        #rm -f $out;
        exit 1;
    fi
    check_files_exists $out
}
export -f calc_zscore

# USAGE:0
run_harvester () {
    # run harvester to get the gwas peaks
    if [ ! -s "$out/99_Misc/fmt_gwas.tsv" ];then
        csvtk cut -tT -f "$info_cols" $gwas -j $threads | sed '1d' | parallel --pipe -j $threads -k -q perl -F"\t" -lane '
            BEGIN{$,="\t";}
            ($vid,$chr,$posi0,$posi,$p,$beta)=@F;
            # $posi0=$posi;$posi++;
            $p=sprintf("%.4e",$p);
            # ? dealing with TRA
            if($p <= $ENV{inlimit}){
                print $chr,$posi0,$posi,$p,$beta,$vid;
            }
        ' | csvtk sort -tTH -k 1 -k 2:n -k 3:n -j $threads -o $out/99_Misc/fmt_gwas.tsv
        if [ $? -ne 0 ];then gst_err "fmt_gwas failed: Non-zero exit";rm -f $out/99_Misc/fmt_gwas.tsv; exit 1;fi
        check_files_exists $out/99_Misc/fmt_gwas.tsv
    fi
    # chr posi0 posi1 P beta variantID
    if [ ! -s "$out/99_Misc/HARVESTER_peak.bed" ];then
        # ? raw HARVESTER scripts
        if [ ! -s "$out/99_Misc/HARVESTER_out.tsv" ];then
            HARVESTER -chrcolumn 1 -lcolumn 2 -pcolumn 4 -header no -file $out/99_Misc/fmt_gwas.tsv -out $out/99_Misc/HARVESTER_out.tsv -delim tab -missing NA -inlimit $inlimit -peak-limit $min_peak_logP -dots 5 -shrink 2
            if [ $? -ne 0 ];then gst_err "HARVESTER failed: Non-zero exit";rm -f $out/99_Misc/HARVESTER_out.tsv; exit 1;fi
        fi
        # ? check if there are no gwas peaks
        local peak_no=$(sed '1d' $out/99_Misc/HARVESTER_out.tsv | wc -l )
        if [ $peak_no -lt 1 ];then
            gst_warn "No GWAS PEAKs found ! Stop running"
            exit 0
        fi
        # ? parse HARVESTER_out, flanking and merge to get all peaks
        sed '1d' $out/99_Misc/HARVESTER_out.tsv | perl -F"\t" -lane '
            BEGIN{$,="\t";}
            $chr=$F[1];
            ($s1,$e)=split(/-/,$F[2]);
            $s=$s1;$e++;
            print $chr,$s,$e,@F[3,11,12,13,19];
        ' | sortBed -i stdin | perl -F"\t" -lane '
            BEGIN{$,="\t";}
            $peak_range=sprintf("%s:%s-%s",$F[0],$F[1],$F[2]);
            # $peak_id=sprintf("Peak%04d",$.);
            print @F,$peak_range;
            # ?format: chr start end max range count spacing GQS peak_range
        ' | perl -F"\t" -lane '
            BEGIN{$,="\t";}
            # ? flanking
            $F[1]-=$ENV{flank};
            $F[2]+=$ENV{flank};
            $F[1]=0 if $F[1] < 0;
            print @F;
        ' | mergeBed -c 4,9,6,7,8 -o max,collapse,sum,sum,mean | perl -F"\t" -lane '
            # add peak id
            BEGIN{$,="\t";}
            $peak_id=sprintf("Peak%04d",$.);
            print @F,$peak_id;
        ' > $out/99_Misc/HARVESTER_peak.bed
        # ? chr start end max Raw_range count spacing GQS Peak_ID
        if [ $? -ne 0 ];then gst_err "parse HARVESTER out failed: Non-zero exit"; exit 1;fi
            ## harvester out:
            #  1  filename
            #  2  chrom
            #  3  range
            #  4  max
            #  5  bestslope0
            #  6  bestslope
            #  7  monot
            #  8  multip
            #  9  reps
            # 10  ratio
            # 11  kolmo
            # 12  range
            # 13  count
            # 14  spacing
            # 15  balance
            # 16  skew
            # 17  vbal1
            # 18  vbal2
            # 19  maxbyMean
            # 20  GQS
    fi

    # ? get variant IDs in the peak regions
    #  1    2     3   4  5     6        7    8    9   10   11   12    13     14    15      16
    # chr posi0 posi1 P beta variantID chr start end max range count spacing GQS peak_id  OVLP_bp
    if [ ! -s "$out/99_Misc/In_Peak_filtered_variant.id" ];then
        bedtools intersect -a $out/99_Misc/fmt_gwas.tsv -b $out/99_Misc/HARVESTER_peak.bed -wo > $out/99_Misc/In_Peak_filtered_variant.tsv &&\
        csvtk cut -tTH -f 6,15 -j $threads -o $out/99_Misc/In_Peak_filtered_variant.id $out/99_Misc/In_Peak_filtered_variant.tsv
        if [ $? -ne 0 ];then gst_err "get in-peak-variant-id failed: Non-zero exit"; exit 1;fi
    fi
    check_files_exists $out/99_Misc/In_Peak_filtered_variant.id $out/99_Misc/In_Peak_filtered_variant.tsv
    gst_log "Subset gwas signals by peaks ..."
    # csvtk join -tTH -f "1;6" $out/99_Misc/In_Peak_filtered_variant.id $out/99_Misc/fmt_gwas.tsv -j $threads -o $out/99_Misc/In_Peak_gwas.tsv
    csvtk cut -tTH -f "6,15,1,2,3,4,5" -j $threads -o $out/99_Misc/In_Peak_gwas.tsv $out/99_Misc/In_Peak_filtered_variant.tsv
    # vid peakid chr s e P beta
    if [ $? -ne 0 ];then gst_err "get in_peak_gwas ffailed: Non-zero exit";rm -f $out/99_Misc/In_Peak_gwas.tsv; exit 1;fi
}
export -f run_harvester

# USAGE:0 peakname outpre
get_z_ld () {
    # generate z-score and ld matrix based on HARVESTER out and gwas input
    local input_gwas=$out/99_Misc/In_Peak_gwas.tsv
    local region=$1
    local peakname=$2
    local outpre=$out/${peakname}
    # vid peakid chr s e P beta
    check_files_exists $input_gwas
    check_var_empty peakname outpre
    csvtk grep -tTH -p "$peakname" -f 2 -j $threads $out/99_Misc/In_Peak_filtered_variant.id | cut -f 1 > ${outpre}_vid.tmp
    # # ? make sure the zscore and the LD vcf have same set and order of variants
    # csvtk grep -tTH -p "$peakname" -f 2 -j $threads $out/99_Misc/In_Peak_filtered_variant.id | cut -f 1 > ${outpre}_zvid.tmp &&\
    # csvtk grep -tTH -p "$peakname" -f 2 -j $threads $input_gwas | cut -f 1 > ${outpre}_gvid.tmp &&\
    # cat ${outpre}_zvid.tmp ${outpre}_gvid.tmp | perl -lane 'print if++$h{$_}>1' > ${outpre}_vid.tmp

    # ? calc z-score
    gst_log "Calculate Z-scores ..."
    if [ ! -s "${outpre}_gwas_addz.tsv.raw" ];then
        csvtk grep -tTH -p "$peakname" -f 2 -j $threads $input_gwas | csvtk grep -tTH -f 1 -P ${outpre}_vid.tmp -j $threads -o ${outpre}_gwas.tmp &&\
        calc_zscore ${outpre}_gwas.tmp ${outpre}_gwas_addz.tsv.raw
        check_files_exists ${outpre}_gwas_addz.tsv.raw
        # rm -f ${outpre}_gwas.tmp
    fi
    region=$( cat ${outpre}_gwas_addz.tsv.raw | perl -lane '($c,$s)=@F[2,3] if $.==1; $e=$F[4] if eof; END{print "$c:$s-$e"}' )
    check_var_empty region
    # ? get subset vcf
    gst_log "get subset vcf (${region}) ..."
    if [ ! -s "${outpre}.vcf.gz" ];then
        check_files_exists $out/99_Misc/In_Peak_filtered_variant.id
        bcftools view -r ${region} $gvcf --threads $threads |\
            bcftools view --threads $threads -i "ID=@${outpre}_vid.tmp" |\
            perl -lane 'print if ++$h{$_}==1' | bcftools view --threads $threads -O z -o ${outpre}.vcf.gz
        if [ $? -ne 0 ];then gst_err "get subset vcf failed: Non-zero exit";rm -f ${outpre}.vcf.gz; exit 1;fi
        # rm -f ${outpre}_vid.tmp
    fi
    gst_log "Get LD matrix (r) ..."
    if [ ! -s "${outpre}.ld" ];then
        plink --vcf ${outpre}.vcf.gz --r --threads $threads --matrix --out ${outpre} &&\
        sed -i 's/nan/0/g' ${outpre}.ld
        if [ $? -ne 0 ];then gst_err "Get LD matrix failed: Non-zero exit";rm -f ${outpre}.ld; exit 1;fi
    fi
    gst_log "Get zscore tsv ..."
    if [ ! -s "${outpre}.zscore" ];then
        bcftools view ${outpre}.vcf.gz -H | cut -f 3 > ${outpre}_vcf_vid.tmp
        csvtk join -tTH -f 1 -j $threads ${outpre}_vcf_vid.tmp <( perl -lane 'print if ++$h{$_}==1' ${outpre}_gwas_addz.tsv.raw ) -o ${outpre}_gwas_addz.tsv &&\
        csvtk cut -tTH -f 1,8 -j $threads ${outpre}_gwas_addz.tsv -o ${outpre}.zscore
        if [ $? -ne 0 ];then gst_err "get addz gwas file failed: Non-zero exit"; rm -f ${outpre}_gwas_addz.tsv ${outpre}.zscore; exit 1;fi
    fi

}
export -f get_z_ld

# USAGE:0
get_credible_set_script () {
    # generate get_credible_set.pl in the output dir
cat << \EOF >$out/99_Misc/get_credible_set.pl
#! /usr/bin/env perl
$file = "";
$prob = 0.95;
for $i  (0..$#ARGV){
    if ($ARGV[$i] eq "-d"){
        $file = $ARGV[++$i];
        next;
    }
    if($ARGV[$i] eq "-p"){
        $prob = $ARGV[++$i];
        next;
    }
    if($ARGV[$i] eq "-g"){
        $gene = $ARGV[++$i];
        next;
    }
}

if($file eq ""){
    print STDERR "Error: dap-g output file is not specified.\n";
    exit(1);
}
open FILE, "grep \\\{ $file | ";
while(<FILE>){
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    next if $data[2]< $prob;
    $data[0] =~/(\d+)/;
    $sig{$1} = {};
    $cum{$1} = 0;
}
if (scalar(keys %sig)==0){
    printf STDERR "No %d%% credible set can be constructed\n", int($prob*100);
    exit;
}
open FILE, "grep \\\(\\\( $file | ";
while(<FILE>){
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    if(defined($sig{$data[4]}) && $cum{$data[4]}< $prob){
        $sig{$data[4]}->{$data[1]} = $data[2];
        $cum{$data[4]} += $data[2];
    }
}
$,="\t";
print "Trait\tSignal_Cluster\tVariant_ID\tPIP\n";
foreach $c (sort {$a <=> $b} keys %sig){
    foreach $snp (sort {$sig{$c}->{$b} <=> $sig{$c}->{$a}} keys %{$sig{$c}}){
        printf "%s\t%d\t%s\t%e\n",$gene,$c, $snp, $sig{$c}->{$snp}
    }
}
EOF

if [ $? -ne 0 ];then gst_err "get credible set script failed: Non-zero exit"; exit 1;fi

}
export -f get_credible_set_script

# USAGE: $0 peakname
run_dapg_zld () {
    # run dapg with z and ld
    local peakname=$1
    local zfile=$out/${peakname}.zscore
    local ldfile=$out/${peakname}.ld
    local outpre=$out/${peakname}
    # ? get trait name
    local tt_name=$trait
    if [ ! -s "${outpre}_dap-g.credible_set" ];then
        check_files_exists $zfile $ldfile
        check_var_empty outpre tt_name

        dap-g -d_z $zfile -d_ld $ldfile -t $threads -ld_control $ld_control --output_all -o ${outpre}_dap-g.out
        # ! dap-g exit status was always non-zero
        # if [ $? -ne 0 ];then gst_err "dap-g failed: Non-zero exit";rm -f ${outpre}_dap-g.out; exit 1;fi
        gst_log "get credible set ..."
        check_var_empty tt_name
        # get_credible_set_script
        check_files_exists $out/99_Misc/get_credible_set.pl ${outpre}_dap-g.out
        perl $out/99_Misc/get_credible_set.pl -g "${tt_name}-${peakname}" -p $credible_prob -d ${outpre}_dap-g.out > ${outpre}_dap-g.credible_set
        if [ $? -ne 0 ];then gst_err "get credible set failed: Non-zero exit";rm -f ${outpre}_dap-g.credible_set; exit 1;fi
    fi
    gst_log "DONE!"
}
export -f run_dapg_zld

get_anno_pip () {
    # get anno and pip, input dapg.out
    local dapg=$1
    export dapg
    local peak=$(echo $dapg | perl -pe 's/.*Each_Peak\/(.*)_dap-g.out/$1/')
    export peak
    local vcf=${dapg%%_dap-g.out}.vcf.gz
    # local trait=$(echo $dapg | perl -pe 's/.*\/(.*?)\/Each_Peak.*$/$1/;')
    local gene=$out/${trait}_${peak}
    if [[ -s "$dapg" ]];then
        mkdir -p ${gene}_tmp_dir
        # ? get anno
        if [[ ! -s "${gene}_anno.vcf" ]];then
            cat $dapg >${gene}_tmp_dir/${trait}_${peak}.dapg &&\
            summarize_dap2enloc.pl -dir ${gene}_tmp_dir -vcf $vcf -tissue $tissue > ${gene}_anno.vcf
            if [ $? -ne 0 ];then echo "anno $gene failed: Non-zero exit"; exit 1;fi
            if [ ! -s "${gene}_anno.vcf" ];then
                rm -f ${gene}_anno.vcf
            fi
            rm -rf ${gene}_tmp_dir
        fi
        # ? get pip
        if [[ ! -s "${gene}_PIPs.pip" ]];then
            grep "((" $dapg | perl -lane '
                BEGIN{
                    open(IN,"grep \"{\" $ENV{dapg} |") or die("Cannot open file: $inputfile");
                    while(<IN>){
                        chomp;s/^\s+//;
                        # $hash{$_}=1;

                        ($sigID,$sigPIP,$sigR2)=(split(/\s+/,$_))[0,1,2];
                        $sigID=~s/\D//g;$sigID=sprintf("sig%04d",$sigID);
                        $hash{$sigID}="$sigPIP\_$sigR2";
                    }
                }
                $sid=sprintf("sig%04d",$F[4]);
                print "$F[1]\t$ENV{peak}_$sid\_$hash{$sid}\t$F[2]" if $F[4] > 0 && $hash{$sid};
            ' >${gene}_PIPs.pip
            if [ $? -ne 0 ];then gst_err "get pip for $gene failed: Non-zero exit"; exit 1;rm -f ${gene}_PIPs.pip;fi
            if [ ! -s "${gene}_PIPs.pip" ];then
                rm -f ${gene}_PIPs.pip
            fi
        fi
    fi
}
export -f get_anno_pip

# USAGE: $0 vcf trait outpre includeID covar
vcf2hsq () {
    # input vcf file, output hsq
    local vcf=$1
    local trait=$2
    local outpre=$3
    local includeID=$4
    local covar=$5
    check_files_exists $vcf $trait $includeID $covar
    # if [ ! -s "${outpre}_hsq.out" ];then
    rm -f ${outpre}_hsq.out
        if [[ ! -s "${outpre}.bed" ]];then
            plink --vcf $vcf --make-bed --extract $includeID --out ${outpre} 1>${outpre}.plinklog 2>&1
            if [ $? -ne 0 ];then gst_err "get plink bed failed: Non-zero exit";rm -f ${outpre}.{bed,bim,fam,nosex,log}; exit 1;fi
        fi
        if [[ ! -s "${outpre}.grm.bin" ]];then
            gcta64 --bfile ${outpre} --make-grm --out ${outpre} 1>${outpre}.grmlog 2>&1
            if [ $? -ne 0 ];then gst_err "get grm failed: Non-zero exit";rm -f ${outpre}.grm.*; exit 1;fi
        fi
        if [[ ! -s "${outpre}.hsq" ]];then
            gcta64 --grm ${outpre} --pheno $trait --reml --reml-maxit 100 --qcovar $covar --out ${outpre} 1>${outpre}_hsq.log 2>&1
            if [ $? -ne 0 ];then gst_warn "GEML for $(basename $outpre) failed: Check ${outpre}_hsq.log";rm -f ${outpre}.hsq;fi
        fi
        if [ -s "${outpre}.hsq" ];then
            rm -f ${outpre}_hsq.out &&\
            mv ${outpre}.hsq ${outpre}_hsq.out
        fi
        if [ $? -ne 0 ];then gst_err "vcf2hsq for $outpre failed: Non-zero exit"; exit 1;fi
        rm -f ${outpre}.{bed,bim,fam,nosex,plinklog,} ${outpre}.grm*
    # fi
}
export -f vcf2hsq

get_each_variant_hsq () {
    #?  get hsq using GCTA for each variant in each peak, input peakXXX.vcf.gz
    local rawvcf=$1
    local peak=$(basename ${rawvcf%%.*})
    local hsq_all="NA";
    local hsq_SNP="NA";
    local hsq_SV="NA";
    local outpre=$out/99_Each_Variant_HSQ/${peak}
    mkdir -p $out/99_Each_Variant_HSQ
    check_files_exists $rawvcf
    check_sftw_path plink gcta64
    check_var_empty outpre
    check_files_exists ${out}/99_Misc/trait.txt
    local invcf=${outpre}_HSQ.vcf.gz
    bcftools view -S <(cut -f 2 ${out}/99_Misc/trait.txt) $rawvcf -O z -o $invcf
    if [ $? -ne 0 ];then gst_err "Get $invcf failed: Non-zero exit"; exit 1;fi
    # ? for each variant, compare it with all other snps and all other variants in the peak
    csvtk cut -tTH -f 3 $invcf > ${outpre}_ALL.id &&\
    grep -P "\.s_" ${outpre}_ALL.id >${outpre}_SNP.id
    if [ $? -ne 0 ];then gst_err "get ID failed: Non-zero exit"; exit 1;fi
    
    check_files_exists ${outpre}_ALL.id ${outpre}_SNP.id
    local nall=$( cat ${outpre}_ALL.id | wc -l )
    local nsnp=$( cat ${outpre}_SNP.id | wc -l )
    # local nsv=$( cat ${outpre}_SV.id | wc -l )
    # ? calc hsq for all
    # gst_log "ALL" #! test
    vcf2hsq $invcf ${out}/99_Misc/trait.txt ${outpre}_ALL ${outpre}_ALL.id ${out}/99_Misc/covar.txt
    if [ -s "${outpre}_ALL_hsq.out" ];then
        hsq_all=$(grep "V(G)/Vp" ${outpre}_ALL_hsq.out | cut -f 2)
    fi
    # ? calc hsq for SNP
    # gst_log "SNP" #! test
    if [ "$nsnp" -ge 1 ];then
        vcf2hsq $invcf ${out}/99_Misc/trait.txt ${outpre}_SNP ${outpre}_SNP.id ${out}/99_Misc/covar.txt
        if [ -s "${outpre}_SNP_hsq.out" ];then
            hsq_SNP=$(grep "V(G)/Vp" ${outpre}_SNP_hsq.out | cut -f 2)
            check_var_empty hsq_SNP
        fi
    fi
    local cur_variant=""
    local cur_type=""
    local hsq_cur=""
    local hsq_rest_snp=""
    local hsq_rest_all=""

    # ? dealing with each variant
    cat ${outpre}_ALL.id | while read cur_variant
    do
        # gst_log "$cur_variant" #! test
        hsq_cur="NA"
        hsq_rest_all="NA"
        hsq_rest_snp="NA"
        # ? calc hsq for each variant
        echo "$cur_variant" > ${outpre}_cur_variant.id

        # ? get Hcur
        vcf2hsq $invcf ${out}/99_Misc/trait.txt ${outpre}_cur_variant ${outpre}_cur_variant.id ${out}/99_Misc/covar.txt
        if [ -s "${outpre}_cur_variant_hsq.out" ];then
            hsq_cur=$(grep "V(G)/Vp" ${outpre}_cur_variant_hsq.out | cut -f 2)
        fi

        # ? get Hsnp and Hrest
        cur_type="";
        cur_type=$(echo "$cur_variant" | perl -lne '
            $type="SVINDEL";
            $type="SNP" if /\.s_/;
            print "$type";
        ' )
        check_var_empty cur_type
        # ? dealing with rest snp

        if [ "$cur_type" == "SNP" ];then
            grep -v -w "$cur_variant" ${outpre}_SNP.id > ${outpre}_rest_snp.id
            if [ ! -s "${outpre}_rest_snp.id" ];then
                gst_warn "No other SNPs for $cur_type"
            else
                rm -f ${outpre}_rest_snp_hsq.out
                vcf2hsq $invcf ${out}/99_Misc/trait.txt ${outpre}_rest_snp ${outpre}_rest_snp.id ${out}/99_Misc/covar.txt
                if [ -s "${outpre}_rest_snp_hsq.out" ];then
                    hsq_rest_snp=$(grep "V(G)/Vp" ${outpre}_rest_snp_hsq.out | cut -f 2)
                else
                    gst_warn "No result for vcf2hsq ${outpre}_rest_snp "
                fi
            fi
        else
            hsq_rest_snp="$hsq_SNP"
        fi
        # ? dealing with rest all
        
        grep -v -w "$cur_variant" ${outpre}_ALL.id > ${outpre}_rest_all.id
        if [ ! -s "${outpre}_rest_all.id" ];then
            gst_err "No other Variants for $cur_type"
            # exit 1
            hsq_rest_all="NA"
        else
            rm -f ${outpre}_rest_all_hsq.out
            vcf2hsq $invcf ${out}/99_Misc/trait.txt ${outpre}_rest_all ${outpre}_rest_all.id ${out}/99_Misc/covar.txt
            if [ -s "${outpre}_rest_all_hsq.out" ];then
                hsq_rest_all=$(grep "V(G)/Vp" ${outpre}_rest_all_hsq.out | cut -f 2)
                if [ -z "$hsq_rest_all" ]; then
                    gst_err "No hsq_rest_all for $cur_variant"
                    # exit 1
                fi
            fi
        fi
        # gst_log "${cur_variant}\t${hsq_cur}" #! test
        # ? get output for each variant
        echo -e "${trait}\t${peak}\t${cur_variant}\t${hsq_cur}\t${hsq_rest_all}\t${hsq_rest_snp}\t${hsq_all}\t${hsq_SNP}"
        rm -f ${outpre}_cur_variant* ${outpre}_rest_*
    done

    # rm -f ${outpre}_*.id ${outpre}_*_hsq.out ${outpre}_*.log
}
export -f get_each_variant_hsq

# >>>>>>>>>>>>>>>>>>>>>>>> main >>>>>>>>>>>>>>>>>>>>>>>>
gst_log ">>>> Get GWAS Peaks ..."
# ? Get GWAS Peaks 
get_credible_set_script
check_files_exists $out/99_Misc/get_credible_set.pl
run_harvester
check_files_exists $out/99_Misc/HARVESTER_peak.bed
peak_no=$(cat $out/99_Misc/HARVESTER_peak.bed | wc -l)
gst_log "Fine mapping within each peak ..."
if [ -s "$out/00_All_credible_sets.tsv" ];then
    gst_warn "Result exists, skip running."
else
    cat $out/99_Misc/HARVESTER_peak.bed | perl -F"\t" -lane '
        $region="$F[0]:$F[1]-$F[2]";$pid=$F[-1];print "$region\t$pid";
    ' | while read cur_region cur_peak
    do
        gst_log ">>>> Dealing with $cur_peak ($cur_region) of $peak_no peaks ..."
        {
        # ? get cur_peak gwas
        get_z_ld $cur_region $cur_peak
        # ? run dap-g
        run_dapg_zld $cur_peak
        } 1> $out/${cur_peak}_dapg.log 2>&1
    done
    if [ $? -ne 0 ];then gst_err "Dealing with each peak failed: Non-zero exit"; exit 1;fi
    gst_log "Merge all credible sets ..."
    cat $out/Peak*_dap-g.credible_set | perl -F"\t" -lane '
        $.==1 && print && next;
        print if $F[0] and $F[0] ne "Trait";
    ' > $out/00_All_credible_sets.tsv
    if [ $? -ne 0 ];then gst_err "get all credible sets failed: Non-zero exit";rm -f $out/00_All_credible_sets.tsv; exit 1;fi

    gst_log "Get merged all PIPs and tar each peak ..."
    # ? mv each peak to a subset dir
    mkdir -p $out/99_Each_Peak &&\
    rm -f $out/Peak*.tmp $out/Peak*.raw &&\
    mv $out/Peak* $out/99_Each_Peak
    # if [ $? -ne 0 ];then gst_err "mv each peak out into $out/99_Each_Peak failed: Non-zero exit"; exit 1;fi
fi

# ? get final outs ...
gst_log "Get final outs ..."
cur_wd=$PWD
ncs=$(cat $out/00_All_credible_sets.tsv | wc -l)
if [ "$ncs" -le 1 ];then
    gst_warn "No Credible Set found !"
else
    # ? gst_log "Get 00_Final_finemap_credible_sets.tsv ..."
    gst_log "Get 00_Final_finemap_credible_sets.tsv ..."
    if [ ! -s "$out/00_Final_finemap_credible_sets.tsv" ];then
        # ? get final peak info and merge all into one
        if [[ "$(ls $out/99_Each_Peak/Peak*_addz.tsv | wc -l)" -eq 1 ]];then
            cp $out/99_Each_Peak/Peak*_addz.tsv $out/merged_cs_gwas.tmp
            if [ $? -ne 0 ];then gst_err "cp to merged_cs_gwas failed: Non-zero exit"; exit 1;fi
        else
            csvtk concat -tTH -j $threads $out/99_Each_Peak/Peak*_addz.tsv -o $out/merged_cs_gwas.tmp
            if [ $? -ne 0 ];then gst_err "concat merged_cs_gwas failed: Non-zero exit"; exit 1;fi
        fi
        cat $out/99_Misc/HARVESTER_peak.bed | perl -F"\t" -lane '
            # chr start end max range count spacing GQS peak_id
            BEGIN{
                $,="\t";
                print "Peak_Region\tPeak_Max-logP\tPeak_Range\tPeak_Counts\tPeak_spacing\tPeak_GQS\tPeak_ID";
            }
            print "$F[0]:$F[1]-$F[2]",@F[3..$#F];
        ' > $out/merged_cs_peak.tmp &&\
        csvtk join -tTH -f "3;1" $out/00_All_credible_sets.tsv $out/merged_cs_gwas.tmp -j $threads | perl -F"\t" -lane '
            BEGIN{
                $,="\t";
                print "Trait\tSignal_Cluster\tVariant_ID\tPIP\tPeak_ID\tChr\tStart\tEnd\tP_value\tEffect_size\tZ_score";
            }
            $F[0]=~s/-Peak\d+$//;
            print @F;
        ' > $out/merged_cs_csgwas.tmp &&\
        csvtk join -tT -f "Peak_ID" $out/merged_cs_csgwas.tmp $out/merged_cs_peak.tmp -j $threads |\
        csvtk cut -tTH -f 1,6-8,3,9-11,4,2,5,12-17 -j $threads -o $out/00_Final_finemap_credible_sets.tsv
        if [ $? -ne 0 ];then gst_err "get final output failed: Non-zero exit";rm -f $out/00_Final_finemap_credible_sets.tsv; exit 1;fi
        rm -f $out/merged_cs_*.tmp
    fi
    # ? gst_log "Get anno and pip for coloc ..."
    gst_log "Get anno and pip for coloc ..."
    if [ -s "$out/${trait}_annote.vcf.gz" -a -s "$out/${trait}.pip.gz" ];then
        gst_warn "Already got anno and pip. Skip."
    else
        # ? get anno and pip for coloc
        gst_log "Get each anno and pip ..."
        dap_g_list=$( ls $out/99_Each_Peak/Peak*_dap-g.out )
        if [ -z "$dap_g_list" ]; then
            gst_warn "No dap-g out. Skip"
        else
            parallel -j $threads -q get_anno_pip :::: <(echo "$dap_g_list")
        fi
        if [ $? -ne 0 ];then gst_err "get each anno_pip failed: Non-zero exit"; exit 1;fi

        # ? merge all anno
        gst_log "Merge all anno ..."
        cat $out/*_anno.vcf | csvtk collapse -tTH -j $threads -f 1,2,3,4,5 -v 6 -s "|" | csvtk sort -tTH -k 1 -k 2:n -j $threads -o $out/${trait}_annote.vcf.gz
        if [ $? -ne 0 ];then gst_err "get merged anno failed: Non-zero exit";rm -f $out/${trait}_annote.vcf.gz; exit 1;fi
        rm -f $out/*_anno.vcf

        # ? merge all pip
        gst_log "Merge all pip ..."
        cat $out/*_PIPs.pip |gzip -c > $out/${trait}.pip.gz
        if [ $? -ne 0 ];then gst_err "get merged pip failed: Non-zero exit";rm -f $out/${trait}.pip.gz; exit 1;fi
        rm -f $out/*_PIPs.pip
    fi

if [ "$heritability" == "TRUE" ];then
    # ? get each Variant HSQ
    if [ -s "$out/99_Each_Peak.tar" ];then
        gst_log "Untar $out/99_Each_Peak.tar ..."
        cd $out &&\
        tar --no-overwrite-dir -xf 99_Each_Peak.tar &&\
        cd $cur_wd
        if [ $? -ne 0 ];then gst_err "Untar $out/99_Each_Peak.tar failed: Non-zero exit"; exit 1;fi
    fi
    # ? fmt trait and covar ...
    gst_log "Fmt trait and covar file ..."
    bcftools view -h $(ls ${out}/99_Each_Peak/*.vcf.gz | head -n 1) | tail -n 1 | cut -f 10- | csvtk transpose -tTH -o ${out}/99_Misc/vcf_samples.txt &&\
    csvtk join -tTH -f 1 ${out}/99_Misc/vcf_samples.txt $trait_file | perl -lane '$,="\t";print @F[0,0,1] if ++$h{$F[0]}==1;' > ${out}/99_Misc/trait.txt &&\
    csvtk join -tTH -f 1 <(cut -f 2 ${out}/99_Misc/trait.txt) $covar_file | perl -lane '$,="\t";print $F[0],@F if ++$h{$F[0]}==1;' > ${out}/99_Misc/covar.txt
    if [ $? -ne 0 ];then gst_err "fmt trait file failed: Non-zero exit"; exit 1;fi
    # ? check if covar matches with trait
    md5covar=$(cut -f 2 ${out}/99_Misc/covar.txt | sort | md5sum)
    md5trait=$(cut -f 2 ${out}/99_Misc/trait.txt | sort |md5sum)
    if [ "$md5covar" != "$md5trait" ];then
        gst_err "samples in $covar_file and $trait_file did not match !"
        exit 1
    fi
    gst_log "Estimate heritability for each variant..."
    if [ -s "$out/00_Each_variant_HSQ.tsv" ];then
        gst_warn "Skip running. Results exits."
    else
        # get trait that has same name with vcf
        vcf_list=$( ls $out/99_Each_Peak/Peak*_dap-g.out | sed 's/_dap-g.out/.vcf.gz/' )
        echo -e "Trait\tPeak_ID\tVID\tHSQ_Variant\tHSQ_rest_all\tHSQ_rest_snp\thsq_all\thsq_snp" > $out/00_Each_variant_HSQ.tsv
        if [ -z "$vcf_list" ]; then
            gst_warn "No avalible vcf. Skip"
            rm -f $out/00_Each_variant_HSQ.tsv
        else
            parallel -j $threads -q get_each_variant_hsq :::: <(echo "$vcf_list") | cat >> $out/00_Each_variant_HSQ.tsv
            if [ $? -ne 0 ];then gst_err "get each_SV_HSQ failed: Non-zero exit";rm -f $out/00_Each_variant_HSQ.tsv; exit 1;fi
        fi
    fi
    # ? tar $out/99_Each_Variant_HSQ
    if [ "$tar_value" == "TRUE" ];then
        gst_log "TAR $out/99_Each_Variant_HSQ ..."
        if [ ! -s "$out/99_Each_Variant_HSQ.tar" ];then
            cd $out &&\
            tar cvf 99_Each_Variant_HSQ.tar 99_Each_Variant_HSQ > 99_Each_Variant_HSQ_file_list.txt &&\
            cd $cur_wd
            if [ $? -ne 0 ];then gst_err "Tar 99_Each_Variant_HSQ failed: Non-zero exit";rm -f 99_Each_Variant_HSQ.tar; exit 1;fi
            # rm -rf $out/99_Each_Variant_HSQ
        fi
        if [ -s "$out/99_Each_Variant_HSQ.tar" ];then
            rm -rf $out/99_Each_Variant_HSQ
            # echo "Done"
        fi
    fi
fi
    # todo :>>> add get final results
    # ? joint gwas summary, pip info, hsq info into one final file
    # todo :<<<
    gst_log "Get final GWAS summary data"

if [ "$tar_value" == "TRUE" ];then
    gst_log "TAR dirs ..."
    # ? tar $out/99_Each_Peak
    gst_log "TAR $out/99_Each_Peak ..."

    if [ ! -s "$out/99_Each_Peak.tar" ];then
        cd $out &&\
        tar cvf 99_Each_Peak.tar 99_Each_Peak > 99_Each_Peak_file_list.txt &&\
        cd $cur_wd
        if [ $? -ne 0 ];then gst_err "Tar 99_Each_Peak failed: Non-zero exit";rm -f 99_Each_Peak.tar; exit 1;fi
    fi
    if [ -s "$out/99_Each_Peak.tar" ];then
        rm -rf $out/99_Each_Peak
        # echo "Done"
    fi
fi
    
gst_log "
---------------------------------------------------------
Merged Credible Set: $out/00_Final_finemap_credible_sets.tsv
Annotation vcf     : $out/${trait}_annote.vcf.gz
PIP file           : $out/${trait}.pip.gz
---------------------------------------------------------
"
fi
gst_log "
---------------------------------------------------------
|                      ALL DONE                         |
---------------------------------------------------------
"
# <<<<<<<<<<<<<<<<<<<<<<<< main <<<<<<<<<<<<<<<<<<<<<<<<


