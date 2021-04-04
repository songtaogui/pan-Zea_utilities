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
export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<

usage="
------------------------------------------------------------
PANZ annotate SV with TE and Gene gff
------------------------------------------------------------
Dependence in PATH: parallel, RepeatMasker, bedtools, bcftools
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIIONS]
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
"
if [[ $# -ne 7 ]]; then 
    echo "$usage" >&2
    exit 1
fi

vcf=$1
gene_gff=$2
te_bed=$3
te_fa=$4
export flank=$5
threads=$6
export output=$7
overlap_TE=0.5

check_sftw_path parallel RepeatMasker bedtools bcftools
check_files_exists $vcf $gene_gff $te_bed $te_fa
check_var_numeric flank threads
check_var_empty output

# ? dealing with TEs
# ? Generate Non-Ins.bed and Ins.fa
gst_log "Dealing with TEs ..."
if [ -s "${output}/01_TE_Anno.tsv" ];then
    gst_warn "Skip running, result file exists."
else
    gst_log "Generate Non-Ins.bed and Ins.fa ..."
    if [ -s "${output}/Non_INS.bed" -a -s "${output}/INS.fa" ];then
        gst_warn "Skip running, result file exists."
    else
        bcftools view --threads $threads $vcf | parallel --pipe -q -k -j $threads perl -F"\t" -lane '
        BEGIN{$,="\t";}
        /^#/ && next;
        if(/SVTYPE=INS/){
            $type="INS";
            $out=join("\t",@F[2,4]);
            # ! RepeatMasker limits the Fasta max id length to 50
            print "$type,$out" if ++$h{$F[2]}==1 and length($F[2]) < 50;
        }else{
            $type="Non_INS";
            ($id,$c1,$p1,$c2,$p2)=split("#",$F[2]);
            if($c1 eq $c2){
                print "$type,$c1",$p1-1,$p2,$F[2] if ++$h{$F[2]}==1;
            }else{
                print "$type,$c1",$p1-1,$p1,$F[2] if ++$h{$F[2]}==1;
                print "$type,$c2",$p2-1,$p2,$F[2] if ++$h{$F[2]}==1;
            }
        }
        ' | csvtk split -H -f 1 -o ${output}
        if [ $? -ne 0 ];then gst_err "get ins/nonins bed fa failed: Non-zero exit";rm -rf ${output}; exit 1;fi
        # ? format bed and fa --> ${output}/Non_INS_Anno.tsv
        check_files_exists ${output}/stdin-INS.csv ${output}/stdin-Non_INS.csv
        csvtk cut -j $threads -f 2 ${output}/stdin-INS.csv | seqkit tab2fx -j $threads -o ${output}/INS.fa &&\
        csvtk cut -j $threads -f 2 ${output}/stdin-Non_INS.csv -o ${output}/Non_INS.bed &&\
        rm -f ${output}/stdin-INS.csv ${output}/stdin-Non_INS.csv
        if [ $? -ne 0 ];then gst_err "format bed and fa failed: Non-zero exit"; exit 1;fi
        check_files_exists ${output}/Non_INS.bed ${output}/INS.fa
        seqkit fx2tab -n -j $threads -i ${output}/INS.fa | parallel --pipe -q -k -j $threads perl -F"#" -lane '
            $,="\t";
            $F[-1]=~s/\s+$//;
            print $F[1],$F[2]-1,$F[4],join("#",@F);
        ' > ${output}/INS.bed
    fi

    gst_log "Check Non_INS SV that >= $overlap_TE overlapped with TEs ..."
    if [ ! -s "${output}/Non_INS_Anno.tsv" ];then
        bedtools intersect -e -f $overlap_TE -F $overlap_TE -a ${output}/Non_INS.bed -b $te_bed -wo |\
            parallel --pipe -j $threads -q -k perl -lane 'print "$F[3]\t$F[7]#OVLP_$F[10]" if $F[10]' |\
            csvtk collapse -j $threads -tTH -f 1 -v 2 -s "," -o ${output}/Non_INS_Anno.tsv
        if [ $? -ne 0 ];then gst_err "Annotate Non_INS failed: Non-zero exit";rm -f ${output}/Non_INS_Anno.tsv; exit 1;fi
    else
        gst_warn "Skip running, result file exists."
    fi

    # ? RM INS.fa with cutoff 250 --> ${output}/INS_Anno.tsv
    gst_log "Identify TE INS using RepeatMasker against $te_fa ..."
    if [ ! -s "${output}/INS_Anno.tsv" ];then
        RepeatMasker -cutoff 250 -nolow -no_is -norna -parallel $threads -lib $te_fa -dir ${output}/INS_RM ${output}/INS.fa
        if [ $? -ne 0 ];then gst_err "RepeatMasker for INS.fa failed: Non-zero exit";rm -rf ${output}/INS_RM; exit 1;fi
        # ? Parsing RM out
        check_files_exists ${output}/INS_RM/INS.fa.out
        sed '1,3d' ${output}/INS_RM/INS.fa.out | parallel --pipe -j $threads -q -k perl -lane 'print "$F[4]\t$F[9]#SW_$F[0]" if ++$h{$F[4]}==1' > ${output}/INS_Anno.tsv
        if [ $? -ne 0 ];then gst_err "Parse RM out failed: Non-zero exit"; rm -f ${output}/INS_Anno.tsv; exit 1;fi
    else
        gst_warn "Skip running, result file exists."
    fi

    gst_log "merge TE annotation"
    cat ${output}/Non_INS_Anno.tsv ${output}/INS_Anno.tsv | csvtk sort -j $threads -tTH -k 1 -o ${output}/01_TE_Anno.tsv
fi

gst_log "Done for TE annotation!
    "

gst_log "Dealing with Genes ..."
# ? get All SV Bed
if [ ! -s "${output}/ALL_SV_flanking.bed" ];then
    check_files_exists ${output}/INS.bed ${output}/Non_INS.bed
    cat ${output}/INS.bed ${output}/Non_INS.bed | parallel --pipe -j $threads -q -k perl -F"\t" -lane '
        $fa=$F[1]-$ENV{flank}; $fa=0 if $fa < 0;
        $fb=$F[2]+$ENV{flank};
        $,="\t";
        print $F[0],$fa,$fb,join("|",@F);
    ' | csvtk sort -tTH -k 1 -k 2:n -k 3:n -j $threads -o ${output}/ALL_SV_flanking.bed
    if [ $? -ne 0 ];then gst_err "get All flank bed failed: Non-zero exit";rm -f ${output}/ALL_SV_flanking.bed; exit 1;fi
fi

# ? filter gene gff to include only col3 =~ /RNA|CDS|UTR/ --> /RNA/
gst_log "Formatting annotation gff ..."

if [ ! -s "${output}/fmtted_gene.bed" ];then
    cat $gene_gff | parallel --pipe -q -k -j $threads perl -F"\t" -lane '
        next unless $F[2]=~/RNA/;
        $id="NA";
        if($F[-1]=~/ID=(.*?);/){
            $id=$1;
        }elsif($F[-1]=~/Parent=([^\s;]+)/){
            $id="$1"."_$F[2]";
        }else{
            next;
        }
        $,="\t";
        #$num=++$h{"$id,$F[6]"};
	#print @F[0,3,4],"$id\#$num",@F[2,6];
	print @F[0,3,4],$id,@F[2,6];
    ' > ${output}/fmtted_gene.bed
    if [ $? -ne 0 ];then gst_err "gt fmtted_gene.bed failed: Non-zero exit"; rm -f ${output}/fmtted_gene.bed; exit 1;fi
fi

# ? overlap with genes
# USAGE:0 <pipe with bedtools intersect output>
# 0        1        2               3             4          5      6            7                         8      9        10
# 1       37099   326679  PZ00001aSV00000102DEL   1       44351   44947   CDS:Zm00001d027230_P001_1       CDS     +       596
fmt_gene_ovlp () {
    # format gene ovlp result from intersectBed to final results
    local line=""
    while read line
    do
        echo "$line" | perl -F"\t" -lane '
            $,="\t";
            ($c,$s,$e,$id)=split(/\|/,$F[3]);
            print $c,$s,$e,$id,@F[4..$#F-1];
        ' |\
        bedtools overlap -i stdin -cols 2,3,6,7 |\
        perl -F"\t" -lane '
            $,="\t";
            $ovlp=abs($F[-1]);
            $tag="";
            if($F[-1] < 0){
                # ? no ovlp
                if($F[2] < $F[5]){
                    # ? upstream
                    $tag="Upstream";
                }else{
                    # ? downstream
                    $tag="Downstream";
                }
            }else{
                # ? ovlp
                $tag="Overlap";
                $tag="embed" if $F[1] >= $F[5] and $F[2] <= $F[6];
                $tag="contain" if $F[1] < $F[5] and $F[2] > $F[6];
            }
            print $F[3],"$F[7]#$tag";
        '
    done
}
export -f fmt_gene_ovlp

gst_log "Get overlapped gene features with flanking $flank bp ..."
if [ ! -s "${output}/02_gene_anno.tsv" ];then
    bedtools intersect -a ${output}/ALL_SV_flanking.bed -b ${output}/fmtted_gene.bed -wo |\
    parallel --pipe -j $threads -q -k fmt_gene_ovlp |\
    csvtk collapse -j $threads -tTH -f 1 -v 2 -s "," -o ${output}/02_gene_anno.tsv
else
    gst_warn "Skip running, result file exists."
fi

gst_log "Done Gene annotate!
"

gst_log "Adding annotations to VCF file ..."
if [ ! -s "${output}/00_Annotated_variants.vcf.gz" ];then
    gst_log "Formatting final vcf ..."
    bcftools view --threads $threads $vcf | perl -F"\t" -lane '
        BEGIN{
            $tefile="$ENV{output}/01_TE_Anno.tsv";
            open(TE,"$tefile") or die("Cannot open file: $tefile");
            while(<TE>){
                chomp;
                # $hash{$_}=1;
                ($var1,$var2)=split(/\t/,$_);
                $TE_hash{$var1}=$var2;
            }
            close TE;
            $genefile="$ENV{output}/02_gene_anno.tsv";
            open(GENE,"$genefile") or die("Cannot open file: $genefile");
            while(<GENE>){
                chomp;
                # $hash{$_}=1;
                ($var1,$var2)=split(/\t/,$_);
                $GENE_hash{$var1}=$var2;
            }
            close GENE;
            $,="\t";
        }
        /^##/ && print && next;
        if(/^#CHROM/){
            print "##INFO=<ID=TE_ANNO,Number=1,Type=String,Description=\"Annotation of overlapped TEs\">\n##INFO=<ID=GENE_ANNO,Number=1,Type=String,Description=\"Annotation of Genes overlapped with flanking $ENV{flank} bp\">";
            print $_;
            next;
        }
        if($TE_hash{$F[2]}){
            $F[7]=$F[7].";TE_ANNO=$TE_hash{$F[2]}";
        }
        if($GENE_hash{$F[2]}){
            $F[7]=$F[7].";GENE_ANNO=$GENE_hash{$F[2]}";
        }
        print @F if ++$h{$F[2]}==1;
    ' | bcftools view --threads $threads -O z -o ${output}/00_Annotated_variants.vcf.gz
    if [ $? -ne 0 ];then gst_err "get annotated vcf failed: Non-zero exit";rm -f ${output}/00_Annotated_variants.vcf.gz; exit 1;fi
fi
if [ ! -s "${output}/00_Annotated_variants.vcf.gz.tbi" ];then
    gst_log "Indexing ..."
    bcftools index --threads $threads -t ${output}/00_Annotated_variants.vcf.gz
    if [ $? -ne 0 ];then gst_err "index failed: Non-zero exit"; exit 1;fi
fi
gst_log "Yes!! All Done !!"
