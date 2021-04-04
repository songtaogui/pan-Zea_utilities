#!/usr/bin/env bash

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

check_files_executable(){
    local num_related_file=1
    local related_file=""
    for related_file in  "$@"
    do
        if [[ ! -x "$related_file" ]]; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> NOT EXECUTABLE: $related_file $" >&2
            let num_related_file++
        fi
    done
    [ "$num_related_file" -ne 1 ] && exit 1
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

check_R_lib () {
    local num_R_lib=1
    local tp_R_lib=""
    Rscript --vanilla --slave -e '
            argv=as.character(commandArgs(TRUE));
            if (all(argv %in% rownames(installed.packages()))) {
                quit(save="no", status=0)
            } else {
                quit(save="no", status=1)
            }
        ' $*
    if [ $? -ne 0 ];then gst_err "One or more of \"$*\" not installed in R"; exit 1;fi
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


export -f check_R_lib gst_log gst_rcd gst_warn gst_err check_files_executable check_var_empty check_var_numeric check_sftw_path check_suffix check_files_exists check_abs_path
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<
usage=$(
cat <<EOF
------------------------------------------------------------
PANZ_Regional_enrich: a wrapper of R/regioneR function, with
some DIY options.
------------------------------------------------------------
Dependence:
    R/regioneR package (http://bioconductor.org/packages/release/bioc/html/regioneR.html)
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -o, --out <str>            [R]      Output prefix.
    -g, --genome <str>         [R]      Genome range in bed format. Used as the bondaries of region manipulating. Eg:
                                            chr1    0   1234567
                                            chr2    0   1345678

    -a, --query <str>          [R]      Query region file in bed format, will do permutation based on this file.

    -b, --feature <str>        [R]      Feature region file in bed format, will count overlaps of each permutation of query region with this file to evaluate the enrichment.

    -r, --ref <str>            [O]      Set the range of query region permutation. Three type of parameters are supported:
                1. [bed:your_bed_file.bed] --> Use a region file as the permutation boundary, so make sure the query region file is subset of the input region file.
                2. [flank:<int>] --> use a flanking of <int> bp length of query region (both side) as the permutation boundary, eg: "flank:10000" will flank left 10kb and right 10kb of query region.
                3. [time:<int>] --> use a flanking of <int> * <length of each query region> bp of query region as the permutation boundary, eg:
                    if your query region is:
                    #    chr1    1000    1100
                    #    chr2    5010    5020
                    and you set ref as "--ref time:5", the permutation boundary would be:
                    #    chr1    500     1600    (query length=100,flanking=100*5)
                    #    chr2    4960    5070    (query length=10,flanking=10*5)
                The default behavior of --ref (if unset) is to use "--genome" file as boundary, that is permutation along the whole genome.

    -n, --ntimes <int>         [O]      Number of permutation times.A large number of permutations will produce more accurate results and a nicer-looking plot but a permutation test can be computationally expensive.(default: 100)

    --cutoff <0-1>             [O]      P-value cutoff to call the result as significantly non-random. (default: 0.05)

    --seed <int>               [O]      Set random seeds to create reproducible results. (default: 1234)

    --plot                     [O]      Plot the permutation results.

    --force_save               [O]      Force save the permutation result rdata. The default behavior is to save only the result that passed the P-value cutoff.
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
export out=""
export genome=""
export query=""
export feature=""
export ref=""
export ntimes=100
export cutoff=0.05
export seed=1234
export plot="FALSE"
export force_save="FALSE"
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -o|--out)
            #echo "set argument \"$1\" with value: $2" >&2
            out=$2
            shift 2
        ;;
        -g|--genome)
            #echo "set argument \"$1\" with value: $2" >&2
            genome=$2
            shift 2
        ;;
        -a|--query)
            #echo "set argument \"$1\" with value: $2" >&2
            query=$2
            shift 2
        ;;
        -b|--feature)
            #echo "set argument \"$1\" with value: $2" >&2
            feature=$2
            shift 2
        ;;
        -r|--ref)
            #echo "set argument \"$1\" with value: $2" >&2
            ref=$2
            shift 2
        ;;
        -n|--ntimes)
            #echo "set argument \"$1\" with value: $2" >&2
            ntimes=$2
            shift 2
        ;;
        --cutoff)
            #echo "set argument \"$1\" with value: $2" >&2
            cutoff=$2
            shift 2
        ;;
        --seed)
            #echo "set argument \"$1\" with value: $2" >&2
            seed=$2
            shift 2
        ;;
        --plot)
            #echo "set argument \"$1\" with value: $2" >&2
            plot="TRUE"
            shift 1
        ;;
        --force_save)
            #echo "set argument \"$1\" with value: $2" >&2
            force_save="TRUE"
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
check_var_empty out genome query feature ntimes cutoff seed plot force_save
check_var_numeric ntimes cutoff seed
check_files_exists $query $genome $feature
check_R_lib regioneR
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# USAGE:$0
# USAGE:0 input1 input2
gst_regioneR () {
    local a=$query
    local b=$feature
    local r=${out}.ref_file.tmp
    local g=$genome
    local s=$seed
    local n=$ntimes
    local c=$cutoff
    local p=$plot
    local f=$force_save
    local o=$out
    check_files_exists $r $g $a $b
    Rscript -e '
        library(regioneR)
        argv=as.character(commandArgs(TRUE));
        print(argv)
        # ? 9 argvs
        G=toGRanges(argv[4])
        A=toGRanges(argv[1]);
        B=toGRanges(argv[2]);
        R=toGRanges(argv[3]);
        set.seed(as.numeric(argv[5]));
        Pt = overlapPermTest(A, B, ntimes=as.numeric(argv[6]), force.parallel=TRUE,genome=G);
        Pt.table=data.frame(A=gsub("^.*/(.*)\\.bed","\\1",argv[1]),B=gsub("^.*/(.*)\\.bed","\\1",argv[2]),P_value=Pt$numOverlaps$pval,N=Pt$numOverlaps$ntimes,Z=Pt$numOverlaps$zscore)
        # ? output results
        write.table(Pt.table,file=paste0(argv[10],".tsv"),quote=F,row.names=F,sep="\t")
        if(argv[9] == "TRUE"){
            # ? force save rdata
            save.image(Pt,Pt.table,file=paste0(argv[10],".rdata"))
        }else{
            if(Pt$numOverlaps$pval < as.numeric(argv[7])){
                # ? save rdata if pass cutoff
                save(Pt,Pt.table,file=paste0(argv[10],".rdata"))
            }
        }
        if(argv[8] == "TRUE"){
            # ? draw plot
            pdf(paste0(argv[10],"_PermTest.pdf"))
            plot(Pt)
            dev.off()
        }
    ' $a $b $r $g $s $n $c $p $f $o
    if [ $? -ne 0 ];then gst_err "Rscript run failed: Non-zero exit"; exit 1;fi
}
export -f gst_regioneR


if [ ! -s "${out}.done" ];then
    # ? parse ref mode and generat ref_file
    gst_log "Parsing options ..."

    export ref_file=$genome
    export flank_bp=""
    export flank_len_time=""
    if [ -z "$ref" ]; then
        gst_warn "No subset ref region provided, use whole genome as ref region ..."
        ref="bed:$genome"
    fi

    case "$ref" in
        bed:*)
            check_files_exists ${ref##bed:}
            cat ${ref##bed:} > ${out}.ref_file.tmp
        ;;
        flank:*)
            flank_bp=${ref##flank:}
            check_var_numeric flank_bp
            # ? generat ref_file
            cat $query | perl -F"\t" -lane '
                BEGIN{
                    $,="\t";
                    $inputfile="$ENV{genome}";
                    open(IN,"$inputfile") or die("Cannot open file: $inputfile");
                    while(<IN>){
                        chomp;
                        # $hash{$_}=1;
                        ($cc,$cs,$ce)=split(/\t/,$_);
                        $chrlen{$cc}=$ce;
                    }
                }
                ($c,$s,$e)=@F[0,1,2];
                $s=$s-$ENV{flank_bp};
                $e=$e+$ENV{flank_bp};
                $s=0 if $s < 0;
                $e=$chrlen{$c} if $chrlen{$c} && $e > $chrlen{$c};
                print $c,$s,$e
            ' > ${out}.ref_file.tmp
        ;;
        time:*)
            flank_len_time=${ref##time:}
            check_var_numeric flank_len_time
            # ? generat ref_file
            cat $query | perl -F"\t" -lane '
                BEGIN{
                    $,="\t";
                    $inputfile="$ENV{genome}";
                    open(IN,"$inputfile") or die("Cannot open file: $inputfile");
                    while(<IN>){
                        chomp;
                        # $hash{$_}=1;
                        ($cc,$cs,$ce)=split(/\t/,$_);
                        $chrlen{$cc}=$ce;
                    }
                }
                ($c,$s,$e)=@F[0,1,2];
                $len=($e-$s)*$ENV{flank_len_time};
                die("Start larger than End for $_ !") if $len < 0;
                $s=$s-$len;
                $e=$e+$len;
                $s=0 if $s < 0;
                $e=$chrlen{$c} if $chrlen{$c} && $e > $chrlen{$c};
                print $c,$s,$e
            ' > ${out}.ref_file.tmp
        ;;
        *)
            gst_err "Wrong --ref option: $ref"
            exit 1
        ;;
    esac

    if [ $? -ne 0 ];then gst_err "parse $ref failed: Non-zero exit"; exit 1;fi
    check_files_exists ${out}.ref_file.tmp

    # ? run main
    gst_log "Run permutation and get outputs ..."

    gst_regioneR 1> ${out}.rlog 2>&1

    if [ $? -ne 0 ];then gst_err "gst_regioneR failed: Non-zero exit"; exit 1;fi

    echo "done" > ${out}.done
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove $PWD/${out}.done to force re-run." >&2
fi

rm -f ${out}.ref_file.tmp ${out}.rlog

gst_log "All done."





