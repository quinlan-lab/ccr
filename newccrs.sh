## folder made by date, in case we make major changes to exac-regions.py or resid-plot.py ##
date=2018_04_20 # default date value
if [ -s coveragecut.txt ]; then
    rm coveragecut.txt segdupscut.txt selfchaincut.txt #deletioncut.txt
fi
q="-x"
nochrom+="Y"
while getopts ":v:x:d:r:q:usegcnwf:" opt; do
    case $opt in
        v)
            echo "-version (date) was triggered, input: $OPTARG" >&2
            date=$OPTARG
            ;;
        x)
            echo "-exclude files passed into model" >&2
            x="-s"
            exclude+=("$OPTARG")
            echo "-exclude BED.gz input: '${exclude[@]}'" >&2
            ;;
        s)
            echo "-synonymous variant density input into the model was triggered" >&2
            syn="-s"
            s="-synonymous"
            ;;
        u)
            echo "-only calculating synonymous density, not adding to the model" >&2
            u="-r"
            ;;
        c)
            echo "-CpG density input into the model was triggered" >&2
            cpg="-c"
            c="-cpg"
            ;;
        n)
            echo "-no singleton input into the model was triggered" >&2
            ns="-n"
            n="-nosingletons"
            ;;
        f)
            echo "-fasta input into the model was specified, input: $OPTARG" >&2
            file=$OPTARG
            ;;
        w)
            echo "-variant flags incorporated into the model" >&2
            var="-w"
            w="-novariant"
            ;;
        d)
            echo "depth and depth cutoff passed into the model" >&2
            depth+=("$OPTARG")
            echo "-depth input: '${depth[@]}'" >&2
            ;;
        r)
            echo "chromosomes to include (all but Y recommended)" >&2
            r="-p"
            chrom+=("$OPTARG")
            echo "-chrom input: '${chrom[@]}'" >&2
            ;;
        q)
            echo "chromosomes to exclude (only X and Y recommended)" >&2
            q="-x"
            nochrom+=("$OPTARG")
            echo "-nochrom input: '${nochrom[@]}'" >&2
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            if [ $OPTARG == "q" ]; then
                nochrom=()
            else
                echo "Option -$OPTARG requires an argument." >&2
                exit 1
            fi
            ;;
    esac
done
echo $nochrom
mkdir -p results/$date/tmp
## generates regions and residuals files
python exac-regions.py $ns $var -c "data/gnomad.exomes.coverage.summary.tsv.bgz" -e data/Homo_sapiens.GRCh37.85.gtf.gz -x data/gnomad2.1-vep-vt.vcf.gz -d "${depth[0]}" -l "${depth[1]}" $x "${exclude[@]}" -f data/grch37.fa > results/$date/exac-regions$n$w.txt
python resid-plot.py $ns $syn $cpg $var $u $r "${chrom[@]}" $q "${nochrom[@]}" -f results/$date/exac-regions$n$w.txt > results/$date/resids$c$s$n$w.txt
cat <(head -1 results/$date/resids$c$s$n$w.txt) <(sed '1d' results/$date/resids$c$s$n$w.txt | sort -k12,12nr) > results/$date/tmp/residsort$c$s$n$w.txt
python weightpercentile.py results/$date/tmp/residsort$c$s$n$w.txt > results/$date/weightedresiduals$c$s$n$w.txt
