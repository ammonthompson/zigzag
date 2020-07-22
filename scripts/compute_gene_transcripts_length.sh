#!/bin/bash

#type compute_gene_transcripts_length.sh -h (--help) for more info

### ALGORITHM:
### Extract exons, transcript id and gene id from columns 4, 5, and 9
### only exons that belong to transcripts which belong to genes will be extracted
### Be careful that all exon lines have a gene and a transcript id

gene_id="gene_id"
transcript_id="transcript_id"
threads=1

while [[ -n $1 ]];do
	case $1 in
		-g|--gene_id)
		shift
		gene_id=$1
		;;
		-t|--transcript_id)
		shift
		transcript_id=$1
		;;
		-f|--gtf_file)
		shift
		gtf_file=$1
		;;
		-c|--threads)
		shift
		threads=$1
		;;
		-h|--help)
		shift
		echo  ' '-g, --gene_id [tag identifier for gene id\'s in column 12] default: \"gene_id\" $'\n' \
		 -t, --transcript_id [tag identifier for transcript id\'s in column 12] default: \"transcript_id\"$'\n' \
		 -f, --gtf_file [annotation file]$'\n' \
		 -c, --threads [number of parallel threads] default: 1
		exit 1
		;;
		*)
		echo -e '\033[0;31m'Error: \"$1\" isnt a recognized parameter
                echo Use:
		echo  ' '-g, --gene_id [tag identifier for gene id\'s in column 12]$'\n' \
                 -t, --transcript_id [tag identifier for transcript id\'s in column 12]$'\n' \
                 -f, --gtf_file [annotation file]$'\n' \
                 -c, --threads [number of parallel threads]
		exit 1
		;;
	esac
	shift
done
	
if [[ ! -n $gtf_file ]];then
	echo -e '\033[0;31m'"You forgot to include the --gtf_file"
	exit 1
fi
	


filename=$(echo $gtf_file |rev |cut -d '/' -f 1 |rev)


### Extract exons, transcript id and gene id from columns 4, 5, and 9
### only exons that belong to transcripts which belong to genes will be extracted
### Be careful that all exon lines have a gene and a transcript id
grep $'\t'exon$'\t' $gtf_file |cut -f 4,5,9 | grep $gene_id |grep $transcript_id > reduced_${filename}.temp


# column 9 contains gene id and transcript id 
# and is assumed to be ; delimmited according to gtf,gff,gff3 format
if [[ ! -f gene_transcript_exonLengths_$filename ]];then

paste <(grep -Eo "${gene_id}[^;]+;?" reduced_${filename}.temp |sed -r "s/${gene_id}[ ,=,\"]+//g" |sed 's/;//g'|sed 's/\"//g') \
<(grep -Eo "${transcript_id}[^;]+;?" reduced_${filename}.temp |sed -r "s/${transcript_id}[ ,=,\"]+//g" |sed 's/;//g'|sed 's/\"//g') \
<(Rscript <(cut -f 1,2 reduced_${filename}.temp| sed 's/\t/-/g'|sed 's/^/abs(/g'|sed 's/$/)/g')|sed 's/\[1\] //g') > gene_transcript_exonLengths_$filename

fi

rm reduced_${filename}.temp

echo gene_id$'\t'transcript_id$'\t'transcript_length > transcript_lengths_$filename

# Function for computing transcript lengths from input transcript intermediate file with one line per transcript

get_transcript_lengths () {

	transcript_array=($(cat $1))
	for gene_transcript in ${transcript_array[@]};do

		transcript=$(echo $gene_transcript |cut -d ':' -f 2|sed -r 's/[ ]+//g')

		if [[ ${#transcript} -eq 0 ]];then
			transcript=$(echo $gene_transcript|cut -d ':' -f 1)
		fi
	
		echo $gene_transcript$'\t'$(echo $(grep -E $transcript[^A-Z,a-z,0-9\.]+ gene_transcript_exonLengths_$filename |cut -f 3) |sed 's/ /+/g' |bc )|sed 's/:/\t/g' >> ${1}_$filename
        
	done
}

# Make transcript length files

echo "Making transcript length file"
cut -f 1,2 gene_transcript_exonLengths_$filename | sed 's/\t/:/g' |sort -V |uniq  > transcript_list.text
split -d -n l/$threads transcript_list.text temp_tscript

for tt in temp_tscript*;do

	get_transcript_lengths $tt &

done
wait

cat temp_tscript*_$filename >> transcript_lengths_$filename
rm temp_tscript* transcript_list.text

# make mean gene length file

echo "Making average gene transcript length file"
echo gene_id$'\t'mean_length > gene_meanLengths_$filename
for gene in $(cut -f 1 transcript_lengths_$filename |sort -V|uniq);do

	mean_length=$(echo \( $(echo $(grep ${gene}$'\t' transcript_lengths_$filename | cut -f 3) |sed 's/ /+/g') \) / $(grep -c ${gene}$'\t' transcript_lengths_$filename) |bc -l | sed -r 's/([0-9]+\.[0-9]{3})[0-9]+$/\1/g')	
	echo $gene$'\t'$mean_length >> gene_meanLengths_$filename

done


rm gene_transcript_exonLengths_$filename

echo Finished
