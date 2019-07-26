#!/bin/bash

# assumes gtf format with 

gene_id=$1
transcript_id=$2
gtf_file=$3
threads=$4

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

#transcript set to sort through
#transcript_array=($(cut -f 1,2 gene_transcript_exonLengths_$filename | sed 's/\t/:/g' |sort -V |uniq ))


echo "Making transcript gene length file"
echo gene_id$'\t'transcript_id$'\t'transcript_length > transcript_lengths_$filename
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

cut -f 1,2 gene_transcript_exonLengths_$filename | sed 's/\t/:/g' |sort -V |uniq  > transcript_list.text
split -d -n l/$threads transcript_list.text temp_tscript

for tt in temp_tscript*;do

	get_transcript_lengths $tt &

done
wait

cat temp_tscript*_$filename >> transcript_lengths_$filename
rm temp_tscript* transcript_list.text


echo "Making average gene transcript length file"
echo gene_id$'\t'mean_length > gene_meanLengths_$filename
for gene in $(cut -f 1 transcript_lengths_$filename |sort -V|uniq);do

	mean_length=$(echo \( $(echo $(grep ${gene}$'\t' transcript_lengths_$filename | cut -f 3) |sed 's/ /+/g') \) / $(grep -c ${gene}$'\t' transcript_lengths_$filename) |bc -l | sed -r 's/([0-9]+\.[0-9]{3})[0-9]+$/\1/g')	
	echo $gene$'\t'$mean_length >> gene_meanLengths_$filename

done


rm gene_transcript_exonLengths_$filename 
