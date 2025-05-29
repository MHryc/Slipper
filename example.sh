grep -v "@" "$file".sam | awk '{if($2==83) print $0}' > ${dir_results}/"$file"_163.sam # R2 +
grep -v "@" "$file".sam | awk '{if($2==147) print $0}' > ${dir_results}/"$file"_99.sam # R2 -
grep -v "@" "$file".sam | awk '{if($2==73) print $0}' > ${dir_results}/"$file"_73.sam # R1 +
grep -v "@" "$file".sam | awk '{if($2==89) print $0}' > ${dir_results}/"$file"_89.sam # R1 -

#Soft-clips with get_softclipped_reads_from_sam.pl script (from https://github.com/smaegol/LINE_1_RACE_seq_analysis)

perl get_softclipped_reads_from_sam.pl \
-input ${dir_results}/"$file"_163.sam  \
-output  ${dir_results}/"$file"_163.softclips

perl get_softclipped_reads_from_sam.pl \
-input ${dir_results}/"$file"_99.sam  \
-output  ${dir_results}/"$file"_99.softclips

perl get_softclipped_reads_from_sam.pl \
-input ${dir_results}/"$file"_73.sam  \
-output  ${dir_results}/"$file"_73.softclips

perl get_softclipped_reads_from_sam.pl \
-input ${dir_results}/"$file"_89.sam  \
-output  ${dir_results}/"$file"_89.softclips

#Parsing results
echo "Parsing results"

grep ">" ${dir_results}/${file}_163.softclips | sed 's/>//g ; s/clip5: //g; s/clip3: //g; s/ref: //g; s/pos: //g' > ${dir_results}/${file}_163.softclips.parsed
grep ">" ${dir_results}/${file}_99.softclips | sed 's/>//g ; s/clip5: //g; s/clip3: //g; s/ref: //g; s/pos: //g' > ${dir_results}/${file}_99.softclips.parsed
grep ">" ${dir_results}/${file}_73.softclips | sed 's/>//g ; s/clip5: //g; s/clip3: //g; s/ref: //g; s/pos: //g' > ${dir_results}/${file}_73.softclips.parsed
grep ">" ${dir_results}/${file}_89.softclips | sed 's/>//g ; s/clip5: //g; s/clip3: //g; s/ref: //g; s/pos: //g' > ${dir_results}/${file}_89.softclips.parsed

#Reformatting soft-clip output
echo "Reformatting soft-clip output"

cut -f 1,3,4,5 ${dir_results}/${file}_163.softclips.parsed > ${dir_results}/${file}_163.softclips.parsed.formatted #clip3
cut -f 1,3,4,5 ${dir_results}/${file}_99.softclips.parsed > ${dir_results}/${file}_99.softclips.parsed.formatted #clip3
cut -f 1,4,5 ${dir_results}/${file}_73.softclips.parsed > ${dir_results}/${file}_73.softclips.parsed.formatted
cut -f 1,4,5 ${dir_results}/${file}_89.softclips.parsed > ${dir_results}/${file}_89.softclips.parsed.formatted

#Adding a column with strand info
echo "Adding a column with strand info"

awk '{print $0, "\t163"}' ${dir_results}/${file}_163.softclips.parsed.formatted > ${dir_results}/${file}_163.softclips.parsed.strand
awk '{print $0, "\t99"}' ${dir_results}/${file}_99.softclips.parsed.formatted > ${dir_results}/${file}_99.softclips.parsed.strand

#Selecting R2 from sam based on mapped R1
echo "Selecting R2 from sam based on mapped R1"

code="73"
cut -f 1 ${dir_results}/${file}_${code}.softclips.parsed.formatted > ${dir_results}/${file}_${code}.names.txt
seqtk subseq "$file"_R2.extracted.reformat.fastq ${dir_results}/${file}_${code}.names.txt \
 >  ${dir_results}/"$file"_R2.extracted.names.fastq.txt
grep "@" ${dir_results}/"$file"_R2.extracted.names.fastq.txt | sed 's/@//g'> ${dir_results}/${file}_${code}.names2.txt
awk '/@/{getline; print}' ${dir_results}/"$file"_R2.extracted.names.fastq.txt | tr ACGTacgt TGCAtgca | rev > ${dir_results}/${file}_${code}.seq.txt
paste ${dir_results}/${file}_${code}.names2.txt ${dir_results}/${file}_${code}.seq.txt > ${dir_results}/${file}_${code}.names.seq.txt

Rscript seq.R \
-f ${dir_results}/${file}_${code}.softclips.parsed.formatted \
-n ${dir_results}/${file}_${code}.names.seq.txt \
-c ${code} \
-o ${dir_results}/${file}_${code}.softclips.parsed.strand

code="89"

cut -f 1 ${dir_results}/${file}_${code}.softclips.parsed.formatted > ${dir_results}/${file}_${code}.names.txt
seqtk subseq "$file"_R2.extracted.reformat.fastq ${dir_results}/${file}_${code}.names.txt \
 >  ${dir_results}/"$file"_R2.extracted.names.fastq.txt
grep "@" ${dir_results}/"$file"_R2.extracted.names.fastq.txt | sed 's/@//g'> ${dir_results}/${file}_${code}.names2.txt
awk '/@/{getline; print}' ${dir_results}/"$file"_R2.extracted.names.fastq.txt  | tr ACGTacgt TGCAtgca | rev > ${dir_results}/${file}_${code}.seq.txt
paste ${dir_results}/${file}_${code}.names2.txt ${dir_results}/${file}_${code}.seq.txt > ${dir_results}/${file}_${code}.names.seq.txt

Rscript seq.R \
-f ${dir_results}/${file}_${code}.softclips.parsed.formatted \
-n ${dir_results}/${file}_${code}.names.seq.txt \
-c ${code} \
-o ${dir_results}/${file}_${code}.softclips.parsed.strand

#Merging all files
echo "Merging files"
cat ${dir_results}/"$file"*.softclips.parsed.strand > ${dir_results}/${file}.softclips.parsed.all.csv

echo "Deleting intermediate files"
#rm ${dir_results}/${file}*.parsed.strand
#rm ${dir_results}/${file}*.parsed.formatted
#rm ${dir_results}/${file}*.parsed
#rm ${dir_results}/${file}*.txt
#rm ${dir_results}/"$file"_163.sam
#rm ${dir_results}/"$file"_99.sam
#rm ${dir_results}/"$file"_73.sam
#rm ${dir_results}/"$file"_89.sam
#rm *positions*

Rscript poli_clip_all.R \
-f ${dir_results}/${file}.softclips.parsed.all.csv \
-o ${dir_results}/${file}.softclips.results.csv