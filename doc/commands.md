# min_N = 1 Samples
./gIMble blocks \
	-o hmel.chr18.min_1_sample \
	-s input/hmel.samples.csv \
	-b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile \
	-t 4 \
	-a 1 && \
./gIMble variants \
	-o hmel.chr18.min_1_sample \
	-s input/hmel.samples.csv \
	-b hmel.chr18.min_1_sample.blocks.h5 \
	-g input/hmel.chr18.genomefile \
	-t 4 \
	-v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks \
	-o hmel.chr18.min_1_sample \
	-s input/hmel.samples.csv \
	-c input/hmel.chr18.chrom_coordinates.txt \
	-b hmel.chr18.min_1_sample.blocks.h5 \
	-g input/hmel.chr18.genomefile && \
./gIMble modify variants \
	-v hmel.chr18.min_1_sample.variants.h5 \
	-o hmel.chr18.min_1_sample \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-o hmel.chr18.min_1_sample \
	-s input/hmel.samples.csv \
	-b hmel.chr18.min_1_sample.modified.blocks.h5 \
	-g input/hmel.chr18.new_coordinates.genomefile \
	-v hmel.chr18.min_1_sample.modified.variants.h5 \
	-w 500 \
	-l 100 \
	-t 4



# min_N = 10 Samples
./gIMble blocks -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4 -a 10 && \
./gIMble variants -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b hmel.chr18.min_10_sample.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18.min_10_sample.blocks.h5 -g input/hmel.chr18.genomefile && \
./gIMble modify variants -v hmel.chr18.min_10_sample.variants.h5 -o hmel.chr18.min_10_sample -m 4 -M 100 && \
./gIMble windows -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b hmel.chr18.min_10_sample.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.min_10_sample.variants.modified.h5 -w 500 -l 100

# min_N = 5 Samples (go ahead for 4th part)
./gIMble blocks -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4 -a 5 && \
./gIMble variants -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b hmel.chr18.min_5_sample.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18.min_5_sample.blocks.h5 -g input/hmel.chr18.genomefile && \
./gIMble modify variants -v hmel.chr18.min_5_sample.variants.h5 -o hmel.chr18.min_5_sample -m 4 -M 100 && \
./gIMble windows -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b hmel.chr18.min_5_sample.modified.blocks.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.min_5_sample.modified.variants.h5 -w 500 -l 100 
./gIMble gridsearch -s ../gIMble/input/hmel.samples.csv -g input/hmel.chr18.genomefile -l models/model.IM.M_D2A.MM_D2A.txt -A "chi" -v hmel.chr18.min_5_sample.variants.modified.h5 -k 2 -o hmel.chr18.min_5_sample. -t 4 --mu 1.9e-9 --block_size 64 --migration_MLE 3.8866e-7 --time_MLE 4e6 --derived_MLE 0.5 --theta_low 0.4 --theta_high 1.2