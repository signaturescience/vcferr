# Set these variables
# Name of sample in vcf you want to mutate
id="sample1"
# Name of VCF you want to mutate a sample in
vcf="example.vcf.gz"
# Proportion of heterozygous sites to mutate
p_het_dropout="0.25"

# example vcf data came from vcflib vcfrandom and some merging/concatenating/bespoke editing
# sample 2 is all hom/ref, sample 3 is all hom-alt
# sample 1: positions 1-2 are hom-ref, pos 3-6 are het, pos 7-9 are hom alt, on both chromosomes
bcftools view ${vcf}

# Extract the sample of interest from the vcf and write a new vcf with just that sample
bcftools view -s $id $vcf -Oz -o tmp1.vcf.gz && tabix -f tmp1.vcf.gz
bcftools view tmp1.vcf.gz

# Get only the heterozygous sites from that sample
bcftools view -i 'GT="het"' tmp1.vcf.gz -Oz -o tmp2.vcf.gz && tabix -f tmp2.vcf.gz
bcftools view tmp2.vcf.gz

# How many hets do you have?
n=$(bcftools index --nrecords tmp2.vcf.gz)
echo $n

# Calculate how many SNPs to mutate: round(p_het_dropout*n)
nmut=$(echo $n $p_het_dropout | awk '{printf "%.0f", $1*$2}')
echo $nmut

# Pull out nmut # regions chromosome and position randomly
shuf -n $nmut <(bcftools query -f '%CHROM\t%POS\n' tmp2.vcf.gz) | sort -nk1,2 > mutated.tsv
cat mutated.tsv

# change those to hom ref and write out a new vcf with just the mutated hets
bcftools view -R mutated.tsv tmp2.vcf.gz | sed 's@0/1@0/0@g' | bcftools sort -Oz -o tmp3.vcf.gz && tabix -f tmp3.vcf.gz
bcftools view tmp3.vcf.gz

# Get the NOT mutated hets
bedtools intersect -header -v -a tmp2.vcf.gz -b tmp3.vcf.gz -wa | bcftools view -Oz -o tmp4.vcf.gz && tabix -f tmp4.vcf.gz
bcftools view tmp4.vcf.gz

# Combine the hets again (some mutated)
bcftools concat -a tmp3.vcf.gz tmp4.vcf.gz | bcftools sort -Oz -o tmp5.vcf.gz && tabix tmp5.vcf.gz
bcftools view tmp5.vcf.gz

# Get the non-hets
bcftools view -i 'GT!="het"' tmp1.vcf.gz -Oz -o tmp6.vcf.gz && tabix -f tmp6.vcf.gz
bcftools view tmp6.vcf.gz

# Cat them together
bcftools concat -a tmp5.vcf.gz tmp6.vcf.gz | bcftools sort -Oz -o mutated.vcf.gz && tabix mutated.vcf.gz
bcftools view mutated.vcf.gz

# Get the samples that are NOT the samples you mutated
bcftools view -s ^${id} ${vcf} -Oz -o tmp7.vcf.gz && tabix -f tmp7.vcf.gz
bcftools view tmp7.vcf.gz

# Put the mutated sample back together with the other samples
bcftools merge mutated.vcf.gz tmp7.vcf.gz | bcftools annotate -x INFO/AN,INFO/AC | bcftools sort -Oz -o mutated-combined.vcf.gz && tabix -f mutated-combined.vcf.gz
bcftools view mutated-combined.vcf.gz

# Clean up
rm -f tmp*.vcf.gz*

# Check your work
bcftools view -H -R mutated.tsv ${vcf}
bcftools view -H -R mutated.tsv mutated-combined.vcf.gz