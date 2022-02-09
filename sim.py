import pysam,random

vcf_in=pysam.VariantFile("example.vcf.gz")
homref = [0,0]
homalt = [1,1]
het = [0,1]
l_homref_het = []
l_homref_het.append(homref)
l_homref_het.append(het)

tmp_sample='sample1'

for variant in vcf_in:
	gt=variant.samples[tmp_sample]['GT']
	gt=list(gt)
	# set all gts to homref (universal het dropout)
	#if gt == het: gt = homref
	# or randomly with probability 
	# (25,75) with homref first would be a 25% chance of het dropout
	# output of random choices is a list so use [0] index to get first element
	if gt == het: gt = random.choices(l_homref_het, weights=(25,75))[0]
	print(gt)

