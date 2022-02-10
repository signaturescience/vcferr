#!/usr/bin/env python

import pysam,random,argparse

def skimsim(sample, input_vcf, p_het_dropout,p_hom_dropout,p_het_dropin,p_hom_dropin):
	vcf_in = pysam.VariantFile(input_vcf)
	vcf_out = pysam.VariantFile("-", 'w', header=vcf_in.header)
	recs = vcf_in.fetch()

	homref = [0,0]
	homalt = [1,1]
	## handle het as either 0,1 or 1,0
	het = [0,1]
	het2 = [1,0]

	homref_het = []
	homref_het.append(homref)
	homref_het.append(het)

	homalt_het = []
	homalt_het.append(homalt)
	homalt_het.append(het)

	homref_homalt = []
	homref_homalt.append(homref)
	homref_homalt.append(homalt)
	
	## get weight for probability of het dropout
	phet_do_w = int(round(p_het_dropout,2)*100)
	## get weight for probability of hom dropout
	phom_do_w = int(round(p_hom_dropout,2)*100)
	## get weight for probability of het dropin
	phet_di_w = int(round(p_het_dropin,2)*100)
        ## get weight for probability of hom dropin
	phom_di_w = int(round(p_hom_dropin,2)*100)
	
	for rec in recs:
		rec.samples[sample]['GT'] = list(rec.samples[sample]['GT'])
		gt=list(rec.samples[sample]['GT'])
		## use weights for random selections by prob
		## for example (25,75) with homref first would be a 25% chance of het dropout
		## output of random choices is a list so use [0] index to get first element
		## het dropout (0,1) to (0,0) or hom dropin (0,1) to (1,1)
		## NOTE: this requires a list with three choices
		if gt == het or gt == het2:
			## alternate way of building list here
			## need this to retain whatever coding (0,1) or (1,0) for gt
			homref_homalt_het=[]
			homref_homalt_het.append(homref)
			homref_homalt_het.append(homalt)
			homref_homalt_het.append(gt) 
			gt = random.choices(homref_homalt_het, weights=(phet_do_w,phom_di_w,100-(phet_do_w+phom_di_w)))[0]
		## hom dropout (1,1) to (0,1)
		elif gt == homalt:
			gt = random.choices(homalt_het, weights=(100-phom_do_w,phom_do_w))[0]
		## het dropin (0,0) to (0,1)
		elif gt == homref:
			gt = random.choices(homref_het, weights=(100-phet_di_w,phet_di_w))[0]
		rec.samples[sample]['GT'] = tuple(gt)
		vcf_out.write(rec)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Simulate error rates in genotype calls")

	parser.add_argument("--sample", help="ID of sample in VCF file to be simulated", default = "")
	parser.add_argument("--input_vcf", help="VCF file to simulate ex: example.vcf.gz", default = "")
	parser.add_argument("--output_vcf", help="Output VCF file containing simulated genotypes ex: example.sim.vcf.gz", default = "")
	parser.add_argument("--p_het_dropout", help="Probability of heterozygous dropout (0,1) to (0,0)", default = 0.1, type=float)
	parser.add_argument("--p_hom_dropout", help="Probability of homozygous dropout (1,1) to (0,1)", default = 0, type=float)
	parser.add_argument("--p_het_dropin", help="Probability of heterozygous dropin (0,0) to (0,1)", default = 0, type=float)
	parser.add_argument("--p_hom_dropin", help="Probability of homozygous dropin (0,1) to (1,1)", default = 0, type=float)

	args = parser.parse_args()
	sample = args.sample
	input_vcf = args.input_vcf
	p_het_dropout = args.p_het_dropout
	p_hom_dropout = args.p_hom_dropout
	p_het_dropin = args.p_het_dropin
	p_hom_dropin = args.p_hom_dropin

	skimsim(sample=sample, input_vcf=input_vcf, p_het_dropout=p_het_dropout, p_hom_dropout=p_hom_dropout, p_het_dropin=p_het_dropin, p_hom_dropin=p_hom_dropin)

