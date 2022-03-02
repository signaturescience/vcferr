
## import third party deps
import argparse

## import packages
from .err import vcferr

if __name__ == '__main__':
	## Argument parsing
#	parser = argparse.ArgumentParser(description = "Simulate error rates in genotype calls")
#
#	parser.add_argument("--sample", help="ID of sample in VCF file to be simulated", default = "", required=True)
#	parser.add_argument("--input_vcf", help="VCF file to simulate ex: example.vcf.gz", default = "", required=True, type=input_file_type)
#	parser.add_argument("--output_vcf", help="Output VCF file containing simulated genotypes ex: example.sim.vcf.gz", default = None)
#	parser.add_argument("--p_rarr", help="Probability of heterozygous dropout (0,1) to (0,0)", default = 0.1, type=drop_rate_type)
#	parser.add_argument("--p_aara", help="Probability of homozygous alt dropout (1,1) to (0,1)", default = 0, type=drop_rate_type)
#	parser.add_argument("--p_rrra", help="Probability of heterozygous dropin (0,0) to (0,1)", default = 0, type=drop_rate_type)
#	parser.add_argument("--p_raaa", help="Probability of homozygous alt dropin (0,1) to (1,1)", default = 0, type=drop_rate_type)
#	parser.add_argument("--p_aarr", help="Probability of double homozygous alt dropout (1,1) to (0,0)", default = 0, type=drop_rate_type)
#	parser.add_argument("--p_rraa", help="Probability of double homozygous alt dropin (0,0) to (1,1)", default = 0, type=drop_rate_type)

#	args = parser.parse_args()
#	sample = args.sample
#	input_vcf = args.input_vcf
#	output_vcf = args.output_vcf
#	p_rarr = args.p_rarr
#	p_aara = args.p_aara
#	p_rrra = args.p_rrra
#	p_raaa = args.p_raaa
#	p_aarr = args.p_aarr
#	p_rraa = args.p_rraa

	## check dropin/dropout rates are < 1
#	if (p_rarr + p_raaa) > 1:
#		parser.error("Heterozygous dropout + homozygous alt dropin cannot be greater than 1")

#	if (p_aara + p_aarr) > 1:
#		parser.error("Homozygous dropout + double homozygous alt dropout cannot be greater than 1")
		
#	if (p_rrra + p_rraa) > 1:
#		parser.error("Heterozygous dropin + double homozygous alt dropin cannot be greater than 1")
		
#	vcferr(sample=sample, input_vcf=input_vcf, output_vcf=output_vcf, p_rarr=p_rarr, p_aara=p_aara, p_rrra=p_rrra, p_raaa=p_raaa, p_aarr=p_aarr, p_rraa=p_rraa)
	vcferr()

