import pysam,random
import click

def error_rate(variable):
	try:
		var_ret = float(variable)
	except:
		print("Drop in and drop out rates must be a number")
		raise ValueError
	if var_ret <= 1 and var_ret >= 0:
		return var_ret
	else:
		print("Drop in and drop out rates must be between 0 and 1")
		raise ValueError

@click.command()
@click.argument('input_vcf',
	metavar='<input_vcf>'
)

@click.option('-s', '--sample',
	help="ID of sample in VCF file to be simulated",
	required=True
)

@click.option('-o', '--output_vcf',
	help="Output VCF file containing simulated genotypes ex: example.sim.vcf.gz",
	default=None
)

@click.option('-p_rarr', '--p_rarr',
	help="Probability of heterozygous dropout (0,1) to (0,0)",
	default=0.1,
	type=error_rate
)

@click.option('-p_aara', '--p_aara',
        help="Probability of homozygous alt dropout (1,1) to (0,1)",
        default=0,
        type=error_rate
)

@click.option('-p_rrra', '--p_rrra',
        help="Probability of heterozygous dropin (0,0) to (0,1)",
        default=0,
        type=error_rate
)

@click.option('-p_raaa', '--p_raaa',
        help="Probability of homozygous alt dropin (0,1) to (1,1)",
        default=0,
        type=error_rate
)

@click.option('-p_aarr', '--p_aarr',
        help="Probability of double homozygous alt dropout (1,1) to (0,0)",
        default=0,
        type=error_rate
)

@click.option('-p_rraa', '--p_rraa',
        help="Probability of double homozygous alt dropin (0,0) to (1,1)",
        default=0,
        type=error_rate
)
@click.pass_context
def vcferr(context,input_vcf,sample,output_vcf,p_rarr,p_aara,p_rrra,p_raaa,p_aarr,p_rraa):
	try:
		vcf_in = pysam.VariantFile(input_vcf)
	except:
		print("Error reading vcf input file "+input_vcf)
		raise ValueError
	## Ensure sampleID is in header of vcf file
	if not sample in list((vcf_in.header.samples)):
		print("VCF file does not appear to contain "+sample)
		raise ValueError
	## switch to allow streaming if no vcf_out is specified
	if output_vcf is None:
		vcf_out = pysam.VariantFile("-", 'w', header=vcf_in.header)
	else:
		vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

	recs = vcf_in.fetch()

	homref = [0,0]
	homalt = [1,1]
	## handle het as either 0,1 or 1,0
	het = [0,1]
	het2 = [1,0]

	homref_homalt_het = []
	homref_homalt_het.append(homref)
	homref_homalt_het.append(homalt)
	homref_homalt_het.append(het)
	
	## get weight for probability of het dropout
	phet_do_w = int(round(p_rarr,2)*100)
	## get weight for probability of hom dropout
	phom_do_w = int(round(p_aara,2)*100)
	## get weight for probability of het dropin
	phet_di_w = int(round(p_rrra,2)*100)
        ## get weight for probability of hom dropin
	phom_di_w = int(round(p_raaa,2)*100)
	## get weight for probability of double hom dropout
	phom_do2_w = int(round(p_aarr,2)*100)
	## get weight for probability of double hom alt dropin
	phom_di2_w = int(round(p_rraa,2)*100)

	
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
		## or double homozygous dropout (1,1) to (0,0)
		elif gt == homalt:
			gt = random.choices(homref_homalt_het, weights=(phom_do2_w, 100-(phom_do2_w+phom_do_w),phom_do_w))[0]
		## het dropin (0,0) to (0,1)
		## or double heterozygous dropin (0,0) to (1,1)
		elif gt == homref:
		  gt = random.choices(homref_homalt_het, weights=(100-(phom_di2_w+phet_di_w), phom_di2_w, phet_di_w))[0]
		rec.samples[sample]['GT'] = tuple(gt)
		vcf_out.write(rec)

