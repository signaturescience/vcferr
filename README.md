# skimsim

The `skimskim` module is a lightweight error simulation framework. The tool operates on an input VCF and can probabilistically simulate the following error models for biallelic SNPs:

- Heterozygous drop out: (0,1) or (1,0) to (0,0)
- Homozygous drop out: (1,1) to (0,1)
- Heterozygous drop in: (0,0) to (0,1)
- Homozygous drop in: (0,1) or (1,0) to (1,1)

## Installation

The `skimsim` tool is delivered as a Python module. To run the `skimsim.py` script you will need:

- Python 3.6.x
- `pysam`
- `random`
- `argparse`

## Usage

The examples below demonstrate basic usage with the `example.vcf.gz` in the `data/` directory of this repository.

The following is a basic example of usage that simulates 20% heterozygous dropout:

```
./skimsim.py --sample 'sample1' --input_vcf=data/example.vcf.gz --p_het_dropout=0.2
```

By default, `skimsim` will stream output of the VCF with errors simulated. However, if an argument is given for `"output_vcf"` then the VCF will be written to disk:

```
./skimsim.py --sample 'sample1' --input_vcf=data/example.vcf.gz --output_vcf=data/sample1sim.example.vcf.gz --p_het_dropout=0.2
```

Note that multiple kinds of error can be simulated simulatenously:

```
./skimsim.py --sample 'sample1' --input_vcf=data/example.vcf.gz --p_het_dropout=0.2 --p_hom_dropin=0.1 --p_het_dropin=0.2 --p_hom_dropout=0
```

