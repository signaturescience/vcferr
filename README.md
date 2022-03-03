# vcferr

The `vcferr` module is a lightweight error simulation framework. The tool operates on an input VCF and can probabilistically simulate the following error models for biallelic SNPs:

- rarr = Heterozygous drop out: (0,1) or (1,0) to (0,0)
- aara = Homozygous alt drop out: (1,1) to (0,1)
- rrra = Heterozygous drop in: (0,0) to (0,1)
- raaa = Homozygous alt drop in: (0,1) or (1,0) to (1,1)
- aarr = Double homozygous alt dropout: (1,1) to (0,0)
- rraa = Double homozygous alt dropin: (0,0) to (1,1)

## Installation

The `vcferr` tool is delivered as a Python module. To install clone this repo and use `pip` from the root of the directory:


```
pip install .
```

Note that the following dependencies are used by `vcferr`:

- Python >=3.6.x
- `pysam`
- `random`
- `click`

## Usage

The examples below demonstrate basic usage with the `example.vcf.gz` in the `data/` directory of this repository.

The following is a basic example of usage that simulates 20% heterozygous dropout:

```
vcferr data/example.vcf.gz --sample='sample1' --p_rarr=0.2
```

By default, `vcferr` will stream output of the VCF with errors simulated. However, if an argument is given for `"output_vcf"` then the VCF will be written to disk:

```
vcferr data/example.vcf.gz --sample='sample1' --output_vcf=data/sample1sim.example.vcf.gz --p_rarr=0.2
```

Note that multiple kinds of error can be simulated simulatenously:

```
vcferr data/example.vcf.gz --sample='sample1' --p_rarr=0.2 --p_raaa=0.1 --p_rrra=0.2 --p_aara=0
```

