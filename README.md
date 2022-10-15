# rlpbwt-test

Pipeline to test PBWT ([richarddurbin/pbwt](https://github.com/richarddurbin/pbwt)) and RLPBWT ([dlcgold/rlpbwt](https://github.com/dlcgold/rlpbwt)) performances.

## Usage

First time only:

```shell
curl -L https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh > miniconda.sh
bash miniconda.sh -b -f -p ./conda/
cd rlpbwt-test
../conda/bin/conda install -n base -c conda-forge mamba
source ../conda/bin/activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
../conda/bin/conda install -n snakemake -c conda-forge libcurl openssl csvkit
source ../conda/bin/deactivate
```

Run tests:

```shell
source ../conda/bin/activate snakemake
snakemake run_after_checkpoint --cores <number_of_cores>
source ../conda/bin/deactivate
```

Tests Edit:

In order to edit input panel, please provide a folder names "samples_x-sites_y" (where x/y are the quantity of sample/sites) containing a `MACS` file ([gchen98/macs](https://github.com/gchen98/macs)). `VCF` is supported for PBWT (change `readMacs` to `readVcfGT` in `Snakefile`). Conversion to `VCF` to `MACS` is possibile using `script/vcf_to_macs.py`.

To test more queries quantity, edit `Snakefile` adding new path in `PANELS` list:

```python
PANELS = [
    "samples_x-sites_y/input.macs/q_z",
    "samples_x-sites_y/input.macs/q_w",
];
```

In `results` folder, the pipeline will provide some `CSV` files regarding `/usr/bin/time --verbose` output.
