# rlpbwt-test

```shell
curl -L https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh > miniconda.sh
bash miniconda.sh -b -f -p ./conda/
cd rlpbwt-test
../conda/bin/conda install -n base -c conda-forge mamba
source ../conda/bin/activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
source ../conda/bin/activate snakemake
../conda/bin/conda install -n snakemake -c conda-forge libcurl  openssl
snakemake run_after_checkpoint --cores 8
source ../conda/bin/deactivate
```
