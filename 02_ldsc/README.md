# LDSC (Bulik-Sullivan et al., 2015)

Make sure you have installed [LDSC](https://github.com/bulik/ldsc) before running this.

## Installation

1. Clone the repository 

```
git clone https://github.com/bulik/ldsc.git
cd ldsc/
```

2. Install the dependencies using `conda`

```
conda env create --file environment.yml
conda activate ldsc
```

3. Move the scripts to somewhere easy to find

```
mkdir ~/.local/bin/
cp -r ldsc ~/.local/bin/
cd ~/.local/bin/
ln -s ldsc/ldpsc.py
ln -s ldsc/munge_sumstats.py
```

You can then add them to the `~/.bashrc` to update the `PATH`, but that's currently optional.

## Run Heritability

Simply run the script as follows. Model can be either `snipar` or `regenie`.

```
chmod +x run_heritability.sh
./run_heritability /path/to/sumstats.gz model /path/to/ldscores_chr@
```