# MetagenomicsV3_Quality
A open access python script to process metagenomics data generated in pipelines of https://doi.org/10.1590/1678-4685-GMB-2023-0015 and https://doi.org/10.3389/fmicb.2022.1002963.

# Quick Tutorial
This Script was built to process the output file of the Slurm Batch script used in the mentioned papers.

To run the analysis, is recommended to have Anaconda (https://www.anaconda.com/download) or Mamba (https://mamba.readthedocs.io/en/latest/) installed in your local machine.

Acess your terminal and activate conda
```bash
conda activate
```
Copy the present repository
```bash
git clone https://github.com/matheus-cosentino/MetagenomicsV3_Quality.git
cd MetagenomicsV3_Quality
```
Run the following code to check the help of the script
```bash
python Illumina_quality_check.py --help
```

```output
usage: Illumina_quality_check.py [-h] --slurm SLURM --output OUTPUT

Process a Slurm output file of Metagenomics.

optional arguments:
  -h, --help       show this help message and exit
  --slurm SLURM    Path to the Slurm output file (.out) generated in MetagenomicV3. Example: MetaGen3_173518.out
  --output OUTPUT  Path to the output file (.csv) with processed data by this script . Example: MetaGen3_173518.csv
```


Run the following code and check all the data you processed in the Output file "MetaGen3_173517.csv"
```bash
python Illumina_quality_check.py --slurm MetaGen3_173517.out --output MetaGen3_173517.csv
```

# Citation
If you find this script useful, please consider citing this GitHub repository.
