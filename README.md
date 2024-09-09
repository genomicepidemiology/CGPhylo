# cgPhylo

cgPhylo is a tool for determining phylogeny between bacterial isolates for which there are core genome MLST schemas ([cgmlst.org](https://https://www.cgmlst.org/)). 
Phylogeny is determined by comparing the multiple alignments of conserved core genes shared by the input samples.

## Installation

To install cgPhylo run:

```bash
conda install -c genomicepidemiology cgphylo
```

## Database download

The most convenient way (by far) to use this tool is to download the pre-indexed database (Current version indexed from december 22nd 2023 cgMLST scheme):

```
wget https://cge.cbs.dtu.dk/services/cgphylo/cgmlst_db.tar.gz
tar -zxvf cgmlst_db.tar.gz
```

Alternatively, if you want the newest cgMLST schemes, you can follow these step and setup the database yourself (The database indexing may take several hours, so consider doing it in a detached screen session):

1.) Download all the cgMLST schemes from the cgMLST website (https://www.cgmlst.org/ncs). 
For Burkholderia mallei only download the FLI scheme and make sure to rename the folder to Burkholderia_mallei_cgMLST_alleles. 
2.) Place all the folders folder and create a tar ball of all the folders.
```
tar -zcvf cgMLST_schemes.tar.gz cgMLST_schemes
```
3.) run the setup_databases.py script found in the scripts folder. This will create the database and index it. Note: If this doesn't work, it may be because of a new update to the cgMLST databases which we have accounted for. Please raies an issue on Github in that case, and we will fix it. 
```
python3 scripts/setup_databases.py -i cgMLST_schemes.tar.gz -o cgMLST_db
```

## Usage

```bash
cgphylo --nanopore /complete/path/to/input_files/*.fastq.gz --o any_output_name --threads <int, default:4> --db_dir /path/to/cgmlst_db
```

## License

Apache License 2.0

## Authors
Malte B. Hallgren