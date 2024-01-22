## Database download

The most convinient way to download the database is to download the pre-indexed database (Current version indexed from december 22nd 2023 cgMLST scheme):

```
wget
tar -zxvf cgMLST_db.tar.gz
```

Alternatively, if you want the newest cgMLST schemes, you can follow these step and setup the database yourself (The database indexing may take several hours, so consider doing it in a detached screen session):

1.) Download all the cgMLST schemes from the cgMLST website (https://www.cgmlst.org/ncs). 
For Burkholderia mallei only download the FLI scheme and make sure to rename the folder to Burkholderia_mallei_cgMLST_alleles. 
2.) Place all the folders folder and create a tar ball of all the folders.
```
tar -zcvf cgMLST_schemes.tar.gz cgMLST_schemes
```
3.) run the setup_databases.py script found in the scripts folder. This will create the database and index it.
```
python3 setup_databases.py -i cgMLST_schemes.tar.gz -o cgMLST_db
```

