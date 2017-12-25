# speciespick

## Summary
Choose which species in RefSeq to analyze.  
Also organize other information on the species.

## Input
0. assembly_summary

    ```
    wget -O ./data/archaea_assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
    wget -O ./data/bacteria_assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    ```

0. NCBI taxonomy

    ```
    cd ${taxonomyDirec}
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip
    sed -e 's/\t|/,/g' nodes.dmp |sed -e 's/\t//g' > nodes.processed
    ```
    
    * the last the command change the field separator of nodes.dmp from "\t|\t" to ","

## Scripts
0. pick.py
    * read 2 assembly_summary file and update species table.
    * species table has a primary key restrictions on taxid
0. download.py
    * download fna, gff, and cds if needed, according to flow table. ftp filepath is already stored in species table.
