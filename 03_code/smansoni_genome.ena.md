# Schistosoma mansoni V10 genome: ENA submission

### author: Stephen Doyle, stephen.doyle[at]sanger.ac.uk

- WSI S.mansoni genome is stored in ENA under the accession number GCA_000237925
    - current version is GCA_000237925.5
- need to update for the V10 genome


- ENA docs for submission
    - https://ena-docs.readthedocs.io/en/latest/update/assembly.html 


```bash






```




# convert GFF to EMBL format.
# using tool "EMBLmyGFF3" from https://github.com/NBISweden/EMBLmyGFF3
# note the locus tag needed to be preregistered with ENA, which I did via DNA Resease at Sanger.

EMBLmyGFF3 \
cercopithifilaria_johnstoni_annotation_SD210906.gff3 \
cercopithifilaria_johnstoni_genome_SD210906.fasta \
--data_class STD \
--accession \
--topology linear \
--molecule_type "genomic DNA" \
--transl_table 1  \
--species 'Cercopithifilaria johnstoni' \
--taxonomy INV \
--locus_tag CJOHNSTONI \
--project_id PRJEB47283 \
--author 'Stephen R. Doyle' \
-o cercopithifilaria_johnstoni_genome_SD210906.embl

# fix lineage - wasnt automatically picked by by tool due to the fact that C.johnstoni is a new species and taxon ID not searchable at the time of submission.
cat  cercopithifilaria_johnstoni_genome_SD210906.embl | sed 's/^OC.*/OC   Eukaryota\; Metazoa\; Ecdysozoa\; Nematoda\; Chromadorea\; Rhabditida;\nOC   Spirurina\; Spiruromorpha\; Filarioidea\; Onchocercidae\; Cercopithifilaria./g' | sed 's/transl_table=5/transl_table=1/g' > tmp; mv tmp cercopithifilaria_johnstoni_genome_SD210906.embl


gzip cercopithifilaria_johnstoni_genome_SD210906.embl