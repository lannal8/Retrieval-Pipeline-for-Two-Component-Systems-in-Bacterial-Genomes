#Replace GCA_000022305.1 with your genome accession 
genome_accession="GCA_000022305.1"

#Searching for HATPase_c, HisKA, DNA_gyrase, Response_reg domains
hmmscan --domtblout 03_hmmscan/${genome_accession}_HATPase.out --cpu 2 --domE 1e-5 -E 1e-5 02_hmm_model_build/HATPase.hmm /data/xzliu/DATA_base/Bacteria_genome_gbk_update/all/${genome_accession}/protein.faa 1>/dev/null 2>&1
hmmscan --domtblout 03_hmmscan/${genome_accession}_HisKA.out --cpu 2 --domE 1e-5 -E 1e-5 02_hmm_model_build/HisKA.hmm /data/xzliu/DATA_base/Bacteria_genome_gbk_update/all/${genome_accession}/protein.faa 1>/dev/null 2>&1
hmmscan --domtblout 03_hmmscan/${genome_accession}_DNA_gyrase.out --cpu 2 --domE 1e-5 -E 1e-5 02_hmm_model_build/DNA_gyrase.hmm /data/xzliu/DATA_base/Bacteria_genome_gbk_update/all/${genome_accession}/protein.faa 1>/dev/null 2>&1
hmmscan --domtblout 03_hmmscan/${genome_accession}_RR.out --cpu 2 --domE 1e-5 -E 1e-5 02_hmm_model_build/RR.hmm /data/xzliu/DATA_base/Bacteria_genome_gbk_update/all/${genome_accession}/protein.faa 1>/dev/null 2>&1

#Filtering matches with i-Evalue < 1e-5
awk '$13 <= 1e-5 {print $0}' 03_hmmscan/${genome_accession}_HATPase.out > 03_hmmscan/${genome_accession}_HATPase.filtered.out
awk '$13 <= 1e-5 {print $0}' 03_hmmscan/${genome_accession}_HisKA.out > 03_hmmscan/${genome_accession}_HisKA.filtered.out
awk '$13 <= 1e-5 {print $0}' 03_hmmscan/${genome_accession}_DNA_gyrase.out > 03_hmmscan/${genome_accession}_DNA_gyrase.filtered.out
awk '$13 <= 1e-5 {print $0}' 03_hmmscan/${genome_accession}_RR.out > 03_hmmscan/${genome_accession}_RR.filtered.out

#Generating lists of protein IDs that passed filtering
awk '!/^#/ {print $4}' 03_hmmscan/${genome_accession}_HATPase.filtered.out |uniq > 04_list_and_overlap/${genome_accession}_HATPase.list
awk '!/^#/ {print $4}' 03_hmmscan/${genome_accession}_HisKA.filtered.out |uniq > 04_list_and_overlap/${genome_accession}_HisKA.list
awk '!/^#/ {print $4}' 03_hmmscan/${genome_accession}_DNA_gyrase.filtered.out |uniq > 04_list_and_overlap/${genome_accession}_DNA_gyrase.list
awk '!/^#/ {print $4}' 03_hmmscan/${genome_accession}_RR.filtered.out |uniq > 04_list_and_overlap/${genome_accession}_RR.list

#First set of candidate Histidine Kinases: Proteins containing both HATPase_c and HisKA domains
grep -F -f 04_list_and_overlap/${genome_accession}_HisKA.list 04_list_and_overlap/${genome_accession}_HATPase.list > 04_list_and_overlap/${genome_accession}_HATPaseANDHisKA.list

#Second set of candidate Histidine Kinases: Proteins containing HATPase_c but lacking HisKA and DNA_gyrase domains
#These require InterProScan analysis to verify they are not single-domain HATPase_c proteins
grep -F -x -v -f 04_list_and_overlap/${genome_accession}_HisKA.list 04_list_and_overlap/${genome_accession}_HATPase.list > 04_list_and_overlap/${genome_accession}_HATPaseANDNOTHisKA.list
grep -F -v -f 04_list_and_overlap/${genome_accession}_DNA_gyrase.list 04_list_and_overlap/${genome_accession}_HATPaseANDNOTHisKA.list > 04_list_and_overlap/${genome_accession}_HATPaseANDNOTHisKAANDNOTDNA_gyrase.list
seqkit grep -f 04_list_and_overlap/${genome_accession}_HATPaseANDNOTHisKAANDNOTDNA_gyrase.list /data/xzliu/DATA_base/Bacteria_genome_gbk_update/all/${genome_accession}/protein.faa > 05_HATPaseANDNOTHisKAANDNOTDNA_gyrase_fa/${genome_accession}.fa
/data/xzliu/Soft/interproscan-5.46-81.0/interproscan.sh -appl Pfam -i 05_HATPaseANDNOTHisKAANDNOTDNA_gyrase_fa/${genome_accession}.fa -o 06_interproscan/${genome_accession}.tsv -cpu 1 -T temp/ -f tsv --disable-precalc
python interproscan_result_get.py 06_interproscan/${genome_accession}.tsv 07_interproscan_list/${genome_accession}.HATPaseANDNOTHisKAANDNOTDNA_gyrase_interproscan.list

#Merging two sets of candidate Histidine Kinases
cat 04_list_and_overlap/${genome_accession}_HATPaseANDHisKA.list 07_interproscan_list/${genome_accession}.HATPaseANDNOTHisKAANDNOTDNA_gyrase_interproscan.list > 08_HK_and_RR_list/${genome_accession}_HK.list

#Candidate Response Regulators: Proteins containing Response_reg domain
cp 04_list_and_overlap/${genome_accession}_RR.list 08_HK_and_RR_list/${genome_accession}_RR.list

#Classifying hybrid proteins based on the sequential order of HATPase_c and Response_reg domains into Hybrid HK and Hybrid RR
python HHK_or_HRR.py -i ${genome_accession}
