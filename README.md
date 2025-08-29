# mupbwt_glimpse

```shell
./GLIMPSE2_simulate_bams_static --input-vcf target_unmasked_100_msprime_sim.bcf --read-length 150 --contig 1 -O reads100_150_default --thread 7

./split_reference/bin/GLIMPSE2_split_reference  --reference 10k/ref_panel_10k_msprime_sim.bcf --input-region 1:77-9999973 --output-region 1:77-9999973 --output 10k/ref_panel_10k

./phase/bin/GLIMPSE2_phase --bam-list unmasked_10k/reads100_150_default/all.txt --reference 10k/ref_panel_10k_1_77_9999973.bin --output unmasked_10k/ref_panel_10k_1_77_9999973_100_mu_mpsc.bcf --main 15 --burnin 5 --thread 5 --mupbwt --mupbwt-mpsc/mupbwt-smems

/concordance/bin/GLIMPSE2_concordance --input unmasked_10k/con_mpsc.txt --output unmasked_10k/concordance_mpsc --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 --thread 4 --gt-val

#1	unmasked_10k/frq.bcf	unmasked_10k/target_unmasked_100_msprime_sim.bcf	unmasked_10k/ref_panel_10k_1_77_9999973_100_orig.bcf

python conc.py unmasked_10k/concordance_orig.rsquare.grp.txt.gz unmasked_10k/orig.pdf
```

