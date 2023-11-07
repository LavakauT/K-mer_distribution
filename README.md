# K-mer_distribution

Both MotifDiscovery or FIMO(source in MEME-Suit) are powerful tool to reveal which kinds of motif(K-mer) existing on custom sequences. For MotifDiscovery (https://github.com/ShiuLab/MotifDiscovery), enriched K-mers will be exported based on promoter regions between positive and negative genes. It means that K-mers are computed with Fisher's exact method and multiple testing correction. The name FIMO stands for 'Find Individual Motif Occurrences' (https://meme-suite.org/meme/doc/fimo.html?man_type=web). A groups of sequences will be sent to search occurrences of known motifs in corresponding database. FIMO apply its unique standard to determine whether one matching event is acceptatble to declare any of known motifs in database.

No matter which kind of tool you use, the preference from strand, similarity cut-off, and existence of reverse complement have not been reveald. Here this repository applied slide window with bin size = 25, window width = 100, and median/mean method to show motif(K-mer) distribution among this vairables.

[At kmer_distri.pdf](https://github.com/LavakauT/K-mer_distribution/files/13276838/At.kmer_distri.pdf)



# KINDLY REMINDER
Once you chose a low cut-off to grab all possible motif matching events, this will make your input dataframe become high dimension. Ｗe suggested not to run on your local computer because it consums lots of resources and times. HPC is a better place to run it.
