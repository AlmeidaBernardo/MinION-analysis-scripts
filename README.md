# MinION-analysis-scripts
This software is licensed under the MIT License, but you must provide proper attribution to the original author (Bernardo Almeida) in any derivative works or commercial use.

This script aims to help in MinION sequences analysis.

Its functions serve 2 purposes:
    → function fastq_gz_merge → join many fastq.gz files in a single fastq.gz file to enable an easier analysis on programs like GenomeDetective, CZ ID, INSaFLU or others;
    → function read_bam → enable an alternative user-friendly way of analysing BAM files for dubious positions in the consensus sequence. Some positions may not be provided directly in the consensus sequence file provided by GenomeDetective and other programs, this happens due to low coverage, leading to low confidence in the correct nucleotide in such a position in the sequence. This positions can still be accessed and it is possible to see how many reads were aligned in such a position and what were the nucleotides in that position. To do so, I created this function read_bam (instructions on how to use it are described below)


How to use read_bam function?
    → Starting from GenomeDetective program, after analysing the merged fastq.gz file created with the function fastq_gz_merge, it is possible to download the BAM file (other analysis programs usually also provide this BAM file);
    → After downloading the BAM file, we need to transform it in a table:
        → To do so, go to the online program Galaxy (https://usegalaxy.org/);
        → Upload BAM file and a file with the reference genome sequence;
        → Seect the tool "Generate Pileup from bam dataset";
        → Then change "use built-in index" to "use one from history" and select the BAM file and the file with the reference genome sequence;
        → Run Tool;
    → It should be created a .tabular file. Dowlnoad it and use the function read_bam to correctly visualize it, without some of the columns that aren't usually needed. To keep those columns change the line "data = data.drop(columns=[0,2,5,6])" and the numbers of the line below accordingly

I am not a programmer, I'm "only" a master's student in Life Sciences so I can't guarantee that some bugs won't occur in specific cases, however, feel free to use this script and adapt it to your needs! Any doubt you might have or suggestion, feel free to contact me through the e-mail bsa.almeida2002@gmail.com or create a discussion post for the comunity.
Hope that I could help, good luck with your analysis of your MinION sequences!!

Thanks!

Bernardo Almeida

This software is licensed under the MIT License, but you must provide proper attribution to the original author (Bernardo Almeida) in any derivative works or commercial use.
