# RepeatMasker (Version 4.2.2) Repeat Annotation

Hello, everyone!
This is a guide to assemble your CycloneSEQ reads using _RepeatMasker_

The tutorial is built for [DCS Cloud](https://cloud.stomics.tech/) users

- Github repository: https://github.com/Nextomics/NextDenovo
- Contact email for [this](https://github.com/keen-laras/CycloneSEQ_NextDeNovo_GenomeAssembly) guide: kinanti.seraphina@gmail.com

## Tutorial
### 1. Image
Use this image to run software in online analysis
<img width="1724" height="402" alt="image" src="https://github.com/user-attachments/assets/8f400a2c-a629-4844-97eb-d47281d4e1b8" />

### 2. Run software
run [run_denovo.sh](https://github.com/keen-laras/RepeatAnnotation/blob/main/run_denovo.sh) for repeat annotation based on your de-novo assembly, or [run_library.sh](https://github.com/keen-laras/RepeatAnnotation/blob/main/run_library.sh) for repeat annotation based on published library for known repeats.

`bash run_denovo.sh` or `bash run_lib.sh`

database for RepeatMasker please see https://www.repeatmasker.org/

### 3. Mask output file
This pipeline results in a lowercase-repeat-fasta file. Run these commands to obtain the 'n-masked' for the genome annotation input.

step 1: create a .bed file containing all the repeat coordinates
`perl extract_lowercase_bed.pl /path/to/repeatmasker/{sample}.fa.masked {output}.nmasked.bed`
step 2: combine and sort both .bed files from libraries
`cat /denovo/{output}.nmasked.bed /dfam_lib/{output}.nmasked.bed > {output}_combined_nmasked.raw.bed
sort -k1,1 -k2,2n {output}_combined_nmasked.raw.bed > {output}_combined_nmasked.sorted.bed`
step 3: merge coordinates from both .bed 
`/share/app/bedtools/2.29.2/bin/mergeBed -i {output}_combined_nmasked.sorted.bed > {output}_combined_nmasked.merged.bed`
step 4: create masked fasta file
`/share/app/bedtools/2.29.2/bin/maskFastaFromBed -mc N -fi /path/to/{sample}.fa -fo {output}.FINAL.masked.fa -bed {output}_combined_nmasked.merged.bed`
