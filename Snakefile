from multiprocessing import cpu_count
full_core=cpu_count()
half_core=int(full_core/2)

configfile: "config.yaml"

rule curing:
    input:
        config["short_reads"]
    output:
        expand("{out_dir}/trim_reads/trim_reads.fastq.gz", out_dir=config["out_dir"])
    threads: 
        cpu_count()
    log: expand("{out_dir}/logs/triming_reads/UrQt.log", out_dir=config["out_dir"])
    message: "triming the short-reads with UrQt"
    shell:
        "(src/UrQt/UrQt --m {threads} --in {input} --gz --min_read_size 75 --out {output} ) 2> {log}"

rule fast_dnapipe:
    input:
        reads=rules.curing.output
    output:
        annoted=expand("{out_dir}/fast_dnapipete_output/Annotation/annoted.fasta", out_dir=config["out_dir"]),
        unannoted=expand("{out_dir}/fast_dnapipete_output/Annotation/unannoted.fasta", out_dir=config["out_dir"])
    threads: cpu_count()
    log: expand("{out_dir}/logs/fast_dnapipete/dnapipete.log", out_dir=config["out_dir"])
    message: "first run of dnapipete to try to catch mitochondrion/contamination sequence and exclud them later"
    params:
        genome_size=config["genome_size"],
        sampling_size=config["sampling_size"],
        out_dir="fast_dnapipete_output",
        direct=expand("{out_dir}",out_dir=config["out_dir"])
    shell:
        "printf '#! /bin/bash\ncd /opt/dnaPipeTE\nexport LIBDIR=/mnt2/RepBase/new_lib\npython3 dnaPipeTE.py -input /mnt/trim_reads/trim_reads.fastq.gz" 
        " -output /mnt/{params.out_dir} -RM_lib /mnt2/RepBase/new_lib/RepeatMasker.lib -genome_size {params.genome_size} -genome_coverage {params.sampling_size} "
        " -sample_number 2 -RM_t 0.2 -cpu {threads}' > {params.direct}/{params.out_dir}/script_dnapipete.sh && chmod u+x {params.direct}/{params.out_dir}/script_dnapipete.sh && "
        " ( singularity exec --bind {params.direct}:/mnt,$PWD:/mnt2 $PWD/src/dnapipete.img /mnt/{params.out_dir}/script_dnapipete.sh ) 2> {log}"

rule concat_trinity:
    input:
        annoted=rules.fast_dnapipe.output.annoted,
        unannoted=rules.fast_dnapipe.output.unannoted
    output:
        expand("{out_dir}/fast_dnapipete_output/Annotation/unannoted_and_annoted.fasta", out_dir=config["out_dir"])
    message: "concatening the two fastdnapipe output annoted and unanotted"
    shell:
        "cat {input.annoted} {input.unannoted} > {output}"

rule extract_mito:
    input:
        DB_mito="database_mito_conta/mitochondrial_sequence_refseq.fasta.gz",
        trinity_out=rules.concat_trinity.output
    output:
        bam=temp(expand("{out_dir}/mitochondrial_map.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/mitochondrial_map_sort.bam", out_dir=config["out_dir"])),
        premito_fasta=temp(expand("{out_dir}/mitochondrial_sequences_extract_pre_rename.fasta", out_dir=config["out_dir"])),
        mito_fasta=expand("{out_dir}/mitochondrial_sequences_extract.fasta", out_dir=config["out_dir"]),
        mito_tab=temp(expand("{out_dir}/correspond_map_ref.tab", out_dir=config["out_dir"]))
    threads: half_core
    log: expand("{out_dir}/logs/extract_mito_seq/mapping.log", out_dir=config["out_dir"])
    message: "mapping the trinity output to the mitoDB to extract the mitochondrial sequence "
    conda:
        "env/minimap2.yaml"
    params:
        species=config["species"],
        acc_num=config["acc_num"],
        tmp_file=expand("{out_dir}/tmp_mito",out_dir=config["out_dir"])
    shell:
        "( minimap2 -a --split-prefix={params.tmp_file} -t {threads} {input.DB_mito} {input.trinity_out} | samtools view  -F 4  -Sb - > {output.bam} && "
        "samtools sort -@ {threads} -o {output.bam_sort} {output.bam} && samtools index {output.bam_sort} && "
        "samtools fasta -@ {threads}  {output.bam_sort} > {output.premito_fasta} && "
        "samtools view -F 2304 {output.bam_sort} | cut -f1,3 > {output.mito_tab} && "
        "python src/rename_seq.py {output.premito_fasta} {output.mito_tab} {params.species} {params.acc_num} {output.mito_fasta} && "
        "cat {output.mito_fasta} >> database_mito_conta/found_mitochondrial_sequence.fasta ) 2> {log}"

rule extract_conta:
    input:
        trinity_out=rules.concat_trinity.output,
        DB_conta="database_mito_conta/conta_sequence_refseq.fasta.gz"
    output:
        bam=temp(expand("{out_dir}/conta_map.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/conta_map_sort.bam", out_dir=config["out_dir"])),
        conta_fasta=expand("{out_dir}/conta_sequences_extract.fasta", out_dir=config["out_dir"])
    threads: half_core
    log: expand("{out_dir}/logs/extract_conta_seq/mapping.log", out_dir=config["out_dir"])
    message: "mapping the trinity output to the contaDB to extract the contamination sequence "
    params:
        tmp_file=expand("{out_dir}/tmp_conta",out_dir=config["out_dir"])
    conda:
        "env/minimap2.yaml"
    shell:
        "( minimap2 -a --split-prefix={params.tmp_file} -t {threads} {input.DB_conta} {input.trinity_out} | samtools view  -F 4  -Sb - > {output.bam} && "
        "samtools sort -@ {threads} -o {output.bam_sort} {output.bam} && samtools index {output.bam_sort} && "
        "samtools fasta -@ {threads}  {output.bam_sort} > {output.conta_fasta} ) 2> {log}"

rule cleaning_reads:
    input:
        mito=rules.extract_mito.output.mito_fasta,
        conta=rules.extract_conta.output.conta_fasta,
        reads=rules.curing.output
    output:
        clean_reads=expand("{out_dir}/trim_reads/clean_trim_reads.fastq", out_dir=config["out_dir"]),
        to_extract=temp(expand("{out_dir}/to_extract.fasta", out_dir=config["out_dir"])),
        bam=temp(expand("{out_dir}/to_extract.bam", out_dir=config["out_dir"])),
        bam_sort=temp(expand("{out_dir}/to_extract_sort.bam", out_dir=config["out_dir"]))

    threads: cpu_count()
    log: expand("{out_dir}/logs/cleaning_reads/mapping.log", out_dir=config["out_dir"])
    message: "mapping the original reads to the mito/conta seq found in the first dnapipete run to exclude them"
    conda:
        "env/minimap2.yaml"
    shell:
        "python src/extract_seq.py {input.mito} {input.conta} {output.to_extract} &&"
        " ( minimap2 -ax sr {output.to_extract} {input.reads} -t {threads} | samtools view -f 4 -Sb - > {output.bam} && "
        " samtools sort -@ {threads} -o {output.bam_sort} {output.bam} && samtools index {output.bam_sort} && "
        " samtools fastq -@ {threads} {output.bam_sort} > {output.clean_reads}) 2> {log}"

rule final_dnapipe:
    input:
        rules.cleaning_reads.output.clean_reads
    output:
        expand("{out_dir}/final_dnapipete_output/Trinity.fasta", out_dir=config["out_dir"])
    threads: cpu_count()
    log: expand("{out_dir}/logs/final_dnapipete/dnapipete.log", out_dir=config["out_dir"])
    message: "final run of dnapipete"
    params:
        genome_size=config["genome_size"],
        sampling_size=config["sampling_size"],
        out_dir="final_dnapipete_output",
        direct=expand("{out_dir}",out_dir=config["out_dir"])
    shell:
        "printf '#! /bin/bash\ncd /opt/dnaPipeTE\nexport LIBDIR=/mnt2/RepBase/new_lib\npython3 dnaPipeTE.py -input /mnt/trim_reads/clean_trim_reads.fastq" 
        " -output /mnt/{params.out_dir} -RM_lib /mnt2/RepBase/new_lib/RepeatMasker.lib -genome_size {params.genome_size} -genome_coverage {params.sampling_size} "
        " -sample_number 2 -RM_t 0.2 -cpu {threads}' > {params.direct}/{params.out_dir}/script_dnapipete.sh && chmod u+x {params.direct}/{params.out_dir}/script_dnapipete.sh && "
        " ( singularity exec --bind {params.direct}:/mnt,$PWD:/mnt2 $PWD/src/dnapipete.img /mnt/{params.out_dir}/script_dnapipete.sh ) 2> {log}"


rule landscape:
    input:
        rules.final_dnapipe.output
    output:
        expand("{out_dir}/final_dnapipete_output/dnaPipeTE_landscapes_subclass.pdf", out_dir=config["out_dir"])
    message: "ploting an histogram of the blastn divergence between raw reads (TE copies in the genomes) and their consensus sequences assembled in the file Trinity.fasta. The script plots only putative TE sequences among the subclasses 'LINE', 'SINE', 'LTR', 'DNA', 'RC' and 'Unknown' (a.k.a. 'NA') "
    log:  expand("{out_dir}/logs/final_dnapipete/dnapipete_landscapes.log", out_dir=config["out_dir"])
    params:
        in_dir=expand("{out_dir}/final_dnapipete_output/", out_dir=config["out_dir"])
    conda:
        "env/dnaPT_utils.yaml"
    shell:
        "(src/dnaPT_utils/dnaPT_landscapes.sh -I {params.in_dir} ) 2> {log}"


rule charts:
    input:
        rules.final_dnapipe.output
    output:
        expand("{out_dir}/final_dnapipete_output/dnaPipeTE_charts.pdf", out_dir=config["out_dir"])
    message: "Processing an output folder of dnaPipeTE and producing 3 graphs \nPiechart 1 : relative proportions of the different repeat categories\nPiechart 2 : proportion of the different repeat categories relative to the total genome\nBarplot : genomic proportion of each dnaPipeTE contig and the associated classification"
    log:  expand("{out_dir}/logs/final_dnapipete/dnapipete_charts.log", out_dir=config["out_dir"])
    params:
        in_dir=expand("{out_dir}/final_dnapipete_output/", out_dir=config["out_dir"])
    conda:
        "env/dnaPT_utils.yaml"
    shell:
        "(src/dnaPT_utils/dnaPT_charts.sh -I {params.in_dir} ) 2> {log}"

rule all:
    input:
        rules.final_dnapipe.output,
        rules.landscape.output,
        rules.charts.output

rule initialise_dnaputil:
    input:
        "Snakefile"
    output:
        "src/dnaPT_utils/dnaPT_charts.sh"
    message: "downloading dnapipe_utils from github"
    shell:
        "cd src && git clone https://github.com/clemgoub/dnaPT_utils.git"

rule initialise_UrQt:
    input:
        "Snakefile"
    output:
        "src/UrQt/UrQt"
    message: "downloading UrQT from github"
    shell:
        "cd src && git clone https://github.com/l-modolo/UrQt && cd UrQt && make"

rule add_RepBase:
    input:
        h5="RepBase/tmp_repbase/Dfam.h5",
        embl="RepBase/tmp_repbase/RMRBSeqs.embl"
    output:
        "RepBase/new_lib/RepeatMaskerLib.h5"
    message: "creating a new RepeatMasker library with RepBase"
    log: "add_RepBase.log"
    params:
        h5="/mnt2/Dfam.h5",
        embl="/mnt2/RMRBSeqs.embl"
    shell:
        "singularity exec --bind RepBase:/mnt2 src/dnapipete.img /mnt2/add_Repbase.sh"

rule NSDPY_database:
    input:
        #conta="database_mito_conta/conta_sequence_refseq.search",
        mito="database_mito_conta/mitochondrial_sequence_refseq.search"
    output:
        #conta="database_mito_conta/conta_sequence_refseq.fasta.gz",
        mito="database_mito_conta/mitochondrial_sequence_refseq.fasta.gz"
    message: "creating mito/conta database"
    params:
        #conta_u="database_mito_conta/conta_sequence_refseq.fasta",
        mito_u="database_mito_conta/mitochondrial_sequence_refseq.fasta",
        #cmd_conta="(Bacteria[orgn] OR Archaea[orgn] OR Fungi[orgn]) AND (reference_genome[filter] OR representative_genome[filter])",
        cmd_mito="srcdb_refseq[PROP] AND (chloroplast[filter] OR mitochondrion[filter])"
    shell:
        #"nsdpy -r '{params.cmd_conta}' && cp NSDPY_results/*/sequences.fasta {params.conta_u} && rm -rf NSDPY_results && "
        " nsdpy -r '{params.cmd_mito}' && cp NSDPY_results/*/sequences.fasta {params.mito_u} && rm -rf NSDPY_results && "
        #" gzip {params.conta_u} &&" 
        " gzip {params.mito_u}"

rule initialise:
    input:
        rules.initialise_dnaputil.output,
        rules.initialise_UrQt.output,
        rules.add_RepBase.output,
        rules.NSDPY_database.output.mito

