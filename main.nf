#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Mandatory Params
params.Samples = "" // Sets the input sample table file 
params.Organism = "" // Sets the organism for analize
// Optional Params
params.Genome = ""
params.Bowtie_index = ""
params.Adapters_fasta = ""
params.mod_FastP = "" // Tell if is necessary a poly G trim or another FastP Flag
params.ConditionsFile = ""
params.Control_name = "Control"
params.Case_name = "Case"
params.PerCHR = false
params.Perl = false
params.DESeq2 = false
// ShortCat Original parameters
params.max_N_alig = 10
params.mismatch = 1
params.cr = 5
params.Strand = "no"

process FASTQC {
    publishDir "${params.outdir}/FASTQC", mode:'copy', pattern: '*.html' // only saves the html files to FASTQC folder    
    tag "FASTQC on $sample_id"
    
    input: // section to declare the input of the process
    tuple val(sample_id), path(reads)

    output: // section to declare the output of the process  
    tuple val(sample_id), path("*.html"), emit: html 
    path "*zip", emit: multiqc_input // sets a special output to Multiqc tool, in order to get a summary of all sequenciation runs

    script: // Section to declare the code in bash or a specific programing language like Python or R
    """ 
    fastqc $reads
    """// Runs Fastqc for each fastq in order to check the quality of sequenciation (bash command)
}
process MULTIQC {   
    publishDir params.outdir, mode:'copy'  

    input:
    file archivos
    val(name)

    output:
    path "*_multiqc_report.html" // defines the html file as the required output of the process
    
    script:
    """    
    multiqc $archivos -n '$name'_multiqc_report.html
    """ // (bash command) calls multiqc tool and uses all the zip files in the folder with the data of quality
}

process Get_adapters {
    cpus 6

    input:
    tuple val(sample_id), path(fastqc_file)

    output:
    path("Adapters.${sample_id}.fasta")

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    Overrepresented_seqs = pd.read_html("${fastqc_file}")[1]
    Overrepresented_seqs = Overrepresented_seqs.loc[~Overrepresented_seqs["Possible Source"].str.lower().str.contains("no hit")].reset_index()
    # Función para convertir la columna a formato FASTA
    def columna_a_fasta(df, nombre_columna,nombre_columna2, archivo_salida):
        with open(archivo_salida, 'w') as f:
            for i in range(0,len(df)):
                f.write(f'>'+df[nombre_columna][i]+" | Sequence N° "+ str(i+1) +'\\n'+df[nombre_columna2][i]+'\\n')

    # Llamar a la función y pasar el nombre de la columna y el archivo de salida
    columna_a_fasta(Overrepresented_seqs, 'Possible Source','Sequence', 'Adapters.${sample_id}.fasta')
    """
}
process Concat_adapters {
    conda "/media/storage2/software/anaconda3/envs/ShortCat"

    input:
    path(Adapters_file)

    output:
    path("all.Adapters.fasta")

    script:
    """
    cat Adapters.*.fasta > all.Adapters.fasta
    """
}
process Remove_duplicated_adapters1 {
    conda "/media/storage2/software/anaconda3/envs/ShortCat"

    input:
    path(Adapters_file)

    output:
    path("Dedup.Adapters.fasta")

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import os
    os.system("cat Adapters.*.fasta > all.Adapters.fasta")
    def eliminar_duplicados_fasta(archivo_entrada, archivo_salida):
        secuencias_unicas = {}
        for record in SeqIO.parse(archivo_entrada, "fasta"):
            # Utiliza la secuencia como clave y guarda la primera cabecera encontrada
            if str(record.seq) not in secuencias_unicas:
                secuencias_unicas[str(record.seq)] = record.id
        with open(archivo_salida, 'w') as f:
            contador=0
            for secuencia, cabecera in secuencias_unicas.items():
                contador+=1
                f.write(f'>{cabecera}_Seq{str(contador)}\\n{secuencia}\\n')

    # Llamar a la función con el nombre de tu archivo FASTA
    eliminar_duplicados_fasta('all.Adapters.fasta', 'Dedup.Adapters.fasta')
    """
}
process Remove_duplicated_adapters2 {
    conda "/media/storage2/software/anaconda3/envs/ShortCat"

    input:
    path(Adapters_file)
    path(External_fasta)

    output:
    path("Dedup.Adapters.fasta")

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    import os
    os.system("cat Adapters.*.fasta ${External_fasta} > all.Adapters.fasta")
    def eliminar_duplicados_fasta(archivo_entrada, archivo_salida):
        secuencias_unicas = {}
        for record in SeqIO.parse(archivo_entrada, "fasta"):
            # Utiliza la secuencia como clave y guarda la primera cabecera encontrada
            if str(record.seq) not in secuencias_unicas:
                secuencias_unicas[str(record.seq)] = record.id
        with open(archivo_salida, 'w') as f:
            contador=0
            for secuencia, cabecera in secuencias_unicas.items():
                contador+=1
                f.write(f'>{cabecera}_Seq{str(contador)}\\n{secuencia}\\n')

    # Llamar a la función con el nombre de tu archivo FASTA
    eliminar_duplicados_fasta('all.Adapters.fasta', 'Dedup.Adapters.fasta')
    """
}
process FASTP {    // Tool for trimming 
    publishDir "${params.outdir}/01-FASTP", mode:'copy' 
    maxForks 6    // allows 6 instances of this process in parallel
    cpus 10
    //publishDir "${params.outdir}/Trimmed", mode:'copy', pattern: '*.fastq.gz' // Saves trimmed fastqs to trimmed folder
    tag "FASTP on $sample_id"    
    
    input:
    tuple val(sample_id), path(reads)
    each path(Adapters)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    
    script:     
    """
    fastp --thread $task.cpus $params.mod_FastP --in1 ${reads[0]} --adapter_fasta ${Adapters} --length_required 16 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 30 --trim_front1=5 --trim_tail1=5 --out1 '$sample_id'_trimmed_R1.fastq.gz
    """
}
process Bowtie {
    publishDir "${params.outdir}/02-Bowtie", mode:"copy", pattern: "*.sam"
    maxForks 6
    cpus 10
    tag "Bowtie2 on $sample_id"  

    input:
    tuple val(sample_id), path(reads)
    each path(Index)

    output:
    path("${sample_id}_modified.sam"), emit: mapped_reads

    script:
    """
    bowtie2 -p $task.cpus --mm -N ${params.mismatch} -x ${Index}/${params.Organism} -U ${reads} -S ${sample_id}_${params.mismatch}_m${params.max_N_alig}.sam 2> BOWTIE_REPORT.${sample_id}.txt
    grep -v "@" ${sample_id}_${params.mismatch}_m${params.max_N_alig}.sam | cut -f2-14 | awk '{print "${sample_id}_"NR"\t"\$LINE}' > ${sample_id}_modified.sam
    """
    // --best --strata -m ${params.max_N_alig}
}
process BowtieIndexBuild {
    storeDir "${params.outdir}/02-Bowtie"
    cpus 60

    input:
    path(genome_fasta)

    output:
    path("${params.Organism}/"), emit: index

    script:
    """
    mkdir ${params.Organism}
    bowtie2-build --threads $task.cpus ${genome_fasta}  ${params.Organism}/${params.Organism}
    """
}
process process_SAMs {
    maxForks 1

    input:
    path(SAMs)

    output:
    path("all.*.sam"), emit: SAM

    script:
    if (params.PerCHR == false)
        """ 
        cat *_modified.sam > all.reads.sam
        """
    else
        """ 
        #!/bin/bash
        cat *_modified.sam > all.reads.sam
        for i in \$(awk '!seen[\$3]++ {print \$3}' all.reads.sam | grep -v "*" | sort -n); do echo "awk '{if(\\\$3 == \\"\$i\\") print \\\$0}' all.reads.sam > all.\$i.reads.sam" >> comands.txt; done
        xargs -a comands.txt -d '\\n' -I {} -P 30 sh -c "{}"
        rm all.reads.sam
        """
} 
process ShortCat {
    publishDir "${params.outdir}/04-ShortCat", mode:"copy", pattern: '*.counts'
    publishDir "${params.outdir}/04-ShortCat/BEDs", mode:"copy", pattern: '*.bed'
    maxForks 10
    cpus 6

    input:
    path(SAM)
    val(Indices)

    output:
    path("all*.counts"), emit: counts
    path("*.bed")

    script:
    """
    #!/usr/bin/env perl
    use strict;
    use warnings;
    use Getopt::Long;
    my \$the_indexes = "${Indices}";
    chomp \$the_indexes;
    \$the_indexes =~ s/ //;
    \$the_indexes =~ s/\\]//;
    \$the_indexes =~ s/\\[//;
    my \$File_name = "${SAM}";
    chomp \$File_name;
    \$File_name =~ s/.sam//;
    \$File_name =~ s/all.//;
    print "\$File_name";
    my \$conditions = \$the_indexes;
    open(PRI,"> all.\$File_name.counts");
    my @condition = split(",",\$conditions);
    system("sam2bed --max-mem 70G < ${SAM} > ShortCat.\$File_name.bed");
    my \$strand_specific = "no";
    if (\$strand_specific eq "yes"){
        system("mergeBed -s -i ShortCat.\$File_name.bed -c 4,6 -o distinct -delim \\";\\" > ShortCat_merge.\$File_name.bed");
    } else {
        system("mergeBed -i ShortCat.\$File_name.bed -c 4,6 -o distinct -delim \\";\\" > ShortCat_merge.\$File_name.bed");
    }
    open(ACE,"ShortCat_merge.\$File_name.bed") or die \$!;
    my \$count = 1;
    print PRI "Chrom\\tStart\\tEnd\\tContig\\tLength\\tstrand\\t#reads";
    foreach my \$condi (@condition) {
        print PRI "\\t\$condi";
    }
    print PRI "\\n";
    my \$indicador = 0;
    while(<ACE>){
        chomp;
        my @col = split(/\\t/,\$_);
        my @infoAf = split(/\\;/,\$col[3]);
        my %counts = ();
        foreach my \$condi (@condition) {
        \$counts{\$condi} = 0;		#cambiar a 1 para sumarle 1 a todas las cuentas por defecto
        }
        my \$cantidad = @infoAf;
        foreach my \$ests (@infoAf) {
            \$indicador = 0;
            foreach my \$condi (@condition) {
                if(\$ests =~ /^\$condi/) {
                    \$counts{\$condi}++;
                    \$indicador = 1;
                }
            }
    #        if (\$indicador == 0) {
    #            die "Tipo no encontrado : \$_\\tline: \$count";
    #        }
        }
        my \$largo = \$col[2] - \$col[1];
        print PRI "\$col[0]\\t\$col[1]\\t\$col[2]\\t\$col[0]_smallRNA_\$count\\t\$largo\\t\$col[4]\\t\$cantidad";
        foreach my \$condi (@condition) {	
            print PRI "\\t\$counts{\$condi}";
        }
        print PRI "\\n";
        \$count++;
    }
    close(PRI); 
    """
}
process ShortCatPython {
    publishDir "${params.outdir}/04-ShortCatPython", mode:"copy", pattern: '*.counts'
    publishDir "${params.outdir}/04-ShortCatPython/BEDs", mode:"copy", pattern: '*.bed'
    maxForks 1
    cpus 6

    input:
    path(SAM)
    val(Indices)

    output:
    path("all*.counts"), emit: counts
    path("*.bed")

    script:
    """
    #!/usr/bin/env python
    import os
    import subprocess

    # Dividir el nombre del archivo de entrada para crear el nombre del archivo de salida si no se proporciona
    filename = "${SAM}".replace(".reads.sam","")
    fileout = filename + ".counts"

    # Dividir las condiciones por coma
    conditions = "${Indices}".replace("[","").replace("]","").replace(" ","").split(",")

    # Ejecutar el comando sam2bed
    subprocess.run("sam2bed --max-mem 70G < ${SAM} > " + filename + ".bed", shell=True)

    # Ejecutar el comando mergeBed dependiendo de la especificación de hebra
    if "${params.Strand}" == 'yes':
        subprocess.run("mergeBed -s -i " + filename + ".bed -c 4,6 -o distinct -delim \\";\\" > " + filename + "_merge.bed", shell=True)
    else:
        subprocess.run("mergeBed -i " + filename + ".bed -c 4,6 -o distinct -delim \\";\\" > " + filename + "_merge.bed", shell=True)


    # Asumiendo que las condiciones y el nombre del archivo de entrada ya están definidos
    filein = filename + "_merge.bed"

    # Abrir el archivo de entrada
    try:
        with open(filein, 'r') as ace:
            # Inicializar el contador y el indicador
            count = 1
            indicador = 0
            # Crear el archivo de salida si no se ha definido
            fileout = filein.split('.')[0] + ".count" if not fileout else fileout
            # Abrir el archivo de salida para escribir
            with open(fileout, 'w') as pri:
                # Escribir los encabezados en el archivo de salida
                headers = ["Chrom", "Start", "End", "Contig", "Length", "strand", "#reads"] + conditions
                pri.write("\\t".join(headers) + "\\n")
                # Leer cada línea del archivo de entrada
                for line in ace:
                    line = line.strip()
                    col = line.split('\\t')
                    infoAf = col[3].split(';')
                    counts = {condi: 0 for condi in conditions}  # cambiar a 1 para sumarle 1 a todas las cuentas por defecto
                    # Contar las ocurrencias de cada condición
                    for ests in infoAf:
                        indicador = 0
                        for condi in conditions:
                            if ests.startswith(condi):
                                counts[condi] += 1
                                indicador = 1
                        if indicador == 0:
                            raise ValueError(f"Tipo no encontrado : {line}\\tline: {count}")
                    # Calcular la longitud y escribir la información en el archivo de salida
                    largo = int(col[2]) - int(col[1])
                    #largo = int(col[2]) - int(col[1]+1)
                    output_line = [col[0], col[1], col[2], col[0] + f"_SmallRNA_{count}", str(largo), col[4], str(len(infoAf))] + \
                                [str(counts[condi]) for condi in conditions]
                    pri.write("\\t".join(output_line) + "\\n")
                    count += 1
    except Exception as e:
        print(e)
    """
}
process ShortCatPythonPerCHR {
    publishDir "${params.outdir}/04-ShortCatPython", mode:"copy", pattern: '*.counts'
    publishDir "${params.outdir}/04-ShortCatPython/BEDs", mode:"copy", pattern: '*.bed'
    maxForks 10
    cpus 6
    memory '11 GB'

    input:
    path(SAM)
    val(Indices)

    output:
    path("all*.counts"), emit: counts
    path("*.bed")

    script:
    """
    #!/usr/bin/env python
    import os
    import subprocess

    # Dividir el nombre del archivo de entrada para crear el nombre del archivo de salida si no se proporciona
    filename = "${SAM}".replace(".reads.sam","")
    fileout = filename + ".counts"

    # Dividir las condiciones por coma
    conditions = "${Indices}".replace("[","").replace("]","").replace(" ","").split(",")

    # Ejecutar el comando sam2bed

    subprocess.run("sam2bed --max-mem 10G < ${SAM} > " + filename + ".bed", shell=True)

    # Ejecutar el comando mergeBed dependiendo de la especificación de hebra
    if "${params.Strand}" == 'yes':
        subprocess.run("mergeBed -s -i " + filename + ".bed -c 4,6 -o distinct -delim \\";\\" > " + filename + "_merge.bed", shell=True)
    else:
        subprocess.run("mergeBed -i " + filename + ".bed -c 4,6 -o distinct -delim \\";\\" > " + filename + "_merge.bed", shell=True)


    # Asumiendo que las condiciones y el nombre del archivo de entrada ya están definidos
    filein = filename + "_merge.bed"

    # Abrir el archivo de entrada
    try:
        with open(filein, 'r') as ace:
            # Inicializar el contador y el indicador
            count = 1
            indicador = 0
            # Crear el archivo de salida si no se ha definido
            fileout = filein.split('.')[0] + ".count" if not fileout else fileout
            # Abrir el archivo de salida para escribir
            with open(fileout, 'w') as pri:
                # Escribir los encabezados en el archivo de salida
                headers = ["Chrom", "Start", "End", "Contig", "Length", "strand", "#reads"] + conditions
                pri.write("\\t".join(headers) + "\\n")
                # Leer cada línea del archivo de entrada
                for line in ace:
                    line = line.strip()
                    col = line.split('\\t')
                    infoAf = col[3].split(';')
                    counts = {condi: 0 for condi in conditions}  # cambiar a 1 para sumarle 1 a todas las cuentas por defecto
                    # Contar las ocurrencias de cada condición
                    for ests in infoAf:
                        indicador = 0
                        for condi in conditions:
                            if ests.startswith(condi):
                                counts[condi] += 1
                                indicador = 1
                        if indicador == 0:
                            raise ValueError(f"Tipo no encontrado : {line}\\tline: {count}")
                    # Calcular la longitud y escribir la información en el archivo de salida
                    largo = int(col[2]) - int(col[1])
                    #largo = int(col[2]) - int(col[1]+1)
                    output_line = [col[0], col[1], col[2], col[0] + f"_SmallRNA_{count}", str(largo), col[4], str(len(infoAf))] + \
                                [str(counts[condi]) for condi in conditions]
                    pri.write("\\t".join(output_line) + "\\n")
                    count += 1
    except Exception as e:
        print(e)
    """
}
process Make_parsed_file { 
    publishDir "${params.outdir}", mode:'copy'

    input:
    path(Counts)

    output:  
    path("All_counts.tsv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    import os 
    from os import walk
    dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2]).rename(columns={0:"file"})
    tab_files = dir_files.loc[dir_files.file.str.contains(".counts")].copy().reset_index()
    All_counts = pd.DataFrame()
    for File in tab_files.file:
        print(File)
        Dataframe = pd.read_csv(File,sep="\\t")
        All_counts = pd.concat([All_counts,Dataframe])
    
    All_counts.to_csv("All_counts.tsv", sep="\\t",header=True, index=False)
    """  
}
process DESeq2 {    
    publishDir "${params.outdir}/05-DEA", mode:"copy"
    cpus 15  

    input:
    path(Matrix)
    path(MetaData)

    output:
    path "*"
    
    script:
    """   
    #!/usr/bin/env Rscript 
    suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(BiocParallel)
    library(stringr)
    library(ggrepel)
    })
    register(MulticoreParam(${task.cpus}))

    Custom_theme <- function(){
    theme_bw() +
        theme(panel.background = element_rect(colour = "black", size=0.1), 
            plot.title = element_text(hjust = 0.5, angle = 0, size = 15, face = "bold", vjust = 1),
            axis.ticks.length=unit(.2, "cm"), axis.text = element_text(size=11), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }


    Deseq2_ShortCatNF <- function(Counts,Sample_data,Min_Reads,Control,Case) {
        expr0 <- read.delim(Counts, sep = "\\t")
        rownames(expr0) <-expr0\$Contig
        expr0 <- expr0[,8:ncol(expr0)]
        keep <- rowSums(expr0>=Min_Reads) >= ncol(expr0)*0.5
        expr0 <- expr0[keep,]
        sample_data <- read.delim(Sample_data, sep = ",")
        sample_data <- sample_data[sample_data\$Condition == Control | sample_data\$Condition == Case,]
        names(sample_data) <- c("column", "Condition")
        coldata <- sample_data 
        matrix <- as.data.frame(expr0)
        RNAs <- "DEA_sRNAs"   
        cts <- matrix[coldata\$column] 
        dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
        keep <- rowSums(counts(dds)) > length(colnames(cts)) ## Pre-Filtering
        dds <- dds[keep,]
        dds\$Condition <- factor(dds\$Condition, levels = c(Control,Case))
        dds\$Condition <- droplevels(dds\$Condition)
        dds <- DESeq(dds, fitType="local", parallel = TRUE)
        DESeq_norm <- counts(dds, normalized=T)
        write.table(DESeq_norm, file = paste0(RNAs,"_",Control,"_vs_",Case,"_DESeq_norm_matrix_WT.csv"), row.names = T, quote = F, col.names = T, sep = ",")  
        vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
        pcaData <- plotPCA(vsd, intgroup="Condition", returnData=TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
        ggtitle(paste(RNAs, sep = "")) +  Custom_theme() +
        theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic"), axis.title.x =element_text(size=12), axis.title.y = element_text(size=12), legend.title = element_text(face = "bold")) +
        geom_point(size=2) +
        geom_text_repel(aes(label=as.character(name)),hjust=0,vjust=0,size=2) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() +
        scale_color_manual(values=c("#56B4E9", "red"))
        ggsave(paste0(RNAs,"_",Control,"_vs_",Case,"_PCA_plot.png"), plot = p, width = 8, height = 8, dpi = 300, units = "in")
        res <- results(dds, parallel = TRUE)
        res2 <- as.data.frame(res)
        res2\$L2FC_ABS <- abs(res2\$log2FoldChange)
        res2 <- res2[order(res2\$L2FC_ABS,decreasing=T),]
        write.csv(na.omit(res2[res2\$padj < 0.05,]), file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE.csv", sep= ""))
        write.csv(res2, file=paste(RNAs,"_",levels(dds\$Condition)[1], "_vs_", levels(dds\$Condition)[2],"_DE_NS.csv", sep= ""))
    }
    Deseq2_ShortCatNF("${Matrix}","${MetaData}",${params.cr},"${params.Control_name}","${params.Case_name}")
    """
}
workflow { // Main workflow to analize RNAseq data
    params.outdir = "ShortCatNFDocker/"
    out_dir = file(params.outdir) // transforms the path of the main directory in a object of the nextflow pipeline
    out_dir.mkdir() // creates the main directory in the system 
    ST= channel.fromFilePairs(params.Samples,size:-1)
    Sample_IDs = ST.map{it[0]}
    if (params.Genome != ""){
        BowtieIndexBuild(Channel.fromPath(params.Genome))
        Index_ch = BowtieIndexBuild.out.index
    } else {
        Index_ch = Channel.fromPath(params.Bowtie_index)
    }
    QC1(ST) // calls fastqc and fastp tools to check quality
    if (params.Adapters_fasta == ""){ 
        Remove_duplicated_adapters1(QC1.out.collect())
        Deduplicated = Remove_duplicated_adapters1.out
    }else {
        Remove_duplicated_adapters2(QC1.out.collect(),channel.fromPath(params.Adapters_fasta))
        Deduplicated = Remove_duplicated_adapters2.out
    }
    Trimm_QC2(ST,Deduplicated) // makes the quality control of the trimming
    Bowtie(Trimm_QC2.out,Index_ch)
    process_SAMs(Bowtie.out.mapped_reads.collect())
    if (params.Perl != false){
        ShortCat(process_SAMs.out.SAM.flatten(),Sample_IDs.toList())
        Make_parsed_file(ShortCat.out.counts.collect())}
    else{
        if (params.PerCHR == false){
            ShortCatPython(process_SAMs.out.SAM.flatten(),Sample_IDs.toList())
            Make_parsed_file(ShortCatPython.out.counts.collect())
        }else{
            ShortCatPythonPerCHR(process_SAMs.out.SAM.flatten(),Sample_IDs.toList())
            Make_parsed_file(ShortCatPythonPerCHR.out.counts.collect())
        }}
             
    if (params.DESeq2 != false){
        DESeq2(Make_parsed_file.out,channel.fromPath(params.ConditionsFile))
    } 
}
workflow QC1{
    take: data
    main:
        FASTQC(data) 
        Get_adapters(FASTQC.out.html)
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC1")
    emit:
        Get_adapters.out
}
workflow Trimm_QC2{
    take: data
          data2
    main:
        FASTP(data,data2)    
        FASTQC(FASTP.out.trimmed_reads)    
        MULTIQC(FASTQC.out.multiqc_input.collect(),"QC2")
    emit:
        FASTP.out.trimmed_reads
}
workflow prueba { // Main workflow to analize RNAseq data
    
}
//perl 1_automatize_shortcat_DROPBOX.pl -t 60 -x Index/GRCh38.p13 --rm 0:0:0 --ANNOT yes -1 /media/storage2/Adolfo2/LIB_works/ShortCat/RawData/hsa.gff3 -i Samples.txt
// nextflow Scripts/ShortCat.nf -with-conda -resume --sample_table=Samples.Seba_Urquiza.txt --mod_FastP=-g --Adapters_fasta=/media/storage1/Adolfo/LIB_works/ShortCat/bin/z_adapters_truseq.fasta
// nextflow GitHub/ShortCat/ShortCat.Docker.nf --Samples="/media/storage2/surquiza_docs/datasets_rat_HC/dataset_NE_cardios/miRNA/*_mi.fq.gz" --mod_FastP='-g' --Adapters_fasta='/media/storage1/Adolfo/LIB_works/ShortCat/bin/z_adapters_truseq.fasta' -c GitHub/ShortCat/nextflow.config  -resume
// '/media/storage1/Adolfo/LIB_works/ShortCat/Scripts/Conditions.csv'