# ShortCat

ShortCat es un flujo de trabajo de Nextflow diseñado para el análisis de secuencias de RNAseq. Este flujo de trabajo incluye varios procesos que se ejecutan en secuencia y en paralelo para optimizar el tiempo de análisis.

## Procesos

1. **FASTQC**: Este proceso realiza un control de calidad inicial en los archivos de secuencias de entrada utilizando la herramienta FastQC.

2. **MULTIQC**: Este proceso recopila los informes de FastQC y genera un informe consolidado utilizando la herramienta MultiQC.

3. **Get_adapters**: Este proceso extrae las secuencias de adaptadores sobre-representadas de los informes de FastQC.

4. **Concat_adapters**: Este proceso concatena todos los archivos de adaptadores en un solo archivo.

5. **Remove_duplicated_adapters1 y Remove_duplicated_adapters2**: Estos procesos eliminan las secuencias de adaptadores duplicadas del archivo de adaptadores.

6. **FASTP**: Este proceso realiza el recorte de las secuencias de entrada utilizando la herramienta FastP.

7. **Bowtie**: Este proceso realiza el mapeo de las secuencias recortadas contra un genoma de referencia utilizando la herramienta Bowtie2.

8. **BowtieIndexBuild**: Este proceso construye el índice del genoma de referencia para el mapeo con Bowtie2.

## Parámetros

El flujo de trabajo ShortCat acepta varios parámetros que se pueden especificar en el archivo de configuración `nextflow.config` o en la línea de comandos al ejecutar el flujo de trabajo. Algunos de estos parámetros incluyen:

- `Samples`: Define el archivo de tabla de muestra de entrada.
- `mod_FastP`: Indica si es necesario un recorte de poli G u otra bandera FastP.
- `Organism`: Define el organismo para analizar.
- `Library_Layout`: Indica si la biblioteca es de extremo único o de extremo pareado.
- `Control_name`: Define el nombre del control.
- `Case_name`: Define el nombre del caso.
- `short_reads`: Establece TRUE si el RNAseq tiene lecturas de menos de 60 - 50 nt.
- `Adapters_fasta`: Define el archivo fasta de adaptadores.
- `Genome`: Define el genoma de referencia.
- `Bowtie_index`: Define el índice de Bowtie2.
- `Index_prefix`: Define el prefijo del índice.
- `Conditions`: Define las condiciones del experimento.

## Ejecución

Para ejecutar el flujo de trabajo ShortCat, se puede utilizar uno de los siguientes comandos:

### Usando archivo YAML
~~~
nextflow GitHub/ShortCat/main.nf -c GitHub/ShortCat/nextflow.config -params-file GitHub/ShortCat/ParamsConfig.ejemplo.yml -resume 
~~~

### Detallando parametro en terminal
~~~
nextflow GitHub/ShortCat/main.nf --Samples="/media/storage2/surquiza_docs/datasets_rat_HC/dataset_NE_cardios/miRNA/*_mi.fq.gz" --mod_FastP='-g' --Adapters_fasta='/media/storage1/Adolfo/LIB_works/ShortCat/bin/z_adapters_truseq.fasta' -c GitHub/ShortCat/nextflow.config  -resume --ConditionsFile="/media/storage1/Adolfo/LIB_works/GitHub/ShortCat/Test/Conditions.csv" --Control_name=Ctrl6 --Case_name=Nora24 --Organism="Rattus_norvegicus" --Genome="/media/storage1/Adolfo/LIB_works/Primary_assembly_ensembl_hard.fa"
~~~