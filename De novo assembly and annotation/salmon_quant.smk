SAMPLES = {'Control_1','Control_2','Control_3',
           'Infected_1','Infected_2','Infected_3',
           'Kana_1','Kana_2','Kana_3',
           'Kana_infected_1','Kana_infected_2','Kana_infected_3'}

index = 'gene_expression/termites_index'
out_dir = 'gene_expression'

rule all:
  input:
    expand("clean_data/{sample}_1.clean.fq.gz",sample=SAMPLES),
    expand("clean_data/{sample}_2.clean.fq.gz",sample=SAMPLES),
    expand("gene_expression/{sample}/quant.sf",sample=SAMPLES)
    
rule salmon_quant:
  input:
    clean_R1 = "clean_data/{sample}_1.clean.fq.gz",
    clean_R2 = "clean_data/{sample}_2.clean.fq.gz"
  output:
    'gene_expression/{sample}/quant.sf'
  threads: 40
  shell:
    "./salmon-1.9.0_linux_x86_64/bin/salmon quant --gcBias --recoverOrphans --validateMappings --minScoreFraction 0.5 -i {index} -l A "
    "-1 {input.clean_R1} -2 {input.clean_R2} --softclip --threads {threads} -o {out_dir}/{wildcards.sample}"
    

  
