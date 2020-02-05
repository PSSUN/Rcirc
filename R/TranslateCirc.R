##' @export
##' @author Peisen Sun

TranslateCirc <- function(out,
                          fastq,
                          adapter,
                          trimmomatic,
                          genome,
                          tmp,
                          circgenome,
                          rRNA) {
  out <- out # output bam
  fastq <- fastq
  adapter <- adapter
  trimmomatic <- trimmomatic
  genome <- geneome
  circgenome <- circgenome
  tmp_file_path <- paste(tmp, "/", sep = "")

  # build genome index
  gindex_command <- paste("bowtie-build --threads 4", genome, tmp_file_path)
  system(gindex_command)

  # build circRNA genome index
  index_command <- paste(
    "STAR --runThreadN 4 --runMode genomeGenerate --genomeDir",
    tmp_file_path,
    paste("--genomeFastaFiles", circgenome)
  )
  system(index_command)

  # build rRNA index
  rindex_command <- paste("bowtie-build --threads 4", rRNA, tmp_file_path)
  system(rindex_command)

  # QC with Trimmomatic for Ribo-seq data
  adapter <- paste("ILLUMINACLIP:", adapter, ":2:30:10", sep = "")
  qc_command <- paste(
    "java -jar",
    trimmomatic,
    "SE -phred33",
    fastq,
    out,
    adapter,
    "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16"
  )
  system(qc_command)

  # remove rRNA
  remove_command <- paste(
    "bowtie -p 4 -norc --un",
    paste(tmp_file_path, "clean.without.rRNA.fastq", sep = ""),
    paste(tmp_file_path, "useless_1", sep = ""),
    paste(tmp_file_path, "clean.fastq", sep = ""),
    ">",
    "map_to_rRNA.sam"
  )
  system(remove_command)

  # remove liner sequence
  remove_liner_command <- paste(
    "bowtie -p 4 -norc --un",
    paste(tmp_file_path, "final_clean.fastq", sep = ""),
    paste(tmp_file_path, "useless_2", sep = ""),
    paste(tmp_file_path, "clean.without.rRNA.fastq", sep = ""),
    ">",
    "map_to_genome.sam"
  )
  system(remove_liner_command)

  # map to circRNA genome
  map_command <- paste(
    "STAR  --runThreadN 4 --outSAMtype BAM SortedByCoordinate --alignIntronMax 10",
    paste("--genomeDir", tmp),
    paste(
      "--readFilesIn",
      paste(tmp_file_path, "final_clean.fastq", sep = "")
    ),
    paste("--outFileNamePrefix", out)
  )
  system(map_command)

  # remove tmp file
  system("rm -rf", tmp)

}
