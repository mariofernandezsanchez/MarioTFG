#' @title Convert a motif list in format ".fasta" to DNAStringSet
#' @description Function that helps us to convert a background sequence
#' in a fasta format to one in DNAStringSet format
#' @param ourpath Path where the sequence will be found
#' @param nFlist The file .fasta (number in your directory)
#' @return The DNAStringSet file
#' @note
#' region <- "peaks" (you need load this value to make the function work)
#' file<-convertToDNAStringSet=function("ourpath", nFlist)
#' Example: seq1<-convertToDNAStringSet=function(".", 1)
#' seq1<-convertToDNAStringSet=function(".", 2)
#'
#' @export
convertToDNAStringSet=function(ourpath, nFlist){

  list_files <- base::list.files(path = ourpath , pattern = region)
  file_seq <- (utils::read.table(list_files[nFlist]))
  file_seq <- file_seq[!BiocGenerics::grepl("^>", file_seq$V1),]
  file_seq <- BiocGenerics::as.vector(file_seq)
  seq_file2 <- Biostrings::DNAStringSet(file_seq)
}

#' @title Compare motifs (list)
#' @description Function that compares all motifs to all other motifs and
#' helps us to observe which are the motifs that are are closer to each other
#' at the sequence level.
#' @param motif Motif in formar list
#' @param type The type of matrix that represents the motif and that will be
#' used for alignment. You can chose between "PPM" or "ICM"
#' @param method Alignment metric between motifs. You can chose one method: PCC,
#' EUCL,SW,KL,ALLR,BHAT,HELL,SEUCL,MAN,ALLR_LL, WEUCL, WPCC. See details
#' @return Motif alignment scoring matrix
#' @note
#' compare_Motifs(motif, "type", "method") Dont't forget use the quotes
#'
#' @export
compare_Motifs=function(motif, type, method){
  universalmotif::compare_motifs(motif, use.type = type, method = method,
  tryRC = F)
}

#' @title Create motifs
#' @description Function that ...
#' @param motifs The list of motif
#' @param name Name of motif.
#' @return A new motif (class universalmotif)
#' @note
#'
#' scan50refined <- scan_Motif (motif,seq, th)
#' 50refined <- create_Motif(scan50refined, "motifp30_50refined")
#'
#' @export
create_MotifRefined=function(scan_target,name){
  universalmotif::create_motif(Biostrings::DNAStringSet(scan_target$match),
                               name = name, pseudocount = 1/3)
}

#' @title Enrich motifs
#' @description Function that helps us to make an enrichment analysis
#' of the refined motif compared to its background sequence.
#' @param motifs The motif to analyze
#' @param seq The target sequence
#' @param bkg The background sequence
#' @param th Threshold that determines whether or not a sequence contains the
#' motif
#' @return The enrichment analysis of this motif
#' @note
#'
#' enrich_Motifs=function(cdk2_100,seq_file_cdk2_100,0.7)
#' @export
enrich_Motifs=function(motifRefined,seq,bkg,th){
  universalmotif::enrich_motifs(motifs = motifRefined, sequences = seq,
                                bkg.sequences = bkg, RC = T, threshold = th,
                                threshold.type = "logodds", qval.method = "BH")
}

#' @title Export motifs to jaspar file
#' @description Function that we use to export the motif in Jaspar format
#' @param motif Motif already merged
#' @param filename Name of the file
#' @return The definitive motif
#' @note
#' export_Jaspar(motifxxx,"file.txt") *don´t forget the ""
#'
#' @export
export_Jaspar=function(motif, filename){
  universalmotif::write_jaspar(motif,filename)
}
#' @title Export motifs to universalmotif file
#' @description Function we use to export the motif in universalmotif format
#' @param motif Motif already merged
#' @param filename Name of the file
#' @return The definitive motif
#' @note
#' export_Jaspar(motifxxx,"file.txt") *don´t forget the ""
#'
#' @export
export_Motif=function(motif, filename){
  universalmotif::write_motifs(motif,filename)
}

#' @title Convert motif (list) in a Data Frame
#' @description Function that helps us convert a motif in list format to
#  data.frame format.
#' @param motif Motif we want to see
#' @return Motif in format data.frame
#' @note
#' listToDF_Motifs(motif)
#'
#' motifDF <- listToDF_Motifs(motif)
#' @export
listToDF_Motifs=function(motif){
  universalmotif::summarise_motifs(motif)
}

#' @title Merge motifs
#' @description Function that helps us to merge all the motifs into one,
#' to obtain the final motif.
#' @param motifs List of motifs we want to merge
#' @return The definitive motif
#' @note
#'
#' merge.motifs <- merge_Motifs=function(motif)
#'
#' @export
merge_Motifs=function(motifs){
  universalmotif::merge_motifs(motifs, method = "PCC", use.type = "PPM",
                               tryRC = F)
}

#' @title Import motifs in format HOMER
#' @description Function that loading y makes a motif understandable to humans
#' @param ourpath Path in our pc of the motif that we want to import
#' @return The motif list with the parameters more importants
#' @note
#' read_Motifs(ourpath = "./example/example.motif")
#'
#' @export
read_Motifs=function(ourpath){
  universalmotif::read_homer(ourpath)
}

#' @title Get a heatmap of the selected motifs (list)
#' @description Function that by means of a graphical plot allows us
#' to visualize the matrix obtained in the comparison of motifs, and that
#' helps us to select the cluster that we are looking for.
#' @param matrix Alignment scoring matrix
#' @return Heatmap of the clouter motifs
#' @note
#'  heatMap_Motifs(compare_Motifs(motif = xxx,"xxx","xxx"))
#'
#' @export
heatMap_Motifs=function(matrix){
  heatmaply::ggheatmap(matrix)
}

#' @title Score of all motifs in a sequence (data.frame)
#' @description Function that obtains the punctuation of all the motifs of a
#' sequence passed as data.frame according to the mathematical formula of
#' (Almagro-Hernández and Fernández-Breis, 2020)
#' @param motif The motif from which you want to obtain the scores
#' @return The score result for that motif / s
#' @note
#' m <- read_Motifs(rutadelmotiv = "./ejemplo/ejemplo.motif")
#'
#' motifDF <- listToDF_Motifs(m)
#'
#' scoreSeq_Motifs_Motifs(motifDF)
#'
#' @export
scoreSeq_Motifs=function(motif){
  for (i in 1:length(motif$name)) {
    print(-log(motif$pval[i])*log(motif$nsites[i]))
  }
}


#' @title Score of motif
#' @description Function that gets the punctuation of a motif specified by
#' parameter according to the mathematical formula of (Almagro-Hernández and
#' Fernández-Breis, 2020)
#' @param motif The motif in data.frame format from which you want to obtain the
#' scores
#' @param x Motif number within the sequence
#' @return The score result for that motif
#' @note
#' score_Motif(motifdf, x)
#'
#' @export
score_Motif=function(motif, x){
  -log(motif$pval[x])*log(motif$nsites[x])
}

#' @title Scan sequences for matches to input motifs
#' @description This function is used to scan the motif of each set
#' in its respective set of peaks sequences, in order to obtain the subsequences
#' of the same length as the motif, that exceed the established threshold.
#' If that score exceeds the predetermined threshold, it is considered that
#' there is a match and that subsequence is one of those that the function
#' returns and will be used to generate the refined motif.
#' @param ClusterMotif The motif in correct format
#' @param seq Sequence obtained in step 1 of each experiment
#' in DNAStringSet format
#' @param th Threshold that determines whether or not a sequence contains the
#' motif
#' @return The analysis of scanned sequences
#' @note
#' scan_Motif(ClusterMotif,seq_file,0,70)
#'
#' @export
scan_Motif=function(ClusterMotif,seq,th){
  universalmotif::scan_sequences(motifs = ClusterMotif,sequences = seq,
                                 threshold = th,threshold.type = "logodds",
                                 RC = TRUE)
}

#' @title Draw a plot with the indicated motif
#' @description Function that draws a graphic of the motif that is indicated as
#' a parameter
#' @param motif Motif we want to see
#' @param n Number of the motif within the sequence
#' @return The graphic with the motif sequence logo
#' @note
#' motif <- read_Motifs("./xxx/xxx.motif")
#'
#' see_Motif(motif)
#' @export
see_Motif=function(motif,n){
  universalmotif::view_motifs(motif[n])
}


#' @title Trim the motifs
#' @description Function that eliminates the positions of those motifs that
#' contain little information
#' @param motif Motif already merged
#' @return The definitive motif
#' @note
#' def.motif <- trim_Motifs(motif)
#'
#' @export
trim_Motifs=function(motif){
  universalmotif::trim_motifs(motifs = motif, min.ic = 0.25)
}

#' @title Etapa 5
#' @description Function that we use to obtain the refined motif
#' @param motifcl The motif obteined in the cluster
#' @param seq_file The target sequences
#' @param seqbkg The background sequences
#' @param name Name of the refined motif
#' @return The refined motif
#' @note
#'
#' example<-etapa5_Motifs(clusterp30[1],seqp30_50,SeqBackground,"Refined_p50")
#' TODO: Añadir extra info y ¿bkg?
#' @export
etapa5_Motifs=function(motifcl,seq_file,seqbkg,name){

    scan<-scan_Motif(motifcl,seq_file,0.7)
    mrefined<-create_MotifRefined(scan, name)
    menrich<-enrich_Motifs(mrefined,seq_file,seqbkg,0.7)
    mrefined@pval <- menrich@listData[["Pval"]]
    mrefined@eval <- menrich@listData[["Eval"]]
    mrefined@qval <- menrich@listData[["Qval"]]
    mrefined
}

#' @title Etapa 1
#' @description Function that returns the genomic sequence of the specified
#' region of the given species.
#' @param species The species name/alias. (EX: human, homo_sapiens, mouse)
#' @param start The start of query region. A maximum of 10Mb is allowed to be
#' requested at any one time
#' @param end The end of query region.
#' @return The sequence in format fasta
#' @note
#' The value of the species argument must be loaded.
#' human <- "human", mouse="mouse" (you need load this value to
#' make the function work)
#' etapa1_Motifsetapa1_Motifs(human,16185,16235)
#' @export
etapa1_Motifs=function(species,start,end) {

  library(httr)
  library(jsonlite)
  library(xml2)
  if (start >= end){
    stop("The start can't be bigger than end")
  }

  s<-format(start, scientific = F) #evitar numeros tipo 1e+6
  e<-format(end, scientific = F)
  server <- "https://rest.ensembl.org"
  #ext <- paste("/sequence/region/human/X:1000000..1000050:?")
  ext2 <- paste0("/sequence/region/",species,"/X:",s,"..",e,":?", collapse = ",")

  r <- GET(paste(server, ext2, sep = ""), content_type("text/x-plain"))

  stop_for_status(r)
  rr <- (content(r))
  res <- cat(paste(rr$id,rr$seq,sep = "\n"));
  return(res)
}
