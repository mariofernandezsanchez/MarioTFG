#' @title Convert a motif list in format ".fasta" to DNAStringSet
#' @description Function that helps us convert ...
#'
#' @param motif Motif we want to see
#' @return Motif in format data.frame
#' @note
#' listToDF_Motifs(motif)
#' @export
convertToDNAStringSet=function(ourpath, nFlist){

  list_files <- base::list.files(path = ourpath , pattern = region)
  file_seq <- (utils::read.table(list_files[nFlist]))
  file_seq <- file_seq[!BiocGenerics::grepl("^>", file_seq$V1),]
  file_seq <- BiocGenerics::as.vector(file_seq)
  seq_file <- Biostrings::DNAStringSet(file_seq)
}

#' @title Compare motifs (list)
#' @description Function that compares all motifs to all other motifs and
#' helps us to observe which are the motifs that are are closer to each other at the sequence level.
#' @param motif Motif in formar list
#' @param type The type of matrix that represents the motif and that will be
#' used for alignment. You can chose between "PPM" or "ICM"
#' @param method Alignment metric between motifs. You can chose one method: PCC,
#' EUCL,SW,KL,ALLR,BHAT,HELL,SEUCL,MAN,ALLR_LL, WEUCL, WPCC. See details
#' @return Motif alignment scoring matrix
#' @note
#'  compare_Motifs(motif, "type", "method") Dont't forget use the quotes in type and method
#'
#' @export
compare_Motifs=function(motif, type, method){
  universalmotif::compare_motifs(motif, use.type = type, method = method, tryRC = F)
}

#' @title Create motifs
#' @description Function that ...
#' @param consensus The consensus of the motif
#' @param name Name of motif
#' @param pseudoc
#'
#' @return A new motif (class universalmotif)
#' @note
#'  alphabet default DNA, type default PPM
#'
#' @export
create_Motif=function(consensus,name){
  universalmotif::create_motif(input = consensus,name = name,pseudocount = 1/3)
}

#' @title Create motifs
#' @description Function that ...
#' @param consensus The consensus of the motif
#' @param name Name of motif
#' @param pseudoc
#'
#' @return A new motif (class universalmotif)
#' @note
#'  alphabet default DNA, type default PPM
#'
#' @export
create_MotifBkg=function(listbkg){
  universalmotif::create_motif(bkg = listbkg)
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


#' @title Import motifs in format HOMER
#' @description Function that loading y makes a motif understandable to humans
#' @param rutadelmotiv Path in our pc of the motif that we want to import
#' @return The motif list with the parameters more importants
#' @note
#' read_Motifs(rutadelmotiv = "./example/example.motif")
#'
#' @export
read_Motifs=function(rutadelmotiv){
  universalmotif::read_homer(rutadelmotiv)
}

#' @title Get a heatmap of the selected motifs (list)
#' @description Function that by means of a graphical plot allows us to visualize the
# 'matrix obtained in the comparison of motifs, and that helps us to select the
# 'cluster that we are looking for.
#' @param motif List of motifs
#' @return Motif alignment scoring matrix
#' @note
#'  heatMap_Motifs(motif)
#'
#' @export
heatMap_Motifs=function(motif){
  heatmaply::ggheatmap(motif)
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
#' @param motif The motif in data.frame format from which you want to obtain the scores
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
#' @description TODO
#' @param motif The motif in data.frame format from which you want to obtain the scores
#' @param x Motif number within the sequence
#' @return The score result for that motif
#' @note
#' score_Motif(motifdf, x)
#'
#' @export
scan_Motif=function(motif,seq,th){
  universalmotif::scan_sequences(motifs = motif,sequences = seq,threshold = th,
                                 threshold.type = "logodds",RC = TRUE)
}

#' @title Draw a plot with the indicated motif
#' @description Function that draws a graphic of the motif that is indicated as a parameter
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
