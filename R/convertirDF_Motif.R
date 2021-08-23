#' @title Convierte un motif (lista) en un Data Frame
#' @description Función que nos ayuda a convertir un motif en formato lista a
#' formato data.frame. Facilitando su visibilidad y manejo en posteriores
#' funciones
#' @param motif Motif que queremos ver
#' @return Motif en formato data.frame
#' @note
#'  convertirDF_Motif(motif)
#'
#'  motifDF <- convertirDF_Motif(motif)
#' @export
convertirDF_Motif=function(motif){
 summarise_motifs(motif)
}
