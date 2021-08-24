#' @title Puntuación de todos los motif de una secuencia (data.frame)
#' @description Función que obtiene la puntuacióin de todos los motifs de una
#' secuencia pasada como data.frame según la formula matematica de
#' (Almagro-Hernández and Fernández-Breis, 2020)
#' @param motif El motiv del que se quiere obtener las puntuaciones
#' @return El resultado de la puntuación para ese/os motif/s
#' @note
#' m <- leer_motif(rutadelmotiv = "./ejemplo/ejemplo.motif")
#'
#' motifDF <- convertirDF_Motifs(m)
#'
#' puntuacionSeq_Motifs(motifDF)
#'
#' @export
puntuacionSeq_Motifs=function(motif){
  for (i in 1:length(motif$name)) {
    print(-log_string_pval(motif$pval[i])*log_string_pval(motif$nsites[i]))
  }
}
