#' @title Puntuación del motif
#' @description Función que obtiene la puntuacióin de un motif especificado por
#' parametro según la formula matematica de (Almagro-Hernández and
#' Fernández-Breis, 2020)
#' @param motif El motiv del que se quiere obtener las puntuaciones
#' #' @param x Numero del motif dentro de la secuencia
#' @return El resultado de la puntuación para ese motif
#' @note
#' puntuacion_Motif(motif)
#'
#' @export
puntuacion_Motif=function(motif, x){
  -log_string_pval(motif$pval[x])*log_string_pval(motif$nsites[x])
}
