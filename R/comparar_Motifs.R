#' @title Compara los motif (lista) entre ellos
#' @description Función que nos ayuda a observar cuáles son los motifs que se
#' encuentran más próximos entre ellos a nivel de secuencia.
#' @param motif Lista de motif
#' @return Matriz de puntuación del alineamiento de los motifs
#' @note
#'  comparar_Motifs(motif)
#'
#' @export
comparar_Motifs=function(motif){
  compare_motifs(motif, use.type = "PPM", method = "PCC", tryRC = F)
}
