#' @title Obtiene un mapa de calor de los motifs seleccionados (lista)
#' @description Función que mediante un plot gráfico nos permite visualizar la
#' matriz obtenida en la comparacion de motifs, y que nos ayuda a seleccionar el
#' cluster que vamos buscando.
#' @param motif Lista de motif
#' @return Matriz de puntuación del alineamiento de los motifs
#' @note
#'  comparar_Motifs(motif)
#'
#' @export
mapaCalor_Motifs=function(motif){
  ggheatmap(comparar_Motifs(motif))
}
