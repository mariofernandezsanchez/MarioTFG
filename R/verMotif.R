#' @title Dibuja un plot con el motif indicado
#' @description Función que dibuja un gráfico del motif que se le índica como parametro
#' @param x Motif que queremos ver
#' @return el gráfico con el logo de secuencia del motif
#  Ejemplos
#  motif <- leerMotif(x)
#  verMotif(motif[1])
#' @export
verMotif=function(x){
  view_motifs(motifs = x)
}
