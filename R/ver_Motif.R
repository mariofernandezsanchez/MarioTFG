#' @title Dibuja un plot con el motif indicado
#' @description Función que dibuja un gráfico del motif que se le índica como parametro
#' @param motif Motif que queremos ver
#' @param n Numero del motif dentro de la secuencia
#' @return El gráfico con el logo de secuencia del motif
#' @note
#' motif <- leer_Motif("./xxx/xxx.motif")
#'
#' ver_Motif(motif[n])
#' @export
ver_Motif=function(motif,n){
  view_motifs(motifs = motif[n])
}
