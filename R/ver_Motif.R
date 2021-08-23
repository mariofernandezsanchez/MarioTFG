#' @title Dibuja un plot con el motif indicado
#' @description Función que dibuja un gráfico del motif que se le índica como parametro
#' @param x Motif que queremos ver
#' @return El gráfico con el logo de secuencia del motif
#' @note
#' motif <- leer_Motif("./xxx/xxx.motif")
#'
#' ver_Motif(motif[1])
#' @export
ver_Motif=function(x){
  view_motifs(motifs = x)
}
