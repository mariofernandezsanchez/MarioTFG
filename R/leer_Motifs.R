#' @title Importa motifs en formato HOMER
#' @description Función que carga y hace comprensible para el humano un motif en formato HOMER
#' @param rutadelmotiv Ruta en nuestro pc del motif que queremos importar
#' @return El motif importado con los parametros más destacables
#' @note
#' leer_Motifs(rutadelmotiv = "./ejemplo/ejemplo.motif")
#'
#' @export
leer_Motifs=function(rutadelmotiv){
  read_homer(rutadelmotiv)
}
