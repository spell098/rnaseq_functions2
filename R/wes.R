#' Simple function to get color palettes from the wesanderson package,
#' which contain palettes from wes anderson movies
#' @import wesanderson
#' @export
wes <- function(){
  col1 = wes_palette("Royal1")
  col2 = wes_palette("Royal2")
  col3 = wes_palette("GrandBudapest")
  col4 = wes_palette("Moonrise1")
  col5 = wes_palette("Moonrise2")
  col6 = wes_palette("Moonrise3")
  col7 = wes_palette("Cavalcanti")
  col8 = wes_palette("Chevalier")
  col9 = wes_palette("BottleRocket")
  col10 = wes_palette("Darjeeling")
  col11 = wes_palette("Darjeeling2")
  col = c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11)
  return(col)
}
