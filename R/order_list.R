#' Order the element of a list by size in decreasing order
#' @param list List to order alphabetically
#' @return The list ordered
#' @examples
#' list <- list(a="toto",b=c("tyro","totor","cipaille","maple_syrup"),c=c("coconut","rdl"))
#' order_list(list)
#' @export
order_list = function(list){
  list = list[order(sapply(list,length),decreasing=T)]
  return(list)
}
