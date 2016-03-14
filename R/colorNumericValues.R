#' Transform intensity values in a gradient of color from green to red.
#' @author Simon J Pelletier
#' @param values Intensity values to be transform in color gradients
#' @return A vector of colors. The darkest green values are the lowest and the darkest red values are the highest.
#' @keywords comparisons
#' @export
colorNumericValues <- function(values){
  if(min(values) < 0){                 #If no value is lower than 0, the values are unchanged
    values2 = values + abs(min(values)) #Takes all the values to a minimum of 0 by adding the absolute of the lowest negative number
  }
  values3 = values2/max(values2) #Transformation of the values in relative values between 0 - 1
  f <- colorRamp(c("green","white", "red"))
  colors <- rgb(f(values3)/255)
  return(colors)
}
