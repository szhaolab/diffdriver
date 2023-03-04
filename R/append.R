#' @param x  new data
#' @param path to old file
#' @return  
append.Rda <- function(x, file) {
     old.objects <- load(file)
     dd=eval(parse(text=old.objects[1]))
     dd[[length(dd)+1]]=x
     save(dd, file = file)
   } 
