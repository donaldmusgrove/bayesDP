#' data_helper
#'
#' data_helper
#'
#' @title data_helper: data_helper
#' @param current dataframe
#' @param historical dataframe
#'
#' @examples
#'
#' @rdname data_helper
#' @export data_helper


setGeneric("data_helper",
           function(current,historical){
             standardGeneric("data_helper")
           })

setMethod("data_helper",
          signature(current = "data.frame", historical = "missing"),
          function(current, historical){

            ok <- current[complete.cases(current),]

            n <- nrow(ok)
            mu <- sapply(ok, mean)
            #sigma <- sapply(ok, sd)
            sigma2 <- sapply(ok, var)
            return(list(n,mu,sigma2))

          })

setMethod("data_helper",
          signature(current = "data.frame", historical = "data.frame"),
          function(current, historical){
            ok <- current[complete.cases(current),]

            n <- nrow(ok)
            mu <- sapply(ok, mean)
            #sigma <- sapply(ok, sd)
            sigma2 <- sapply(ok, var)

            ok0 <- current[complete.cases(historical),]

            n0 <- nrow(ok0)
            mu0 <- sapply(ok0, mean)
            #sigma <- sapply(ok, sd)
            sigma02 <- sapply(ok0, var)

            return(list(n,mu,sigma2,n0,mu0,sigma02))
          })


