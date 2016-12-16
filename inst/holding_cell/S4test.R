setClass("opcbinomial",
         slots = c(grid = "GridTopology",
                   server = ""),
         validity = function(object) {
           return(TRUE)
         }
)
