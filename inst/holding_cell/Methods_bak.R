setGeneric("mytest",
           function(x){
             standardGeneric("mytest")
           })

setMethod("mytest", signature(x="ANY"), function(x){
  setClass("numbers", representation(a = "numeric", b = "numeric"))
  num1 = new("numbers", a = 12, b = 42)
  return(num1)
})

setGeneric("pprint",
           function(x){
             standardGeneric("pprint")
           })


setGeneric("aprint",
           function(x){
             standardGeneric("aprint")
           })

setMethod("aprint", signature(x = "numbers"), function(x){
  print(x@a)
})

setGeneric("bprint",
           function(x){
             standardGeneric("bprint")
           })

setMethod("bprint", signature(x = "numbers"), function(x){
  print(x@b)
})

setMethod("print", signature(x = "numbers"), function(x){
  print(x@b+x@a)
})
