# not in function
'%nin%' <- Negate('%in%')
# is nothing
isnothing = function(x) {
  # test against NULL NA NaN
  any(is.null(x))  || any(is.na(x))  || any(is.nan(x)) 
}

## left hand side vs right hand side
'%||%' <- function(lhs, rhs) {
  # if both true
  if ( !isnothing(lhs) & !isnothing(rhs)){ print("Both sides are free from NULL, NA and NaN") }
  # else
  else{  
    if ( !isnothing(lhs) ) {  lhs } 
    else if ( !isnothing(rhs) ){  rhs }
    else { print("Both sides contain NULL, NA or NaN!")}
    }
}
