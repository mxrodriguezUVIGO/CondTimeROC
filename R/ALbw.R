ALbw <-
function(type_kernel = "n", vec_data)
#######################################################################
#               REQUIRED INPUTS FOR THE FUNCTION        		          #
#######################################################################
#   "type_kernel" kernel function:  "e" Epanechnikov,	"n" Normal, 
#   "vec_data" sample data to estimate distribution function
{
  n <- length(vec_data)
  orderr <- 0
  # pilot badwidth
  bp <- (n^(-0.3)) * sd(vec_data)
 
  # weight function
  ss <- quantile(vec_data, c(orderr, 1-orderr))
  w <- (ss[1]<=vec_data) & (vec_data<=ss[2])

  # plug-in estimator of the asymptotically optimal bandwidth
  aux <- outer(vec_data, vec_data, "-")/bp
  aux <- kernel_function(type_kernel, aux)
  diag(aux) <- 0
  D2_F <- sum(aux*w)/(bp*n*(n-1))
  V2 <- 2 * A1_k(type_kernel) * D2_F
  aux1 <- outer(vec_data, vec_data, "-")/bp
  aux1 <- derivative_kernel_function(type_kernel, aux1)
  aux2 <- apply(t(aux1)%*%(aux1*w), 1, sum)
  D3_F <- sum(aux2)/((n^3)*(bp^4))
  B3 <- 0.25 * (A2_k(type_kernel))^2 * D3_F

  ALbw <- (((0.25*V2)/B3)^(1/3)) * (n^(-1/3))
  
  return(ALbw)

}
kernel_function <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,  "n" Normal, 
#                                  "b" Biweight, "t" Triweight            
#   "u" array or single value where the kernel is evaluated

{
  if(type_kernel == "e") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- 0.75 * (1 - (Uval^2))
    return(result)
  } else if(type_kernel == "n") {
    result <- dnorm(u)
    return(result)
  } else if(type_kernel == "b") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- (15/16) * ((1 - (Uval^2)))^2
    return(result)
  } else if(type_kernel == "t") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- (35/32) * ((1 - (Uval^2)))^3
    return(result)
  }
}

derivative_kernel_function <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,  "n" Normal, 
#                                  "b" Biweight, "t" Triweight         
#   "u" array or single value where the kernel is evaluated

{
  if(type_kernel == "e") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- -1.5 * Uval
    return(result)
  } else if(type_kernel == "n") {
    result <- (1/(sqrt(2 * pi))) * exp(-(u^2)/2) * (-u)
    return(result)
  } else if(type_kernel == "b") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- -(15/4) * Uval * (1 - (Uval^2))
    return(result)
  } else if(type_kernel == "t") {
    result <- u
    logic0 <- (u <= -1 | u >= 1)
    logic1 <- (u > -1 & u < 1)
    result[logic0] <- 0
    Uval <- result[logic1]
    result[logic1] <- -(105/16) * Uval * ((1 - (Uval^2)))^2
    return(result)
  }
}
A1_k <-
function(type_kernel)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,  "n" Normal, 
#                                  "b" Biweight, "t" Triweight         

{
  if(type_kernel == "e")  
    result <- 0.12857 
  else 
    if(type_kernel == "n")  
      result <- 0.28209    
  else 
    if(type_kernel == "b")  
      result <- 0.10823  
  else 
    if(type_kernel == "t")  
      result <- 0.095183  
  return(result)
}

A2_k <-
function(type_kernel)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,  "n" Normal, 
#                                  "b" Biweight, "t" Triweight         

{
  if(type_kernel == "e")  
    result <- 1/5 
  else 
    if(type_kernel == "n")  
      result <- 1    
  else 
    if(type_kernel == "b")  
      result <- 1/7  
  else 
    if(type_kernel == "t")  
      result <- 1/9 
  return(result)
}