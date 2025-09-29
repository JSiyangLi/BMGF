library(caracas)
library(Rcpp)
sourceCpp("logsum.cpp")

log_terms <- function(expr, evaluation_list) {
  # symbolically expand log(expression)
  
  #---------------------------
  # split and log-expand terms
  #---------------------------
  if (!inherits(expr, "caracas_symbol")) {
    stop("expr must be in the class of caracas_symbol")
  }
  
  # Split into additive terms
  terms <- py_to_r(expr$pyobj$as_ordered_terms())
  
  # Apply log to each term and expand
  terms_list <- lapply(terms, function(t) {
    caracas::expand_log(log(as_sym(t)))
  })
  
  #---------------------------
  # evaluate SymPy expressions
  #---------------------------
  term_sgn <- term_value <- rep(NA, length(terms_list))
  for (i in 1:length(terms_list)) {
    current_term <- terms_list[[i]]
    if (grepl("1i*", as_character(current_term))) { # this term is negative
      # label this term as negative
      term_sgn[i] <- -1
      
      # evaluate this term
      evaluation_list <- c(list(pi = pi, zInf = -Inf), evaluation_list)
      expr_eval <- subs(current_term, evaluation_list)
      cat("expr_eval: ", as_character(expr_eval), "\n")
      term_value[i] <- Re(as_func(expr_eval)()) # the log of negative terms come out with log(Real) + pi*1i
    } else {
      # label this term as positive
      term_sgn[i] <- 1
      
      # evaluate this term
      evaluation_list <- c(list(zInf = -Inf), evaluation_list)
      expr_eval <- subs(current_term, evaluation_list)
      cat("expr_eval: ", as_character(expr_eval), "\n")
      term_value[i] <- as_func(expr_eval)()
    }
  }
  
  #----------------------
  # sum-minus combination
  #----------------------
  if (sum(term_sgn == 1) == 0) {
    stop("The expression to be logged has no positive term.")
  } else if (sum(term_sgn == -1) == 0) {
    return(logplusvec(term_value))
  } else {
    negative_logterms <- logplusvec(term_value[term_sgn == -1])
    positive_logterms <- logplusvec(term_value[term_sgn == 1])
    if (negative_logterms > positive_logterms) {
      warning("the negative part of the expression is larger in magnitude than the positive part")
    }
    return(logminus(positive_logterms, negative_logterms))
  }
}
