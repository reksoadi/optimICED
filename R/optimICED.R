gss <- function(f, lower, upper, tol = 0.0001, maxiter = 1000,
                showiter = FALSE, makinfo = TRUE) {
  # A bracketing function to peform golden-section search.
  #
  # Args:
  #   f: the function to be optimized.
  #   lower: the lower end point of the interval to be searched.
  #   upper: the upper end point of the interval to be searched.
  #   tol: the desired accuracy.
  #   maxiter: maximum iteration taken to run the bracketing function. Should
  #            be defined to prevent the case of infinite iteration.
  #   showiter: give description or information for each iteration.
  #   makinfo: give general information on function's results.
  #
  # Returns:
  #   A list of minimum (or maximum), total number of iteration, and objective
  #   which give the location of the minimum (or maximum). The description of
  #   result is also given.

  # Set the variable needed for the analysis.
  golden.ratio <- (sqrt(5) - 1)/2             # Define golden ratio.
  x1 <- upper - golden.ratio*(upper - lower)  # Set x1 as upper test point.
  x2 <- lower + golden.ratio*(upper - lower)  # Set x2 as lower test point.
  f1 <- f(x = x1, d = 0)                    # Evaluate x1 position in f(x).
  f2 <- f(x = x2, d = 0)                    # Evaluate x2 position in f(x).
  old.value <- (upper + lower)/2     # Give value of current (then old) minima.
  iteration <- 0                              # Set iteration counter.
  if (makinfo == TRUE)
    cat( "Golden Section Search\n")           # Name the optimization technique.

  # Put golden section search technique into loop. There are two conditions
  # that needs to be evaluate. First, if (f2 > f1), x2 becomes the new upper
  # bound, x1 becomes new x2, and the new upper test point becomes new x1.
  # Second, reversely, for (f2 < f1), x1 becomes the new upper bound,
  # x2 becomes new x1, and the new lower test point becomes new x2.
  repeat {
    iteration <- iteration + 1                    # Add iteration count.

    if (f2 > f1) {                                # First condition (f2 > f1).
      upper <- x2
      x2 <- x1
      f2 <- f1
      x1 <- upper - golden.ratio*(upper - lower)
      f1 <- f(x = x1)
    } else {                                      # Second condition (f2 < f1).
      lower <- x1
      x1 <- x2
      f1 <- f2
      x2 <- lower + golden.ratio*(upper - lower)
      f2 <- f(x = x2)
    }

    # Set initial optimization value (the midpoint between upper and
    # lower test point).
    initial.value <- (upper + lower)/2

    # We might want to keep the values in each iteration for function check.
    # This code will only run when 'showiter' is specified as TRUE.
    if ( showiter == TRUE ) {
      cat("Iteration ", iteration, ", lower = ", lower, ", upper = ",
          upper, ", value = ", initial.value, "\n")
    }

    # Put break to 'repeat' when optimization value (the difference between
    # upper and lower test point) is lower than predetermined tolerance 0.0001
    # (or any other value defined in 'tol').
    if (abs(upper - lower) < tol ) break

    # Also, put break when the number of iterations has achieved certain number
    # (in default case, 1000 iterations).
    if (iteration > maxiter ) break

    # Put the initial value as old value before the function repeats.
    old.value <- initial.value
  }
  # Give information on results.
  if (makinfo == TRUE) {
    cat("total number of iteration	= ", iteration, "\n")
    cat("last difference	= ", abs(initial.value - old.value), "\n")
    cat("final value		= ", initial.value, "\n")
    cat("final function value 	= ", f(x = initial.value), "\n")
  }

  # Return a list containing number of iteration and minima and/or maxima value
  # in the list.
  list("iter" = iteration, "value" = initial.value, "obj" = f(initial.value))
}


