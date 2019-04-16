summary.dptee =
  function (object,
            digits = max(3L, getOption("digits") - 3L),
            signif.stars = getOption("show.signif.stars"),
            ...)
  {
    obj = object
    if (!inherits(obj, "dptee"))
      stop("argument 'object' needs to be class 'dptee'")
    eq = obj$formula
    if (terms(Formula(eq), rhs = 3)[[3]] != 0) {
      ki = 1
      kl = length(attr(terms(eq, rhs = 3), "term.labels"))
    } else{
      ki = 0
      kl = 0
    }
    con = obj$con
    eff = obj$effects
    prop = obj$result[[4]]
    n = obj$n
    t = obj$t
    nt = n * t
    et = obj$et
    k = obj$result[[6]]
    res = c(obj$result[[3]])
    ntu = length(res)
    moments = obj$mom
    fresid = summary(res, digits = digits)
    statistics = obj$result[[1]]
    standev = obj$result[[2]]
    zval = statistics / standev
    pval = 2 * pnorm(-abs(zval))
    coefs <- matrix(NA, 0L, 4L)
    colnames(coefs) <-
      c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    coefs = cbind(
      `Estimate` = statistics,
      `Std. Error` = standev,
      `z value` = zval,
      `Pr(>|z|)` = pval
    )
    Jhat = obj$Jhat
    pJhat = pchisq(Jhat, (moments - length(standev)), lower.tail = F)
    LinT = obj$LinT
    wtc = obj$result[7]
    dfwtc = obj$result[8]
    pwtc = obj$result[9]
    wtd = obj$result[10]
    dfwtd = obj$result[11]
    pwtd = obj$result[12]
    #cat("Dynamic Panels with Threshold Effect and Endogeneity (Seo & Shin, 2016) \n\n")
    cat("\nCall:\n",
        paste(deparse(obj$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = "")
    sa = obj$result[[5]]
    pvalsa = pchisq(sa, df = (moments - k), lower.tail = FALSE)
    cat(paste(c("Panel:", " n = ", n, ", T = ", t, ", nT = ", nt), collapse =
                ""), "\n")
    cat(paste(c("# obs. Used: ", ntu), collapse = ""), "\n")
    cat(paste(c("# of Moments: ", moments), collapse = ""), "\n\n")
    cat("Residuals \n")
    print(fresid)
    cat("\nCoefficients:\n")
    printCoefmat(
      coefs,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
    cat("\n")
    cat(
      "Proportion of upper and lower regimes:",
      formatC(prop, format = 'f', digits = digits),
      "and" ,
      formatC(1 - prop, format = 'f', digits = digits),
      "\n"
    )
    cat(paste(
      c(
        "J-statistic: ",
        formatC(Jhat, format = 'f', digits = digits),
        " (p.value = ",
        formatC(pJhat, format = 'f', digits = digits),
        ")"
      ),
      collapse = ""
    ), "\n")
    cat(paste(c(
      "Linearity test (p.value): ",
      formatC(LinT, format = 'f', digits = digits)
    ), collapse = ""), "\n")
    cat(
      paste0(
        "Sargan test: chisq(",
        as.integer(moments - k),
        ") = ",
        formatC(sa, format = 'f', digits = digits),
        " (p.value=",
        formatC(pvalsa, format = 'f', digits = digits),
        ") \n"
      )
    )
    cat(
      paste0(
        "Wald test (coefficients): chisq(",
        unlist(dfwtc),
        ") = ",
        formatC(unlist(wtc), format = 'f', digits = digits),
        " (p.value=",
        formatC(unlist(pwtc), format = 'f', digits = digits),
        ") \n"
      )
    )
    if (eff == "twoways") {
      cat(
        paste0(
          "Wald test (time dummies): chisq(",
          unlist(dfwtd),
          ") = ",
          formatC(unlist(wtd), format = 'f', digits = digits),
          " (p.value=",
          formatC(unlist(pwtd), format = 'f', digits = digits),
          ") \n"
        )
      )
    }
    class(obj) = "summary.dptee"
  }
