summary.dptee =
  function (x,
            method = "analytic",
            digits = max(3L, getOption("digits") - 3L),
            signif.stars = getOption("show.signif.stars"),
            ...)
  {
    obj = x
    if (!inherits(obj, "dptee"))
      stop("argument 'object' needs to be class 'dptee'")
    met = method
    if (method == "averaging") {
      mod = 2
    } else if (method == "analytic") {
      mod = 1
    } else{
      mod = 0
      stop("Choose the method.")
    }
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
    prop = obj$result[[mod]][[4]]
    n = obj$n
    t = obj$t
    nt = n * t
    et = obj$et
    k = obj$result[[mod]][[4]]
    res = c(obj$result[[mod]][[3]])
    ntu = length(res)
    tu = ntu / n
    moments = obj$mom
    fresid = summary(res, digits = digits)
    if (ki == 0 & eff == "individual") {
      statistics = obj$result[[mod]][[1]]
      standev = obj$result[[mod]][[2]]
    } else if (ki == 1 & eff == "individual") {
      statistics = obj$result[[mod]][[1]]
      standev = obj$result[[mod]][[2]]
    } else if (ki == 0 & eff == "twoways") {
      statistics = obj$result[[mod]][[1]]
      statistics = statistics[-c((length(statistics) - (tu - 2)):length(statistics))]
      standev = obj$result[[mod]][[2]]
      standev = standev[-c((length(standev) - (tu - 2)):length(standev))]
    } else if (ki == 1 & eff == "twoways") {
      statistics = obj$result[[mod]][[1]]
      statistics = statistics[-c((length(statistics) - (tu - 2)):length(statistics))]
      standev = obj$result[[mod]][[2]]
      standev = standev[-c((length(standev) - (tu - 2)):length(standev))]
    }
    zval = statistics / standev
    pval = 2 * pnorm(-abs(zval))
    coefs <- matrix(NA, 0L, 4L)
    colnames(coefs) <-
      c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    coefs = cbind(
      `Estimate` = statistics, 3,
      `Std. Error` = standev,
      `z value` = zval,
      `Pr(>|z|)` = pval
    )
    P1 = gsub("\\s", "\\1", attr(terms(eq, rhs = 1), "term.labels"))
    P2 = gsub("\\s", "\\1", attr(terms(eq, rhs = 2), "term.labels"))
    P3 = gsub("\\s", "\\1", attr(terms(eq, rhs = 3), "term.labels"))
    T1 = paste0("Threshold (", P2, ")")
    T2 = paste0(P1, " (-)")
    if (con == 0)
      T3 = paste0(P1, " (d)")
    else
      T3 = c("constant", paste0(P1, " (d)"))
    if (ki == 0 & eff == "individual") {
      T4 = paste0(P1, " (+)")
    } else if (ki == 1 & eff == "individual") {
      T4 = c(paste0(P1, " (+)"), P3)
    } else if (ki == 0 & eff == "twoways") {
      T4 = paste0(P1, " (+)")
    } else if (ki == 1 & eff == "twoways") {
      T4 = c(paste0(P1, " (+)"), P3)
    }
    tabnames = c(T1, T2, T3, T4)
    rownames(coefs) = tabnames
    Jhat = obj$Jhat
    Jhat1 = obj$Jhat1
    pJhat = pchisq(Jhat, (moments - length(standev)), lower.tail = F)
    pJhat1 = pchisq(Jhat1, (moments - length(standev)), lower.tail = F)
    LinT = obj$LinT
    LinT1 = obj$LinT1
    wtc = obj$result[[mod]][7]
    dfwtc = obj$result[[mod]][8]
    pwtc = obj$result[[mod]][9]
    wtd = obj$result[[mod]][10]
    dfwtd = obj$result[[mod]][11]
    pwtd = obj$result[[mod]][12]
    #cat("Dynamic Panels with Threshold Effect and Endogeneity (Seo & Shin, 2016) \n\n")
    textcon =
      ifelse(con == 1,
             "(constant in the upper regime):",
             "(no constant in the upper regime):")
    textmet =
      ifelse(
        method == "analytic",
        "Output based on analytic method ",
        "Output based on averaging method "
      )
    cat("\nCall:\n",
        paste(deparse(obj$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = "")
    jstat = ifelse(method == "analytic", Jhat1, Jhat)
    pjstat = ifelse(method == "analytic", pJhat1, pJhat)
    lintest = ifelse(method == "analytic", LinT, LinT1)
    sa = obj$result[[mod]][[5]]
    pvalsa = pchisq(sa, df = (moments - k), lower.tail = FALSE)
    cat(textmet, textcon, "\n\n")
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
        formatC(jstat, format = 'f', digits = digits),
        " (p.value = ",
        formatC(pjstat, format = 'f', digits = digits),
        ")"
      ),
      collapse = ""
    ), "\n")
    cat(paste(c(
      "Linearity test (p.value): ",
      formatC(lintest, format = 'f', digits = digits)
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
