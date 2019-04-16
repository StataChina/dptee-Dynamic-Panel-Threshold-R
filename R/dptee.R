dptee =
  function (formula,
            data,
            con = 0,
            trim = 0.15,
            qn = 200,
            ns = 300,
            jj = 1000,
            model = "analytical",
            effects = "individual",
            normal.inst = T,
            display.dummies = F,
            ...)
  {
    cl <- match.call()
    # Comp. time --------------------------------------------------------------
    start_time = Sys.time()
    # 1. Check packages needed ------------------------------------------------
    packlist = list("MASS", "plm", "Formula")
    for (i in 1:3) {
      if (is.element(packlist[[i]], installed.packages()[, 1]) == 0)
        stop("Please install ", packlist[[i]], " package.")
    }
    # 2. Check model formula and data -----------------------------------------
    eq = Formula(formula)
    data = data
    if (!inherits(eq, "Formula"))
      stop("Need a Formula equation. See Formula package.")
    if (!inherits(data, "pdata.frame"))
      data = pdata.frame(data)
    if (is.pbalanced(data) == 0)
      stop("Only balanced panels can be used!")
    if (length(eq)[2] < 2)
      stop("Incomplete model formulae!")
    if (terms(eq, rhs = 3)[[3]] == 0) {
      ki = 0
    } else{
      ki = 1
    }
    P1 = gsub("\\s", "\\1", attr(terms(eq, rhs = 1), "term.labels"))
    nP1 = paste(unlist(strsplit(P1, "\\,|\\:|\\)| ")))
    ml = max(suppressWarnings(as.numeric(nP1)), na.rm = T) + 2
    # 3. Extract data, variables, and parameters ------------------------------
    n = pdim(data)$nT[[1]]
    t = pdim(data)$nT[[2]]
    nt = n * t
    if (effects == "twoways") {
      tt = table(1:nt, factor(rep(1:t, n)))
      tt = tt[,-c(1:(ml + 1))]
      colnames(tt) = paste0("D", c(1:ncol(tt)))
    } else {
      tt = NULL
    }
    mf = makeMF(eq, data, effects, normal.inst = normal.inst)
    z2 = mf$mf[[1]] #covariates
    q = as.vector(mf$mf[[2]]) #threshold
    iv = mf$mf[[4]] #add. instruments
    if (effects == "twoways")
      iv = cbind(iv, trimrmatrix(n, t, tt, ml))
    y = mf$mf[[5]]
    ly = lagnmatrix(n, t, y, 1)
    dy = trimrmatrix(n, t, (y - ly), ml)
    if (isTRUE(anyNA(iv)) | isTRUE(anyNA(dy)))
      stop("Data must not contain NA values.")
    if (con == 0) {
      z1 = z2
    } else if (con == 1) {
      z1 = cbind(1, z2)
      colnames(z1) = c("constant", colnames(z2))
    } else {
      stop("Choose a model with or without constant.")
    }
    if (effects == "twoways" & con == 1)
      stop("The model does not make sense.")
    k1 = dimen(z1)$nc
    k2 = dimen(z2)$nc
    k = k1 + k2
    qq1 =
      quantile(unique(q),  probs = seq(trim, (1 - trim), by = (1 / qn)), names =
                 F, na.rm = T)
    qn1 = length(qq1)
    dz = NULL
    if (effects == "individual" & ki == 0) {
      xk = 0
      kxk = 0
      ktt = 0
      for (i in 1:length(qq1)) {
        z = cbind(z2, z1 * (q > qq1[i]))
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      }
    } else if (effects == "individual" & ki == 1) {
      xk = mf$mf[[3]] # kink vars.
      kxk = dimen(xk)$nc
      ktt = 0
      for (i in 1:length(qq1)) {
        z =
          cbind(z2, z1 * (q - qq1[i]) * (q > qq1[i]), xk)
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      }
    } else if (effects == "twoways" & ki == 0) {
      xk = 0
      kxk = 0
      ktt = dimen(tt)$nc
      for (i in 1:length(qq1)) {
        z =
          cbind(z2, z1 * (q - qq1[i]) * (q > qq1[i]), tt)
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      }
    } else if (effects == "twoways" & ki == 1) {
      xk = mf$mf[[3]] # kink vars.
      kxk = dimen(xk)$nc
      ktt = dimen(tt)$nc
      for (i in 1:length(qq1)) {
        z =
          cbind(z2, z1 * (q - qq1[i]) * (q > qq1[i]), xk, tt)
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      }
    }
    npar = dimen(dz[[1]])$nc
    npar1 = 2 + npar + 1
    p1 = cbind(cbind(matrix(0, k1, k2)), diag(k1))
    ntu = dimen(iv)$nr
    mom = dimen(iv)$nc
    zst = replicate(mom, rnorm(jj))
    coefs = matrix(0, ns, (1 + npar))
    Js = colp = rep(0, ns)
    wn_A = wnst_A = 0
    # 4. Start loop for models ------------------------------------------------
    for (it in 1:ns) {
      if (it == 1) {
        w1 = makeW(n, iv)
      } else{
        w1 = makeW(n, iv, u = rnorm(ntu))
      }
      col1 = matrix(0, qn1, (2 + npar + ntu))
      col2 = col1
      for (i in 1:qn1) {
        gmm1 = gmm(n, dy, dz[[i]], iv, w1)
        col1[i, ] = c(gmm1$s, qq1[i], gmm1$b, gmm1$resid)
      }
      ccol1 = dimen(col1)$nc
      col1hat = col1[order(col1[, 1]), ]
      u1 = col1hat[1, npar1:ccol1]
      w2 = makeW(n, iv, u = u1)
      wnst = matrix(0, qn1, jj)
      wn = matrix(0, qn1, 1)
      for (i in 1:qn1) {
        gmm2 = gmm(n, dy, dz[[i]], iv, w2)
        col2[i, ] = c(gmm2$s, qq1[i], gmm2$b, gmm2$resid)
        #-- supW statistic
        gmm2b = gmm2$b
        if (effects == "individual" & ki == 0) {
          p1b = p1 %*% gmm2b
          ivdz = crossprod(iv, dz[[i]]) / n
        } else if (effects == "individual" & ki == 1) {
          p1b = p1 %*% gmm2b[-(k + 1):-(k + kxk)]
          ivdz = crossprod(iv, dz[[i]][, -(k + 1):-(k + kxk)]) / n
        } else if (effects == "twoways" & ki == 0) {
          p1b = p1 %*% gmm2b[-(k + 1):-(k + ktt)]
          ivdz = crossprod(iv, dz[[i]][, -(k + 1):-(k + ktt)]) / n
        } else if (effects == "twoways" & ki == 1) {
          p1b = p1 %*% gmm2b[-(k + 1):-(k + kxk + ktt)]
          ivdz = crossprod(iv, dz[[i]][, -(k + 1):-(k + kxk + ktt)]) / n
        }
        u2 = col1[i, npar1:ccol1]
        iwi = makeW(n, iv, u = u2)
        ciwi = chol(iwi)
        ivdzt = crossprod(ivdz, iwi) %*% ivdz
        if (min(eigen(ivdzt)$values) < 1e-8) {
          ivdzi = ginv(ivdzt)
        } else {
          ivdzi = solve(ivdzt)
        }
        idtt = p1 %*% tcrossprod(ivdzi, p1)
        if (min(eigen(idtt)$values) < 1e-8) {
          idti = ginv(idtt)
        } else {
          idti = solve(idtt)
        }
        p1idti = crossprod(p1, idti) %*% p1
        wn[i] = n * (crossprod(p1b, idti) %*% p1b)
        mt = ivdz %*% ivdzi %*% p1idti %*% tcrossprod(ivdzi, ivdz)
        for (j in 1:jj) {
          wnst[i, j] =
            crossprod(zst[j, ], ciwi) %*% mt %*% crossprod(ciwi, zst[j, ])
        }
      }
      col2hat = col2[order(col2[, 1]), ]
      dehat = col2hat[1, npar1:dimen(col2hat)$nc]
      coefs[it, ] = col2hat[1, 2:(2 + npar)]
      Js[it] = n * col2hat[1, 1]
      colp[it] =
        1 - length(which(apply(wnst, 2, function(x)
          max(x, na.rm = TRUE)) <= max(wn))) / jj
      wnst_A = wnst_A + wnst
      wn_A = wn_A + wn
    }
    # 5. Extract and format results -------------------------------------------
    if (model == "analytical") {
      coef = coefs[1, ]
      Jhat = Js[1]
      LinT = colp[1]
    } else if (model == "averaging") {
      coef = colMeans(coefs)
      Jhat = mean(Js)
      LinT =
        1 - length(which(apply(wnst_A, 2, function(x)
          max(x, na.rm = TRUE)) <= max(wn_A))) / jj
    } else {
      stop("Choose a model.")
    }
    names(coef) = c(colnames(mf$mf[[2]]),
                    colnames(z2),
                    colnames(z1),
                    colnames(xk),
                    colnames(tt))
    result = estimate(coef,
                      q,
                      eq,
                      z2,
                      z1,
                      n,
                      t,
                      ml,
                      dy,
                      iv,
                      con,
                      xk,
                      kxk,
                      effects,
                      tt,
                      ktt,
                      ki,
                      display.dummies)
    # Comp. time --------------------------------------------------------------
    end_time = Sys.time()
    et = end_time - start_time
    # 6. Export results -------------------------------------------------------
    dptee =
      list(
        "formula" = eq,
        "data" = data,
        "n" = n,
        "t" = t,
        "result" = result,
        "et" = et,
        "Jhat" = Jhat,
        "LinT" = LinT,
        "mom" = mom,
        "con" = con,
        "effects" = effects,
        "call" = cl
      )
    class(dptee) = "dptee"
    return(dptee)
  }
