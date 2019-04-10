dptee =
  function (formula,
            data,
            con = 0,
            trim = 0.1,
            qn = 20,
            ns = 20,
            jj = 20,
            effects = "individual",
            ...)
  {
    # Comp. time --------------------------------------------------------------
    start_time = Sys.time()
    # 1. Check packages needed ------------------------------------------------
    packlist = list("MASS", "plm", "Formula")
    for (i in 1:3) {
      if (is.element(packlist[[i]], installed.packages()[, 1]) == 0)
        stop("Please install ", packlist[[i]], " package.")
    }
    # 2. Check model formula and data -----------------------------------------
    data = pdata.frame(data)
    if (!inherits(eq, "Formula"))
      stop('Need a Formula equation. See Formula package.')
    if (!inherits(data, "pdata.frame"))
      stop('Need a pdata.frame objetc. See plm package.')
    if (is.pbalanced(data) == 0)
      stop("Only balanced panels can be used!")
    eq = Formula(eq)
    if (length(eq)[2] < 2)
      stop("Incomplete model formulae!")
    if (terms(Formula(eq), rhs = 3)[[3]] != 0) {
      ki = 1
    } else{
      ki = 0
    }
    P1 = gsub("\\s", "\\1", attr(terms(eq, rhs = 1), "term.labels"))
    nP1 = unlist(strsplit(P1, "\\,|\\)| "))[2]
    ml = max(as.numeric(unlist(strsplit(nP1, ":")))) + 2
    # 3. Extract data, variables, and parameters ------------------------------
    n = pdim(data)$nT[[1]]
    t = pdim(data)$nT[[2]]
    nt = n * t
    tt = table(1:nt, factor(rep(1:t, n)))
    model = makeMF(eq, data)
    z2 = model$mf[[1]] #covariates
    q = model$mf[[2]] #threshold
    q = as.vector(q)
    if (ki == 0) {
      xk = 0
      kxk = 0
    } else if (ki == 1) {
      xk = model$mf[[3]] # kink vars.
      kxk = dimen(xk)$nc
    }
    if (effects == "twoways") {
      ki = 1
      if (xk == 0) {
        xk = tt[, -c(1:(ml + 1))]
        kxk = dimen(tt)$nc
      } else {
        xk = cbind(xk, tt[, -c(1:(ml + 1))])
        kxk = dimen(tt)$nc
      }
    }
    xx = model$mf[[4]] #add. instruments
    y = attr(terms(Formula(eq)), "variables")[[2]]
    y = data$y
    ly = lagnmatrix(n, t, y, 1)
    dy = trimrmatrix(n, t, (y - ly), ml)
    iv = makeIV(eq, data)
    if (con == 0) {
      z1 = z2
    } else if (con == 1) {
      z1 = cbind(z2, 1)
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
                 F)
    qn1 = length(qq1)
    dz = NULL
    for (i in 1:length(qq1)) {
      if (ki == 0) {
        z = cbind(z2, z1 * (q > qq1[i]))
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      } else if (ki == 1) {
        z =
          cbind(z2, z1 * (q - qq1[i]) * q > qq1[i], xk)
        lz = lagnmatrix(n, t, z, 1)
        dz[[i]] = trimrmatrix(n, t, (z - lz), ml)
      }
    }
    npar = dimen(dz[[1]])$nc
    npar1 = 2 + npar + 1
    ntu = dimen(iv)$nr
    mom = dimen(iv)$nc
    zst = replicate(mom, rnorm(jj))
    coefs = matrix(0, ns, (1 + npar))
    Js = colp = rep(0, ns)
    wn_A = wnst_A = 0
    # 4. Start loop for averaging method --------------------------------------
    for (it in 1:ns) {
      if (it == 1) {
        w = makeW(n, iv)
      } else{
        w = makeW(n, iv, u = rnorm(ntu))
      }
      col1 = matrix(0, qn1, (2 + npar + ntu))
      col2 = col1
      for (i in 1:qn1) {
        gmm1 = gmm(n, dy, dz[[i]], iv, w)
        col1[i,] = c(gmm1$s, qq1[i], gmm1$b, gmm1$resid)
      }
      ccol1 = dimen(col1)$nc
      col1hat = col1[order(col1[, 1]),]
      u1 = col1hat[1, npar1:ccol1]
      w = makeW(n, iv, u = u1)
      wnst = matrix(0, qn1, jj)
      wn = matrix(0, qn1, 1)
      for (i in 1:qn1) {
        gmm2 = gmm(n, dy, dz[[i]], iv, w)
        col2[i,] = c(gmm2$s, qq1[i], gmm2$b, gmm2$resid)
        #-- supW statistic
        gmm2b = gmm2$b
        if (ki == 0) {
          p1 = cbind(cbind(matrix(0, k1, k2)), diag(k1))
          p1b = p1 %*% gmm2b
          ivdz = crossprod(iv, dz[[i]]) / n
        } else if (ki == 1) {
          p1 = cbind(cbind(matrix(0, k1, k2)), diag(k1))
          p1b = p1 %*% gmm2b[-(k + 1):-(k + kxk)]
          ivdz = crossprod(iv, dz[[i]][,-(k + 1):-(k + kxk)]) / n
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
            crossprod(zst[j,], ciwi) %*% mt %*% crossprod(ciwi, zst[j,])
        }
      }
      col2hat = col2[order(col2[, 1]),]
      dehat = col2hat[1, npar1:dimen(col2hat)$nc]
      coefs[it,] = col2hat[1, 2:(2 + npar)]
      Js[it] = n * col2hat[1, 1]
      # LinT1 = colp[1]
      colp[it] =
        1 - length(which(apply(wnst, 2, function(x)
          max(x, na.rm = TRUE)) <= max(wn))) / jj
      wnst_A = wnst_A + wnst
      wn_A = wn_A + wn
    }
    # 5. Extract and format results -------------------------------------------
    coef1 = coefs[1,]
    coef = colMeans(coefs)
    Jhat1 = Js[1]
    Jhat = mean(Js)
    # LinT
    ltp_A =
      1 - length(which(apply(wnst_A, 2, function(x)
        max(x, na.rm = TRUE)) <= max(wn_A))) / jj
    result = list()
    method =
      list(
        an = list(coef1, q, eq, z2, z1, n, t, ml, dy, iv, con, xk, kxk, effects),
        av = list(coef, q, eq, z2, z1, n, t, ml, dy, iv, con, xk, kxk, effects)
      )
    for (type in 1:2) {
      result[[type]] = estimate(
        method[[type]][[1]],
        method[[type]][[2]],
        method[[type]][[3]],
        method[[type]][[4]],
        method[[type]][[5]],
        method[[type]][[6]],
        method[[type]][[7]],
        method[[type]][[8]],
        method[[type]][[9]],
        method[[type]][[10]],
        method[[type]][[11]],
        method[[type]][[12]],
        method[[type]][[13]],
        method[[type]][[14]]
      )
    }
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
        "Jhat1" = Jhat1,
        "LinT" = ltp_A,
        "LinT1" = colp[1],
        "mom" = mom,
        "con" = con,
        "effects" = effects
      )
    class(dptee) = "dptee"
    return(dptee)
  }
