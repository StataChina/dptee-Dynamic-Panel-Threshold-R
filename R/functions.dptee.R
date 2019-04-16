#' dptee internal functions
#' @description {Functions are: (1) lagnmatrix, (2) trimrmatrix, (3) dimen, (4) gmm, (5) estimate, (6) makeMF, and (7) makeW.}
#' @rdname functions.dptee
#' @aliases functions
#' @keywords internal
#' @export
functions.dptee = (dptee)

#' @rdname functions.dptee
#' @description Function (1)
#' @aliases functions
#' @keywords internal
#' @export
lagnmatrix = function(n, t, x, p) {
  nt = n * t
  k = dimen(x)$nc
  x = matrix(x, nt, k)
  k = ncol(x)
  lx = matrix(NA, nt, k)
  if (t == p) {
    i1 = i2 = c(outer(c(t:t), seq(0, nt - t, t), `+`))
    lx[i2,] = x[i1, ]
    return(lx)
  } else {
    i1 = c(outer(c(1:(t - p)), seq(0, nt - t, t), `+`))
    i2 = c(outer(c((p + 1):t), seq(0, nt - t, t), `+`))
    lx[i2,] = x[i1, ]
    return(lx)
  }
}

#' @rdname functions.dptee
#' @description Function (2)
#' @aliases functions
#' @keywords internal
#' @export
trimrmatrix  = function(n, t, x, p) {
  as.matrix(x)
  idx = -c(outer(c(1:p), seq(0, (n * t), t), `+`))
  y = x[idx, ]
  return(y)
}

#' @rdname functions.dptee
#' @description Function (3)
#' @aliases functions
#' @keywords internal
#' @export
dimen = function(x) {
  if (is.null(ncol(x))) {
    dimc = 1
    dimr = length(x)
  } else {
    dimc = ncol(x)
    dimr = nrow(x)
  }
  dim = list("nr" = dimr, "nc" = dimc)
  return(dim)
}

#' @rdname functions.dptee
#' @description Function (4)
#' @aliases functions
#' @keywords internal
#' @export
gmm = function(n, y, x, iv, w) {
  ivy = crossprod(iv, y) / n
  ivx = crossprod(iv, x) / n
  ivxw = crossprod(ivx, w)
  if (min(eigen(ivxw %*% ivx)$values) < 1e-8) {
    b = ginv(ivxw %*% ivx) %*% ivxw %*% ivy
  } else {
    b = solve(ivxw %*% ivx) %*% ivxw %*% ivy
  }
  u = ivy - ivx %*% b
  s = crossprod(u, w) %*% u
  resid = y - x %*% b
  gmmest = list("b" = c(b),
                "resid" = c(resid),
                "s" = c(s))
  return(gmmest)
}

#' @rdname functions.dptee
#' @description Function (5)
#' @aliases functions
#' @keywords internal
#' @export
estimate =
  function(coef,
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
           display.dummies) {
    k1 = dimen(z1)$nc
    k2 = dimen(z2)$nc
    k = k1 + k2
    kt = length(coef)
    id = diag(k2)
    bTh = coef[1]
    bU = coef[2:(k2 + 1)]
    bD = coef[(k2 + 2):(k + 1)]
    prop = mean(q > bTh)
    if (con == 0) {
      p1 = cbind(id, id)
      bL = bU + bD
      nP = c(
        paste0("Threshold(", names(bTh), ")"),
        paste0(names(bU), "(-)"),
        paste0(names(bD), "(d)"),
        paste0(names(bU), "(+)")
      )
    } else if (con == 1) {
      p1 = cbind(id, 0, id)
      bL = bU + bD[-1]
      nP = c(
        paste0("Threshold(", names(bTh), ")"),
        paste0(names(bU), "(-)"),
        paste0(names(bD), "(d)"),
        paste0(names(bU), "(+)")
      )
    }
    if (effects == "individual" & ki == 0) {
      bK = bT = NULL
      zU = cbind(z2, z1 * (q > bTh))
      nPk = NULL
    } else if (effects == "individual" & ki == 1) {
      bK = coef[(k + 2):kt]
      nPk = names(coef[(k + 2):kt])
      bT = NULL
      zU = cbind(z2, z1 * (q - bTh) * (q > bTh), xk)
    } else if (effects == "twoways" & ki == 0) {
      bL = bU + bD
      bK = NULL
      nPk = NULL
      bT = coef[((kt + 1) - ktt):kt]
      names(bT) = colnames(tt)
      zU = cbind(z2, z1 * (q - bTh) * (q > bTh), tt)
    } else if (effects == "twoways" & ki == 1) {
      bL = bU + bD
      bK = coef[(k + 2):(k + kxk + 1)]
      nPk = names(coef[(k + 2):(k + kxk + 1)])
      bT = coef[((kt + 1) - ktt):kt]
      names(bT) = colnames(tt)
      zU = cbind(z2, z1 * (q - bTh) * (q > bTh), xk, tt)
    }
    lzU = lagnmatrix(n, t, zU, 1)
    dzU = trimrmatrix(n, t, cbind(zU - lzU), ml)
    u = dy - dzU %*% coef[-1]
    w = makeW(n, iv, u = u)
    h = 1.06 * sd(q, na.rm = T) * n ^ (-0.2)
    kn = dnorm((bTh - q) / h)
    s = z1 * as.vector(kn)
    ls = lagnmatrix(n, t, s, 1)
    ds = trimrmatrix(n, t, (s - ls), ml)
    GbU = (crossprod(-iv, dzU)) / n
    Gr = (crossprod(-iv, ds) / (n * h)) %*% bD
    GU = cbind(Gr, GbU)
    vU = solve(crossprod(GU, w) %*% GU) / n
    sU = sqrt(diag(vU))
    sdTh = sU[1]
    sdU = sU[2:(k2 + 1)]
    sdD = sU[(k2 + 2):(k + 1)]
    sdL = sqrt(diag(p1 %*% vU[c(2:(k + 1)), c(2:(k + 1))] %*% t(p1)))
    if (effects == "individual" & ki == 0) {
      sdK = sdT = NULL
      wtc = as.numeric(crossprod(coef[-1], crossprod(solve(vU[-1, -1]), coef[-1])))
      dfwtc = length(coef[-1])
      pwtc = pchisq(wtc, df = dfwtc, lower.tail = FALSE)
      wtd = dfwtd = pwtd = NULL
    } else if (effects == "individual" & ki == 1) {
      sdK = sU[(k + 2):kt]
      sdT = NULL
      wtc = as.numeric(crossprod(coef[-1], crossprod(solve(vU[-1, -1]), coef[-1])))
      dfwtc = length(coef[-1])
      pwtc = pchisq(wtc, df = dfwtc, lower.tail = FALSE)
      wtd = dfwtd = pwtd = NULL
    } else if (effects == "twoways" & ki == 0) {
      sdK = NULL
      sdT = sU[(k + 2):kt]
      Rc = c(2:(k + 1))
      wtc = as.numeric(crossprod(coef[Rc], crossprod(solve(vU[Rc, Rc]), coef[Rc])))
      dfwtc = length(Rc)
      pwtc = pchisq(wtc, df = dfwtc, lower.tail = FALSE)
      Rd = c((kt - ktt + 1):kt)
      wtd = as.numeric(crossprod(coef[Rd], crossprod(solve(vU[Rd, Rd]), coef[Rd])))
      dfwtd = length(Rd)
      pwtd = pchisq(wtd, df = dfwtd, lower.tail = FALSE)
    } else if (effects == "twoways" & ki == 1) {
      sdK = sU[(k + 2):(kt - ktt)]
      sdT = sU[((kt - ktt) + 1):kt]
      Rc = c(2:((k + 1) + kxk))
      wtc = as.numeric(crossprod(coef[Rc], crossprod(solve(vU[Rc, Rc]), coef[Rc])))
      dfwtc = length(Rc)
      pwtc = pchisq(wtc, df = dfwtc, lower.tail = FALSE)
      Rd = c((kt - ktt + 1):kt)
      wtd = as.numeric(crossprod(coef[Rd], crossprod(solve(vU[Rd, Rd]), coef[Rd])))
      dfwtd = length(Rd)
      pwtd = pchisq(wtd, df = dfwtd, lower.tail = FALSE)
    }
    ivu = crossprod(iv, u)
    sa = as.numeric(crossprod(ivu, w) %*% ivu) / n
    if (display.dummies == T){
      b = c(bTh, bU, bD, bL, bK, bT)
      names(b) = (c(nP, nPk, colnames(tt)))
      sd = c(sdTh, sdU, sdD, sdL, sdK, sdT)
    } else {
      b = c(bTh, bU, bD, bL, bK)
      names(b) = (c(nP, nPk))
      sd = c(sdTh, sdU, sdD, sdL, sdK)
    }
    names(sd) = (c(nP, nPk))
    estim =
      list(
        "b" = b,
        "sd" = sd,
        "resid" = u,
        "prop" = prop,
        "sargan" = sa,
        "npar" = k,
        "wtc" = wtc,
        "dfwtc" = dfwtc,
        "pwtc" = pwtc,
        "wtd" = wtd,
        "dfwtc" = dfwtd,
        "pwtd" = pwtd
      )
    return(estim)
  }

#' @rdname functions.dptee
#' @description Function (6)
#' @aliases functions
#' @keywords internal
#' @export
makeMF = function(eq, data, effects, normal.inst) {
  if (!inherits(eq, "Formula"))
    stop('Need a Formula object. See Formula package.')
  if (!inherits(data, "pdata.frame"))
    stop('Need a pdata.frame objetc. See plm package.')
  n = pdim(data)$nT[[1]]
  t = pdim(data)$nT[[2]]
  nt = n * t
  mf = list()
  # Variables ---------------------------------------------------------------
  for (k in 1:3) {
    P = attr(terms(eq, rhs = k), "term.labels")
    if (!identical(P, character(0))) {
      nP = gsub("\\s", "\\1", P)
      Ps = list()
      for (i in 1:length(P)) {
        tvr = gsub("\\lag|\\(|", "", unlist(strsplit(nP[i], "\\,|\\)| ")))
        if (length(tvr) == 1 & substr(tvr[1], 1, 3) != "log") {
          vr = tvr[1]
          ln = NA
          il = 0
        } else if (length(tvr) == 1 &
                   substr(tvr[1], 1, 3) == "log") {
          vr = substr(tvr[1], 4, nchar(tvr[1]))
          ln = NA
          il = 1
        } else if (length(tvr) > 1 &
                   substr(tvr[1], 1, 3) != "log") {
          vr = tvr[1]
          ln = tvr[length(tvr)]
          il = 0
        } else if (length(tvr) > 1 &
                   substr(tvr[1], 1, 3) == "log") {
          vr = substr(tvr[1], 4, nchar(tvr[1]))
          ln = tvr[length(tvr)]
          il = 1
        } else {
          stop("Could not extract model frame.")
        }
        if (is.na(ln)) {
          if (il == 0) {
            Ps[[i]] = matrix(data[, vr])
            colnames(Ps[[i]]) = vr
          } else {
            Ps[[i]] = matrix(log(data[, vr]))
            colnames(Ps[[i]]) = paste0("log(", vr, ")")
          }
        }
        else if (!is.na(ln) & nchar(ln) == 1) {
          ll = as.numeric(unlist(strsplit(ln, ":")))
          if (il == 0) {
            Ps[[i]] = lagnmatrix(n, t, data[, vr], ll)
            colnames(Ps[[i]]) = paste0("lag(", vr, "," , ll, ")")
          } else{
            Ps[[i]] = lagnmatrix(n, t, log(data[, vr]), ll)
            colnames(Ps[[i]]) = paste0("lag(log(", vr, ")," , ll, ")")
          }
        }
        else if (!is.na(ln) & nchar(ln) > 1) {
          ll = as.numeric(unlist(strsplit(ln, ":")))
          ml1 = ll[1]
          ml2 = ll[2]
          Psi = matrix(0, nt, (ml2 - ml1 + 1))
          ii = 1
          if (il == 0) {
            for (j in ml1:ml2) {
              Psi[, ii] = lagnmatrix(n, t, data[, vr], j)
              ii = ii + 1
            }
            nPsi = c(ml1:ml2)
            for (jj in 1:(ml2 - ml1 + 1)) {
              nPsi[jj] = paste0("lag(", vr, "," , (jj + ml1 - 1), ")")
            }
            Ps[[i]] = Psi
            colnames(Ps[[i]]) = nPsi
          } else {
            for (j in ml1:ml2) {
              Psi[, ii] = lagnmatrix(n, t, log(data[, vr]), j)
              ii = ii + 1
            }
            nPsi = c(ml1:ml2)
            for (jj in 1:(ml2 - ml1 + 1)) {
              nPsi[jj] = paste0("lag(log(", vr, ")," , (jj + ml1 - 1), ")")
            }
            Ps[[i]] = Psi
            colnames(Ps[[i]]) = nPsi
          }
        }
      }
      tmp = do.call(cbind, Ps)
    } else {
      tmp = NULL
    }
    mf[[k]] = tmp
  }
  #--------------------------------------------------------------------------
  P1 = attr(terms(eq, rhs = 1), "term.labels")
  P4 = attr(terms(eq, rhs = 4), "term.labels")
  lP1 = length(P1)
  lP4 = length(P4)
  nP1 = gsub("\\s", "\\1", P1)
  nP4 = gsub("\\s", "\\1", P4)
  ln1 = ln4 = list()
  vr4 = rep(NA, lP4)
  for (i in 1:lP1) {
    ls1 = unlist(strsplit(nP1[i], "\\,|\\)| "))
    ln1[[i]] = suppressWarnings(as.numeric(unlist(strsplit(ls1[length(ls1)], ":"))))
  }
  l0 = max(unlist(ln1), na.rm = T) + 2   #Initial lag (l1)
  xx = matrix(0, n * (t - l0), lP4)
  x2 = list()
  dx2 = 0
  # Matrix of instruments ---------------------------------------------------
  for (i in 1:lP4) {
    l4 = unlist(strsplit(nP4[i], "\\,|\\)| "))
    ls4 = l4[length(l4)]
    ln4[[i]] = as.numeric(unlist(strsplit(ls4, ":")))
    vr4 = gsub("\\lag|\\(|", "", unlist(strsplit(nP4[i], "\\,|\\)| "))[1])
    if (substr(vr4, 1, 3) == "log")
      xl = lagnmatrix(n, t, log(data[, substr(vr4, 4, 100L)]), l0)
    else if (substr(vr4, 1, 3) != "log")
      xl = lagnmatrix(n, t, data[, substr(vr4, 1, 100L)], l0)
    xt = trimrmatrix(n, t, xl, l0)
    xi = lapply(seq((t - l0), nt, by = (t - l0)), function(j)
      xt[(j - (t - l0) + 1):j])
    if (is.na(ln4[[i]][1]))
      ln4[[i]][1] = 0
    if (ln4[[i]][1] <= l0) {
      minl = 1
    } else {
      minl = ln4[[i]][1] - l0
    }
    if (is.na(ln4[[i]][2]))
      ln4[[i]][2] = minl
    if (ln4[[i]][2] <= l0) {
      maxl = 1
    } else {
      maxl = ln4[[i]][2] - l0
    }
    x1 = list()
    for (ii in 1:n) {
      x0 = list()
      for (it in 1:(t - l0)) {
        if (minl != maxl & it < minl) {
          x0[[it]] = xi[[ii]][1]
        } else if (minl != maxl & it >= minl & it <= maxl) {
          x0[[it]] = xi[[ii]][1:(it - minl + 1)]
        } else if (minl != maxl & it > maxl) {
          x0[[it]] = xi[[ii]][(it - minl):(it - minl + 1)]
        } else if (minl == maxl) {
          x0[[it]] = xi[[ii]][it]
        }
      }
      x1[[ii]] = do.call(rbind, lapply(x0, `[`, seq_len(max(
        sapply(x0, length)
      ))))
    }
    x2[[i]] = do.call(rbind, x1)
    if (ncol(x2[[i]]) > dx2)
      dx2 = ncol(x2[[i]])
  }
  if (attr(terms(eq, rhs = 4), "intercept") & effects == "individual")
    x = 1
  else
    x = NULL
  for (i in 1:dx2) {
    for (j in 1:length(x2)) {
      exc = tryCatch(
        x2[[j]][, i],
        error = function(e)
          NULL
      )
      if (!is.null(exc))
        x = cbind(x, x2[[j]][, i])
    }
  }
  mm = length(which(!is.na(x[1:(t - l0), ])))
  sx = x[1:(t - l0), ]
  si = matrix(1, t - l0, 2)
  # Only way? ---------------------------------------------------------------
  if (dimen(sx)$nc == 1){
    si[1, 2] = length(na.omit(sx[1]))
    for (i1 in 2:(t - l0)) {
      si[i1, ]  = c(si[i1 - 1, 2] + 1, length(na.omit(sx[i1])) + si[i1 - 1, 2])
    }
    ixi = lapply(seq((t - l0), n * (t - l0), by = (t - l0)), function(i0) x[(i0 - (t - l0) + 1):i0,])
    iv = list()
    for (i2 in 1:n) {
      ivi = matrix(0, (t - l0), mm)
      for (t2 in 1:(t - l0)) {
        ivi[t2, c(si[t2, 1]:si[t2, 2])] = c(na.omit(ixi[[i2]][t2]))
        iv[[i2]] = ivi
      }
    }
  } else {
    si[1, 2] = length(na.omit(sx[1, ]))
    for (i1 in 2:(t - l0)) {
      si[i1, ]  = c(si[i1 - 1, 2] + 1, length(na.omit(sx[i1, ])) + si[i1 - 1, 2])
    }
    ixi = lapply(seq((t - l0), n * (t - l0), by = (t - l0)), function(i0) x[(i0 - (t - l0) + 1):i0,])
    iv = list()
    for (i2 in 1:n) {
      ivi = matrix(0, (t - l0), mm)
      for (t2 in 1:(t - l0)) {
        ivi[t2, c(si[t2, 1]:si[t2, 2])] = c(na.omit(ixi[[i2]][t2, ]))
        iv[[i2]] = ivi
      }
    }
  }
  # -------------------------------------------------------------------------
  iv = do.call(rbind, iv)
  niv = NULL
  for (i in 1:length(nP1)) {
    if (!any(gsub("\\lag\\(|\\,.*", "", nP4) == gsub("\\lag\\(|\\,.*", "", nP1[i]))) {
      nPI = gsub("\\lag\\(|\\,.*", "", nP1[i])
      if (substr(nPI, 1, 3) == "log")
        niv[[i]] = lagnmatrix(n, t, log(data[, gsub("\\log\\(|\\)|", "", nPI)]), l0)
      else
        niv[[i]] = lagnmatrix(n, t, data[, gsub("\\log\\(|\\)|", "", nPI)], l0)
    }
  }
  niv = do.call(cbind, niv)
  niv = trimrmatrix(n, t, niv, l0)
  if (normal.inst == T)
    mf[[4]] = cbind(iv, niv)
  else
    mf[[4]] = iv
  # Response ----------------------------------------------------------------
  yn = unlist(attr(eq, "lhs")[[1]])
  vr5 = gsub("\\lag|\\(|", "", yn)
  if (vr5[1] == "log") {
    y = log(data[, vr5[2]])
  } else if (substr(yn, 1, 3) != "log") {
    y = data[, vr5]
  }
  #--------------------------------------------------------------------------
  mf[[5]] = y
  return(list('mf' = mf))
}

#' @rdname functions.dptee
#' @description Function (7)
#' @aliases functions
#' @keywords internal
#' @export
makeW = function(n, iv, u = NULL) {
  nt = nrow(iv)
  t = nt / n
  v0 = v1 = v2 = 0
  liv = lagnmatrix(n, t, iv, 1)
  div = iv - liv
  snt = seq(1, nt, by = t)
  w1 = lapply(snt, function(i)
    iv[i:(i + t - 1),])
  u1 = lapply(snt, function(i)
    u[i:(i + t - 1)])
  if (is.null(u)) {
    for (i in 1:n) {
      for (j in 1:(t + 1)) {
        w0 = cbind(0, rbind(0, w1[[i]], 0), 0)
        wi = w0[j, ] - w0[j + 1, ]
        v0 = v0 + tcrossprod(wi, wi)
      }
    }
    v = (v0[-c(1, nrow(v0)),-c(1, ncol(v0))]) / n
    if (min(eigen(v)$values) < 1e-8) {
      w = ginv(v)
    } else {
      w = solve(v)
      return(w)
    }
  } else if (!is.null(u)) {
    for (i in 1:n) {
      wi = crossprod(w1[[i]], u1[[i]])
      v1 = v1 + tcrossprod(wi, wi)
      v2 = v2 + wi
    }
    v0 = (v1 / n) - tcrossprod((v2 / n), (v2 / n))
    if (min(eigen(v0)$values) < 1e-8) {
      w = ginv(v0)
    } else {
      w = solve(v0)
    }
  } else
    stop("Incorrect usage of the function.")
}
