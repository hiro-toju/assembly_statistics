
# ssid-wise energy profile on de (factor) grid: NA where absent
.build_energy_by_ssid <- function(gela, ssids) {
  de <- gela[[2]]
  E <- lapply(ssids, function(j) {
    vapply(seq_along(de), function(i) {
      ss_i <- gela[[1]][[i]][[1]][[1]]
      idx  <- match(j, ss_i)
      if (is.na(idx)) return(NA_real_)
      as.numeric(gela[[1]][[i]][[1]][[2]][idx])
    }, numeric(1))
  })
  names(E) <- ssids
  E
}

# Longest consecutive run in a logical vector
.longest_run <- function(x) {
  r <- rle(x)
  if (!length(r$lengths)) return(0L)
  max(r$lengths[r$values], 0L)
}

# Helper: first/last index per ssid based on finite entries in E[[ssid]]
.ssid_ranges_from_E <- function(E) {
  P <- lapply(E, function(e) is.finite(e))
  rng <- lapply(P, function(p) {
    idx <- which(p)
    if (!length(idx)) c(Inf, -Inf) else c(min(idx), max(idx))
  })
  names(rng) <- names(E)
  rng
}

# Build directed "adjacent in y" edges: last(from)+1 == first(to)
.build_adj_edges <- function(ssids, rng) {
  n <- length(ssids)
  from <- character(0)
  to   <- character(0)

  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) next
    a <- ssids[i]; b <- ssids[j]
    ra <- rng[[a]]; rb <- rng[[b]]
    if (is.finite(ra[2]) && is.finite(rb[1]) && (ra[2] + 1L == rb[1])) {
      from <- c(from, a)
      to   <- c(to, b)
    }
  }
  data.frame(from = from, to = to, stringsAsFactors = FALSE)
}

# Topological order that respects adjacency direction (left -> right).
# Tie-break: smaller first-index goes more left.
.topo_order <- function(nodes, edges, rng) {
  nodes <- unique(nodes)

  # adjacency list and indegree
  adj <- setNames(vector("list", length(nodes)), nodes)
  indeg <- setNames(rep(0L, length(nodes)), nodes)

  if (nrow(edges)) {
    for (k in seq_len(nrow(edges))) {
      f <- edges$from[k]; t <- edges$to[k]
      if (!(f %in% nodes) || !(t %in% nodes)) next
      adj[[f]] <- c(adj[[f]], t)
      indeg[[t]] <- indeg[[t]] + 1L
    }
  }

  first_idx <- setNames(vapply(nodes, function(x) rng[[x]][1], integer(1)), nodes)

  out <- character(0)
  remaining <- nodes

  while (length(remaining)) {
    cand <- remaining[indeg[remaining] == 0L]
    if (!length(cand)) {
      # Cycle or inconsistent constraints: fallback to sorting by first index
      cand <- remaining
    }
    pick <- cand[order(first_idx[cand], cand)][1]
    out <- c(out, pick)

    remaining <- setdiff(remaining, pick)
    for (t in adj[[pick]]) {
      if (t %in% remaining) indeg[[t]] <- max(0L, indeg[[t]] - 1L)
    }
    indeg[[pick]] <- 0L
  }

  out
}

# Reorder ssids so that adjacency chains always go to the right:
# if last(a)+1==first(b), then a must appear left of b.
# Components are ordered by earliest first-index.
.order_ssid_rightward_by_y <- function(ssids, E) {
  rng <- .ssid_ranges_from_E(E)
  edges <- .build_adj_edges(ssids, rng)

  # If no adjacency edges, just sort by first index
  if (!nrow(edges)) {
    first_idx <- vapply(ssids, function(x) rng[[x]][1], integer(1))
    return(ssids[order(first_idx, ssids)])
  }

  # Undirected adjacency for weak components
  und <- setNames(vector("list", length(ssids)), ssids)
  for (k in seq_len(nrow(edges))) {
    a <- edges$from[k]; b <- edges$to[k]
    und[[a]] <- c(und[[a]], b)
    und[[b]] <- c(und[[b]], a)
  }

  # BFS components
  seen <- setNames(rep(FALSE, length(ssids)), ssids)
  comps <- list()
  for (v in ssids) {
    if (seen[[v]]) next
    q <- c(v); seen[[v]] <- TRUE
    comp <- character(0)
    while (length(q)) {
      cur <- q[1]; q <- q[-1]
      comp <- c(comp, cur)
      nei <- unique(und[[cur]])
      for (u in nei) {
        if (!seen[[u]]) { seen[[u]] <- TRUE; q <- c(q, u) }
      }
    }
    comps[[length(comps) + 1]] <- comp
  }

  # Order components by earliest first-index
  comp_first <- vapply(comps, function(cc) {
    min(vapply(cc, function(x) rng[[x]][1], integer(1)))
  }, integer(1))
  comps <- comps[order(comp_first)]

  # Topo-order within each component
  ordered <- unlist(lapply(comps, function(cc) {
    sub_edges <- edges[edges$from %in% cc & edges$to %in% cc, , drop = FALSE]
    .topo_order(cc, sub_edges, rng)
  }))

  # Append any missing nodes (safety)
  missing <- setdiff(ssids, ordered)
  if (length(missing)) {
    first_idx <- vapply(missing, function(x) rng[[x]][1], integer(1))
    ordered <- c(ordered, missing[order(first_idx, missing)])
  }

  ordered
}

# Distance matrix using adjacency (gap==0) reward + boundary energy + membership
.build_D_ssid <- function(ssids, E, nspecies,
                          w_adj = 10.0,     # weight for adjacency reward (bigger = more force)
                          w_energy = 2.0,
                          w_mem = 1.0,
                          w_gap = 2.0,      # penalty for larger gaps (non-adjacent)
                          adj_energy_scale = c("median","iqr","sd"),
                          adj_only = TRUE   # if TRUE: adjacency dominates; if FALSE: also uses overlap when exists
) {
  adj_energy_scale <- match.arg(adj_energy_scale)
  n <- length(ssids)
  D <- matrix(0, n, n, dimnames = list(ssids, ssids))

  # Presence and interval endpoints
  P <- lapply(E, function(e) is.finite(e))
  rng <- lapply(P, function(p) {
    idx <- which(p)
    if (!length(idx)) c(Inf, -Inf) else c(min(idx), max(idx))
  })

  # Membership vectors
  B <- lapply(ssids, function(j) id2bin(j, nspecies))
  names(B) <- ssids

  # Collect boundary energy diffs for scaling
  bd_diffs <- c()

  # Helper: boundary gap and boundary energy diff
  boundary_metrics <- function(sa, sb) {
    ra <- rng[[sa]]; rb <- rng[[sb]]
    ea <- E[[sa]];   eb <- E[[sb]]

    # Determine temporal order by interval position (on index axis)
    # If disjoint: earlier-last < later-first
    if (ra[2] < rb[1]) {
      gap <- rb[1] - ra[2] - 1
      dE  <- abs(ea[ra[2]] - eb[rb[1]])
      return(list(gap = gap, dE = dE))
    }
    if (rb[2] < ra[1]) {
      gap <- ra[1] - rb[2] - 1
      dE  <- abs(eb[rb[2]] - ea[ra[1]])
      return(list(gap = gap, dE = dE))
    }

    # Overlapping intervals (or touching in a complex way): use nearest boundary-like points
    ia <- which(is.finite(ea)); ib <- which(is.finite(eb))
    if (!length(ia) || !length(ib)) return(list(gap = Inf, dE = Inf))

    da <- outer(ia, ib, function(x, y) abs(x - y))
    k  <- which(da == min(da), arr.ind = TRUE)[1,]
    gap <- max(0, abs(ia[k[1]] - ib[k[2]]) - 1)
    dE  <- abs(ea[ia[k[1]]] - eb[ib[k[2]]])
    list(gap = gap, dE = dE)
  }

  # First pass: compute boundary metrics + membership diffs
  gap_mat <- matrix(Inf, n, n)
  dE_mat  <- matrix(Inf, n, n)
  dM_mat  <- matrix(0,   n, n)

  for (a in seq_len(n)) for (b in a:n) {
    if (a == b) next
    sa <- ssids[a]; sb <- ssids[b]

    bm <- boundary_metrics(sa, sb)
    gap_mat[a,b] <- gap_mat[b,a] <- bm$gap
    dE_mat[a,b]  <- dE_mat[b,a]  <- bm$dE
    bd_diffs <- c(bd_diffs, bm$dE)

    dM <- sum(B[[sa]] != B[[sb]]) / nspecies
    dM_mat[a,b] <- dM_mat[b,a] <- dM
  }

  # Scale boundary energy diffs (so energy doesn't dominate)
  sc <- if (!length(bd_diffs)) 1 else {
    bd <- bd_diffs[is.finite(bd_diffs)]
    if (!length(bd)) 1
    else if (adj_energy_scale == "median") stats::median(bd)
    else if (adj_energy_scale == "iqr") stats::IQR(bd)
    else stats::sd(bd)
  }
  if (!is.finite(sc) || sc <= 0) sc <- 1
  dE_scaled <- dE_mat / sc

  # Build final distance:
  for (a in seq_len(n)) for (b in a:n) {
    if (a == b) next
    gap <- gap_mat[a,b]
    dE  <- dE_scaled[a,b]
    dM  <- dM_mat[a,b]

    # adjacency reward: gap==0 should be strongly favored
    adj_reward <- if (is.finite(gap) && gap == 0) w_adj else 0

    # gap penalty (0 if adjacent)
    gap_pen <- if (is.finite(gap)) w_gap * (gap / 1) else w_gap * 100

    # base distance components
    d <- w_energy * dE + w_mem * dM + gap_pen

    # apply reward (make distance smaller)
    d <- d - adj_reward

    # if adjacency should dominate, clip adjacent pairs to be extremely small (tie-break by dE/dM)
    if (adj_only && is.finite(gap) && gap == 0) {
      d <- min(d, 1e-6 + w_energy * dE + w_mem * dM)
    }

    D[a,b] <- D[b,a] <- d
  }

  # Ensure finite
  finite_vals <- D[is.finite(D) & D > 0]
  big <- if (length(finite_vals)) max(finite_vals) * 10 else 1e6
  D[!is.finite(D)] <- big

  D
}

# 1D ordering: nearest-neighbor greedy + simple 2-opt refinement
.order_by_D <- function(D) {
  n <- nrow(D)

  # Start from the node with smallest mean distance
  start <- which.min(rowMeans(D))
  remaining <- setdiff(seq_len(n), start)
  ord <- start

  # Greedy chain construction
  while (length(remaining)) {
    last <- ord[length(ord)]
    nxt <- remaining[which.min(D[last, remaining])]
    ord <- c(ord, nxt)
    remaining <- setdiff(remaining, nxt)
  }

  # 2-opt refinement to reduce total adjacent cost
  cost <- function(o) sum(D[cbind(o[-length(o)], o[-1])])
  best <- ord
  best_cost <- cost(best)

  improved <- TRUE
  while (improved) {
    improved <- FALSE
    for (i in 2:(n-2)) for (k in (i+1):(n-1)) {
      cand <- best
      cand[i:k] <- rev(cand[i:k])
      cst <- cost(cand)
      if (cst + 1e-12 < best_cost) {
        best <- cand
        best_cost <- cst
        improved <- TRUE
      }
    }
  }

  best
}

# Build an ssid order that:
#  (1) Uses distance-based ordering for "who should be adjacent"
#  (2) Enforces that y-adjacent chains go to the right (larger y => larger x direction)
.order_ssid_with_rightward_adjacency <- function(gela, ssid, nspecies,
                                                w_adj = 10,
                                                w_energy = 2,
                                                w_mem = 1,
                                                w_gap = 2,
                                                adj_only = TRUE,
                                                adj_energy_scale = "median") {
  # energy profiles
  E <- .build_energy_by_ssid(gela, ssid)

  # distance-based adjacency preference
  D <- .build_D_ssid(ssid, E, nspecies,
                     w_adj = w_adj,
                     w_energy = w_energy,
                     w_mem = w_mem,
                     w_gap = w_gap,
                     adj_energy_scale = adj_energy_scale,
                     adj_only = adj_only)

  ord_idx <- .order_by_D(D)
  ssid_nn <- ssid[ord_idx]

  # enforce "rightward by y" within adjacency chains (gap==0 edges)
  ssid_final <- .order_ssid_rightward_by_y(ssid_nn, E)

  ssid_final
}

    
GELSObj2 <- function (gela, sa, threads = 1)
{
    itr <- 10000
    de <- gela[[2]]
    s <- length(sa[1][[1]][, 1])
    refenv <- gela[[3]]
    elsurf <- foreach(x = gela[[1]]) %do% {
        GraphObj(x[[1]])
    }
    ssp <- lapply(elsurf, function(x) {
        x[which(sapply(x$ccc_str, function(y) length(eval(parse(text = y))) ==
            1)), ]
    })
    sssp <- foreach(x = ssp, y = de) %do% {
        mat <- matrix(c(rep(y, nrow(x)), x[, 3], x[, 5]), nrow(x))
        split(mat, seq_len(nrow(mat)))
    }
    sssp <- unlist(sssp, recursive = FALSE)
    sslist <- unique(unlist(lapply(ssp, function(df) df[, 5])))
    if (length(sslist) > 1) {
        ss_combi <- combn(sslist, 2)
        mm <- apply(ss_combi, 2, function(pair) {
            dist <- hamming_distance(id2bin(pair[1], s), id2bin(pair[2],
                s))
            if (dist == 1)
                return(pair)
        })
        mm <- Filter(function(x) !is.null(x), mm)
    }
    else {
        mm <- NULL
    }
    concompo <- gen_concompo(mm)
    grp <- c(concompo, sapply(setdiff(sslist, unlist(concompo)),
        function(x) list(x)))
    ssspa <- lapply(grp, function(gg) {
        lapply(sssp, function(ss) {
            ss[grepl(paste(gg, collapse = "|"), ss[3])]
        })
    })
    ssspa <- lapply(ssspa, function(x) x[sapply(x, function(y) length(y) >
        0)])
    if (length(grp) == 1) {
        ssspa <- lapply(ssspa, function(x) x[order(as.numeric(sapply(x,
            "[[", 1)))])
    }
    else {
        ssspa <- sapply(ssspa, function(x) x[order(as.numeric(sapply(x,
            "[[", 1)))])
    }
    ssspa <- ssspa[order(sapply(ssspa, length))]
    ssspa <- ssspa[order(sapply(ssspa, function(x) {
        as.numeric(x[[1]][[1]])
    }))]
    ssspa_unlist <- unlist(ssspa, recursive = FALSE)
    spl <- con_split(ssspa_unlist, sapply(ssspa_unlist, "[[",
        3))
    tssspastb <- lapply(spl, function(x) ssextend(x, de))
    ssid <- sapply(tssspastb, function(x) {
        x[[2]][[1]]
    })

    nspecies <- length(sa[1][[1]][, 1])

    ssid_ordered <- .order_ssid_with_rightward_adjacency(
      gela, ssid, nspecies,
      w_adj = 20,         # increase if adjacency must be very strong
      w_energy = 2,
      w_mem = 1,
      w_gap = 5,
      adj_only = TRUE,
      adj_energy_scale = "median"
    )

    ord_idx <- match(ssid_ordered, ssid)
    ssid <- ssid[ord_idx]
    tssspastb <- tssspastb[ord_idx]



    if (length(ssid) > 1) {
        hh <- apply(matrix(stats::embed(ssid, 2)[, 2:1], ncol = 2),
            1, function(x) hamming_distance(id2bin(x[1], s),
                id2bin(x[2], s)))
        ml <- max(7, max(hh))
        posit <- Reduce(function(x, y) x + y, hh, accumulate = TRUE,
            init = ml + 1)
    }
    else {
        ml <- 7
        posit <- 8
    }
    cluster = makeCluster(threads)
    clusterCall(cluster, function(x) .libPaths(x), .libPaths())
    registerDoParallel(cluster)
    on.exit(stopCluster(cluster))
    filledsurf <- foreach(i = seq(length(de)), .packages = c("rELA",
        "tidyverse", "purrr", "foreach", "gtools"), .export = c("tssspastb",
        "refenv", "sa", "posit", "ml", "MinTippingPath", "MarginePath",
        "foldList", "replace_bin")) %dopar% {
        v <- lapply(tssspastb, function(x) x[i, ])
        v <- do.call(rbind, lapply(v, function(x) as.matrix(x[,
            ])))
        xde <- as.numeric(v[1, 1])
        basenv <- refenv
        basenv[is.na(basenv)] <- xde
        basenv <- t(basenv)
        sa2p <- sa2params(sa, c(basenv))
        he <- sa2p[[1]]
        je <- sa2p[[2]]
        ge <- sa2p[[3]]
        hge <- sa2p[[4]]
        pstb <- as.numeric(v[, 3])
        sid <- ssid[as.numeric(v[, 3]) == 1]
        zx <- posit[as.numeric(v[, 3]) == 1]
        scape <- unlist(map(seq_along(sid)[-length(sid)], function(j) {
            MinTippingPath(sid[j], sid[j + 1], hge, je, itr)
        }))
        if (length(scape) == 0) {
            scaledscp <- data.frame(z = zx, energy = MinTippingPath(sid[[1]],
                sid[[1]], hge, je, itr))
        }
        else {
            scaledscp <- data.frame(z = seq(zx[1], zx[length(zx)],
                length.out = length(scape) - 1), energy = scape[-length(scape)])
        }
        leftmargine <- MarginePath(ssid[1], zx[1] - 1, hge, je,
            itr)
        leftmargine <- data.frame(z = seq_along(leftmargine[-1]),
            energy = rev(leftmargine[-1]))
        rightmargine <- MarginePath(ssid[length(ssid)], ml +
            (posit[length(posit)] - zx[length(zx)]), hge, je,
            itr)
        rightmargine <- data.frame(z = seq(zx[length(zx)] + 1,
            posit[length(posit)] + 1 + ml), energy = rightmargine)
        rbind(leftmargine, scaledscp, rightmargine)
    }
    cc <- 0
    intpol <- lapply(filledsurf, function(df) {
        cc <<- cc + 1
        spline_fit <- smooth.spline(df$z, df$energy, spar = 0.2)
        pol <- predict(spline_fit, seq(df$z[1], df$z[length(df$z)],
            length.out = 100))$y
        dx <- seq(1, df$z[length(df$z)], length.out = 100)
        dex <- rep(de[cc], length(pol))
        data.frame(dx = dx, dex = dex, pol = pol)
    })
    surf3d <- lapply(intpol, function(df) {
        spline_fit <- smooth.spline(df$dx, df$pol, spar = 0.2)
        smoothed_pol <- predict(spline_fit, df$dx)$y
        filtered_pol <- stats::filter(smoothed_pol, filter = rep(1/5,
            5), sides = 2, circular = TRUE)
        data.frame(x = df$dx, y = df$dex, z = filtered_pol)
    })
    surf3dx <- dplyr::bind_rows(surf3d)
    l3d <- lapply(tssspastb, function(y) {
        si <- which(y$stb == 1)
        xe <- posit[which(ssid == y[1, 2])]
        xval <- surf3d[[1]][, 1]
        nearest_x <- which.min(abs(xval - xe))
        data.frame(t(sapply(si, function(x) surf3d[[x]][nearest_x,
            ])))
    })
    l3dx <- dplyr::bind_rows(lapply(seq_along(l3d), function(i) {
        df <- l3d[[i]]
        df$Index <- i
        df$ssid <- sapply(tssspastb, function(x) {
            unique(x[, 2])
        })[i]
        df <- df[, c(4, 5, 1, 2, 3)]
        return(df)
    }))
    return(list(surf3dx, l3dx))
}