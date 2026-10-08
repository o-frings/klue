# Alternative availability: wide columns av_1..av_J mask unavailable
# alternatives in the per-task softmax. All-ones availability must be
# byte-identical to no availability; a chosen-yet-unavailable alternative is an
# error. External agreement (vs Apollo with a real avail list) is checked in
# dev/benchmark_availability.R.

test_that("all-ones availability is byte-identical to no availability", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 150, T_tasks = 10, true_K = 2,
                     separation = 2.0, heterogeneity = 0.1, seed = 101, dgp = dgp)
  st <- klue_starts(d$database, 2, "kmeans", dgp = dgp)
  db1 <- d$database
  for (j in 1:3) db1[[paste0("av_", j)]] <- 1L
  f0 <- estimate_lcmnl(d$database, 2, start_betas = st$betas, dgp = dgp)
  f1 <- estimate_lcmnl(db1,        2, start_betas = st$betas, dgp = dgp)
  expect_identical(f1$LL, f0$LL)
  expect_identical(unname(f1$betas), unname(f0$betas))
  expect_identical(unname(f1$vcov), unname(f0$vcov))
})

test_that("masking an alternative in some tasks estimates cleanly", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 200, T_tasks = 10, true_K = 2,
                     separation = 2.0, heterogeneity = 0.1, seed = 5, dgp = dgp)
  db <- d$database
  for (j in 1:3) db[[paste0("av_", j)]] <- 1L
  # mark a non-chosen alternative unavailable in ~25% of tasks
  set.seed(5)
  for (r in seq_len(nrow(db))) if (runif(1) < 0.25) {
    rem <- setdiff(1:3, db$CHOICE[r]); db[[paste0("av_", rem[1])]][r] <- 0L
  }
  f <- klue_lcmnl(db, 2, dgp = dgp)
  expect_true(f$converged)
  expect_false(is.null(f$vcov))
  expect_true(all(eigen(f$vcov, only.values = TRUE)$values > 0))
  expect_true(all(sqrt(diag(f$vcov)) > 0))
})

test_that("klue_database keep_unavailable emits av columns and retains tasks", {
  set.seed(1); N <- 30; Tt <- 6; J <- 3
  rows <- list()
  for (id in 1:N) for (tk in 1:Tt) {
    av <- c(1L, 1L, 1L)
    if (runif(1) < 0.4) av[sample(1:2, 1)] <- 0L   # drop a non-reference alt
    ch <- sample(which(av == 1), 1)                # chosen must be available
    for (j in 1:J) rows[[length(rows) + 1]] <- data.frame(
      id = id, task = tk, alt = j, choice = as.integer(j == ch),
      X1 = rnorm(1), X2 = rnorm(1), price = runif(1), avail = av[j])
  }
  raw <- do.call(rbind, rows)
  common <- list(data = raw, format = "long", id_col = "id", task_col = "task",
                 alt_col = "alt", choice_col = "choice", choice_format = "indicator",
                 attribute_cols = c("X1", "X2"), price_col = "price",
                 avail_col = "avail", verbose = FALSE)
  d_drop <- do.call(klue_database, common)
  d_keep <- do.call(klue_database, c(common, list(keep_unavailable = TRUE)))

  expect_false(any(grepl("^av_", names(d_drop))))          # default drops, no av cols
  expect_true(all(c("av_1", "av_2", "av_3") %in% names(d_keep)))
  expect_gt(nrow(d_keep), nrow(d_drop))                    # keeps more tasks
  expect_true(all(vapply(seq_len(nrow(d_keep)),            # chosen always available
    function(r) d_keep[[paste0("av_", d_keep$CHOICE[r])]][r] == 1, logical(1))))
  f <- klue_lcmnl(d_keep, 1, dgp = klue_dgp(2, 3))         # estimator runs on it
  expect_true(f$converged)
})

test_that("a chosen-yet-unavailable alternative is rejected", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 80, T_tasks = 6, true_K = 1,
                     separation = 1.0, heterogeneity = 0.1, seed = 9, dgp = dgp)
  db <- d$database
  for (j in 1:3) db[[paste0("av_", j)]] <- 1L
  db[[paste0("av_", db$CHOICE[1])]][1] <- 0L    # chosen alt marked unavailable
  expect_error(estimate_lcmnl(db, 1, dgp = dgp), "unavailable")
})
