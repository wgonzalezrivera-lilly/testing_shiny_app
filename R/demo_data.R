make_demo_finemap_data <- function(n = 5000, seed = 22) {
	set.seed(seed)
	position <- sort(sample(seq_len(51e6), n))
	pip <- rbeta(n, shape1 = 0.35, shape2 = 8)
	credible_signal <- sample(seq_len(n), 18)
	pip[credible_signal] <- sort(runif(length(credible_signal), 0.65, 0.99), decreasing = TRUE)
	data.frame(
		CHR = 22L,
		POS = position,
		SNP = sprintf("rs%07d", seq_len(n)),
		PIP = pip,
		stringsAsFactors = FALSE
	)
}
