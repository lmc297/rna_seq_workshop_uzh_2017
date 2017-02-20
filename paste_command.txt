ls *.cov | paste -sd'\t' - >> coverage_combined.cov | paste -d"\t" *sense.cov >> coverage_combined.cov
