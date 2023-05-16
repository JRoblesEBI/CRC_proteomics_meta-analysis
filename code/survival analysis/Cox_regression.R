# Extract the survival time and status columns
time <- data$survival_time
status <- data$survival_status

# Create a new data frame to store the results
results <- data.frame(gene=character(), hazard_ratio=numeric(), p_value=numeric())

# Perform the Cox regression analysis for each gene


for (i in 2:(ncol(data)-2)) {
  gene <- colnames(data)[i]
  fit <- coxph(Surv(time, status) ~ data[,i])
  sum <- summary(fit)
  results <- rbind(results, data.frame(gene, sum$coefficients[2], sum$s[3]))
}

write.table(results, "results.txt", sep="\t", row.names=FALSE, col.names=TRUE)