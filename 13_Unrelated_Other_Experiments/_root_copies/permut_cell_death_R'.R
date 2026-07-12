}

gc()
# LT-HSC
# Sample data (replace with your actual data)
group_A <- c(346.56,411.36,584.16,416.16,440.16,232.32,456.96,497.76,831.6,895.2,658.8,796.8,518.4,696,564,693.6,778.8,1302,963.6,778.8,662.4,873.6,754.8,922.8,635.4,769.8,976.2,790.2,743.4,697.8,564.6,544.2)  # DOX group
group_B <- c(309.36, 221.28, 173.88, 227.88, 284.16, 246.84, 454.56, 246.24, 657.03, 217.95, 503.43, 449.43, 300.63, 387.03, 453.03, 459.03, 538.5, 860.1, 902.1, 598.5, 489.3, 576.9, 711.3, 744.9, 735, 603, 983.4, 773.4, 430.08, 597, 448.56, 446.04)     # 5FU group

# Calculate observed difference in sums
observed_diff <- sum(group_A) - sum(group_B)
cat("Observed difference in sums:", observed_diff, "\n")


# Number of permutations
num_permutations <- 100000

# Initialize vector to store permuted differences
perm_diffs <- numeric(num_permutations)

# Permutation test
# Permutation test
for (i in 1:num_permutations) {
  # Randomly shuffle group labels (permutation)
  permuted_groups <- sample(c(group_A, group_B))
  
  # Calculate permuted sums for each group
  perm_group_A <- permuted_groups[1:length(group_A)]
  perm_group_B <- permuted_groups[(length(group_A) + 1):length(permuted_groups)]
  
  # Calculate permuted difference in sums
  perm_diff <- sum(perm_group_A) - sum(perm_group_B)
  
  # Store permuted difference
  perm_diffs[i] <- perm_diff
}

# Calculate p-value
p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
cat("P-value:", p_value, "\n")

# Plotting histogram of permuted differences
hist(perm_diffs, breaks = 30, main = "Histogram of Permuted Differences",
     xlab = "Permuted Differences", ylab = "Frequency")

# Add observed differences as colored lines
abline(v = observed_diff, col = "red", lwd = 2)

# Add p-value to the plot
text(x = observed_diff, y = 10000, labels = paste("P-value =", round(p_value, 4)), col = "red")

# Show the plot

gc()
# Lin-_negative
# Sample data (replace with your actual data)
group_A <- c(225.012,244.212,268.212,228.612,279.012,199.932,263.412,269.412,339.012,399.6,349.2,343.2,306,328.8,294,312,410.4,621.6,481.2,400.8,330,367.2,344.4,390,373.8,433.8,456.6,448.2,384.6,400.2,315.96,337.8)  # DOX group
group_B <- c(216.972, 173.652, 175.092, 183.132, 235.812, 199.572, 246.612, 195.852, 372.3, 206.94, 291.9, 309.9, 242.7, 260.7, 267.9, 275.1, 408, 540, 566.4, 369.6, 269.76, 350.4, 372, 385.2, 401.4, 353.4, 437.4, 462.6, 358.2, 445.8, 359.4, 345)     # 5FU group

# Calculate observed difference in sums
observed_diff <- sum(group_A) - sum(group_B)
cat("Observed difference in sums:", observed_diff, "\n")


# Number of permutations
num_permutations <- 100000

# Initialize vector to store permuted differences
perm_diffs <- numeric(num_permutations)

# Permutation test
# Permutation test
for (i in 1:num_permutations) {
  # Randomly shuffle group labels (permutation)
  permuted_groups <- sample(c(group_A, group_B))
  
  # Calculate permuted sums for each group
  perm_group_A <- permuted_groups[1:length(group_A)]
  perm_group_B <- permuted_groups[(length(group_A) + 1):length(permuted_groups)]
  
  # Calculate permuted difference in sums
  perm_diff <- sum(perm_group_A) - sum(perm_group_B)
  
  # Store permuted difference
  perm_diffs[i] <- perm_diff
}

# Calculate p-value
p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
cat("P-value:", p_value, "\n")

# Plotting histogram of permuted differences
hist(perm_diffs, breaks = 30, main = "Histogram of Permuted Differences",
     xlab = "Permuted Differences", ylab = "Frequency")

# Add observed differences as colored lines
abline(v = observed_diff, col = "red", lwd = 2)

# Add p-value to the plot
text(x = observed_diff, y = 10000, labels = paste("P-value =", round(p_value, 4)), col = "red")

# Show the plot
