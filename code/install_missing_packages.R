###############################################################################
# File: install_missing_packages.R
# Purpose: Install missing spatial packages required for the analysis
# Author: Claude Code
# Date: 2025-10-15
###############################################################################

message("Installing missing spatial packages...")

# List of missing packages that need to be installed
missing_packages <- c("rnaturalearth", "rnaturalearthdata", "ggspatial")

# Install each package
for (pkg in missing_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s...", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(sprintf("%s is already installed.", pkg))
  }
}

# Verify installation
message("\nVerifying installations:")
for (pkg in missing_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("✓ %s installed successfully", pkg))
  } else {
    warning(sprintf("✗ %s installation failed", pkg))
  }
}

message("\nNow updating renv snapshot...")
if (requireNamespace("renv", quietly = TRUE)) {
  renv::snapshot(prompt = FALSE)
  message("✓ renv snapshot updated")
} else {
  message("Note: renv not available to update snapshot")
}

message("\nDone! You can now run source('code/0_libraries.R')")
