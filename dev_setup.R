# Setup script for SynMut development
# This script installs devtools if needed and loads the package in development mode.

if (!requireNamespace("devtools", quietly = TRUE)) {
    message("Installing devtools...")
    install.packages("devtools", repos = "https://cloud.r-project.org")
}

message("Loading SynMut in development mode...")
devtools::load_all(".")

message("SynMut loaded! You can now use the package functions.")
message("Any changes you make to the code will be available after running devtools::load_all('.') again.")
