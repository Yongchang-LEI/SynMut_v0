#!/bin/bash

# Pull the latest changes from GitHub
echo "Pulling latest changes from GitHub..."
git pull

# Install the package
echo "Installing SynMutDev..."
R CMD INSTALL .

echo "Done! You can now use library('SynMutDev') in R."
