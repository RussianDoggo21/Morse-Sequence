#!/bin/bash

# Stop on any error
set -e

# Check that $VIRTUAL_ENV is set
if [ -z "$VIRTUAL_ENV" ]; then
  echo "Error: VIRTUAL_ENV is not set. Please activate your virtual environment first."
  exit 1
fi

# Remove existing build directory
echo "Removing build directory..."
rm -rf build

# Set up Meson build
echo "Setting up Meson..."
meson setup build --prefix="$VIRTUAL_ENV"

# Compile the project
echo "Compiling..."
meson compile -C build

# Install the project
echo "Installing..."
meson install -C build

echo "✅ Done."
