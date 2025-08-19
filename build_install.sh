#!/bin/bash

# Stop on any error
set -e

# Ensure VIRTUAL_ENV is set
if [ -z "$VIRTUAL_ENV" ]; then
  echo "❌ Error: VIRTUAL_ENV is not set. Please activate your virtual environment first."
  exit 1
fi

# Remove previous build
echo "🧹 Removing build directory..."
rm -rf build

# Set up Meson for full rebuild with prefix pointing to the venv
echo "♻️ Setting up Meson with venv prefix..."
meson setup build --wipe --prefix="$VIRTUAL_ENV"

# Compile
echo "🛠️ Compiling..."
meson compile -C build

# Install into the venv
echo "📦 Installing into venv..."
meson install -C build

# Verify installation
SITEPKG=$(python3 -c "import site; print(site.getsitepackages()[0])")
echo "📌 Python site-packages detected at: $SITEPKG"

echo "✅ Done. morse_sequence and _core.so are installed in your venv."
