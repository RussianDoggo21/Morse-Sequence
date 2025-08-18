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

# Set up Meson
echo "♻️ Forcing full rebuild..."
meson setup build --wipe --prefix=/usr/local

# Compile
echo "🛠️ Compiling..."
meson compile -C build

# Install into a temporary directory inside the venv
echo "📦 Installing to venv using --destdir..."
meson install -C build --destdir="$VIRTUAL_ENV"

# Find Python's site-packages in venv
SITEPKG=$(python3 -c "import site; print(site.getsitepackages()[0])")

# Copy installed files from venv/usr/local/... to site-packages
echo "📁 Copying installed files to site-packages..."
cp -r "$VIRTUAL_ENV/usr/local/lib/python3.13/dist-packages/morse_sequence" "$SITEPKG" 2>/dev/null || true
cp "$VIRTUAL_ENV/usr/local/lib/python3.13/dist-packages/"_core*.so "$SITEPKG/morse_sequence/" 2>/dev/null || true

echo "✅ Done. The package is now installed in your virtual environment."
