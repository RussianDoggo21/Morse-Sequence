#!/bin/bash

# Stop on any error
set -e

# Ensure VIRTUAL_ENV is set
if [ -z "$VIRTUAL_ENV" ]; then
  echo "❌ Error: VIRTUAL_ENV is not set. Please activate your virtual environment first."
  exit 1
fi

# === Step 1: Clean previous installation ===
echo "🧹 Removing previous build directory..."
rm -rf build

echo "🧹 Removing previously installed morse_sequence from venv..."
rm -rf "$VIRTUAL_ENV/lib/python3.13/site-packages/morse_sequence"
rm -f "$VIRTUAL_ENV/lib/python3.13/site-packages/_core"*.so

# === Step 2: Set up Meson ===
echo "♻️ Setting up Meson with venv prefix..."
meson setup build --wipe --prefix="$VIRTUAL_ENV"

# === Step 3: Compile ===
echo "🛠️ Compiling..."
meson compile -C build

# === Step 4: Install ===
echo "📦 Installing into venv..."
meson install -C build

# === Step 5: Verify installation ===
SITEPKG=$(python3 -c "import site; print(site.getsitepackages()[0])")
echo "📌 Python site-packages detected at: $SITEPKG"

echo "✅ Done. morse_sequence and _core.so are installed in your venv."
