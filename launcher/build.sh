#!/bin/sh
# Build the one-click launchers from any OS with Go >= 1.22 installed:
#   dist/kerncount.exe       Windows, 64-bit
#   dist/kerncount-mac.zip   kerncount.app for Mac (Apple silicon)
#   sh launcher/build.sh
# The exe embeds the current kerncount source (pyproject.toml + src/kerncount, no model weights or tests).
set -e
cd "$(dirname "$0")/.."
python3 - <<'EOF'
import zipfile, pathlib
root = pathlib.Path(".")
files = [root / "pyproject.toml", root / "README.md"] + [
    p for p in (root / "src").rglob("*")
    if p.is_file() and "__pycache__" not in p.parts and "models" not in p.parts and ".egg-info" not in str(p)]
with zipfile.ZipFile("launcher/app.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(files):
        z.write(p, p.as_posix())
print(f"launcher/app.zip: {len(files)} files")
EOF
mkdir -p dist
cd launcher
GOOS=windows GOARCH=amd64 CGO_ENABLED=0 go build -trimpath -ldflags "-s -w" -o ../dist/kerncount.exe .
echo "built dist/kerncount.exe"

# Mac: a minimal .app bundle around the arm64 binary
APP=../dist/kerncount.app
rm -rf "$APP" ../dist/kerncount-mac.zip
mkdir -p "$APP/Contents/MacOS"
GOOS=darwin GOARCH=arm64 CGO_ENABLED=0 go build -trimpath -ldflags "-s -w" -o "$APP/Contents/MacOS/kerncount" .
VERSION=$(sed -n 's/^version = "\(.*\)"/\1/p' ../pyproject.toml)
cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>kerncount</string>
  <key>CFBundleDisplayName</key><string>kerncount</string>
  <key>CFBundleIdentifier</key><string>dev.kerncount.launcher</string>
  <key>CFBundleExecutable</key><string>kerncount</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleShortVersionString</key><string>${VERSION}</string>
  <key>CFBundleVersion</key><string>${VERSION}</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>LSArchitecturePriority</key><array><string>arm64</string></array>
</dict></plist>
PLIST
if command -v codesign >/dev/null 2>&1; then
  xattr -cr "$APP"                          # no extended attributes: they become ._ files in the zip
  codesign --force --deep --sign - "$APP"   # ad-hoc signature (Apple silicon refuses unsigned code)
fi
cd ../dist
if command -v ditto >/dev/null 2>&1; then
  # no resource forks / attributes, so any unzip tool (not only Finder) gives a validly signed app
  ditto -c -k --norsrc --noextattr --noacl --keepParent kerncount.app kerncount-mac.zip
else
  zip -qr kerncount-mac.zip kerncount.app
fi
echo "built dist/kerncount.app and dist/kerncount-mac.zip"
