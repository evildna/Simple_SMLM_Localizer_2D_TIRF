#!/bin/bash
# Builds SMLM_Localization.jar against the ij.jar of a local Fiji (no Maven needed) and regenerates the Groovy script.
# Usage: ./build.sh [/path/to/Fiji.app]      (default: /Applications/Fiji.app)
# Compiles for Java 8, so the jar runs in both Java 8 and Java 21 Fiji.
# Environment overrides: JAVAC=/path/to/javac  IJ_JAR=/path/to/ij.jar
set -e
cd "$(dirname "$0")"
FIJI="${1:-/Applications/Fiji.app}"
if [ -z "$JAVAC" ]; then
  JAVAC="$(find "$FIJI/java" -type f -name javac 2>/dev/null | head -1)"; JAVAC="${JAVAC:-javac}"
fi
JAR="$(dirname "$(command -v "$JAVAC")")/jar"; [ -x "$JAR" ] || JAR=jar
CP="${IJ_JAR:-$(ls "$FIJI"/jars/ij-*.jar 2>/dev/null | head -1)}"
[ -f "$CP" ] || { echo "ij.jar not found; pass the Fiji.app path or set IJ_JAR" >&2; exit 1; }
SRC=src/main/java/no/uio/imaging/SMLM_Localization.java

# 1. Groovy script = Java source without the package line + a call to run()
{
  echo "// SMLM_Localization_Fiji.groovy - GENERATED from $SRC by build.sh; do not edit by hand."
  echo "// Install: copy to Fiji.app/plugins/Scripts/Plugins/ (menu Plugins > SMLM Localization Fiji), or open in the"
  echo "// Script Editor and press Run. Same code as the jar plugin; needs an open image."
  echo
  # @CompileStatic: without it Groovy 3/4 runs this numeric code ~25x slower than Groovy 2.4 or Java
  grep -v '^package no\.uio\.imaging;' "$SRC" | sed 's/^public class SMLM_Localization implements PlugIn {/@groovy.transform.CompileStatic\
public class SMLM_Localization implements PlugIn {/'
  echo
  echo 'new SMLM_Localization().run("")'
} > SMLM_Localization_Fiji.groovy

# 2. Jar
RELEASE="--release 8"; "$JAVAC" -version 2>&1 | grep -q ' 1\.8' && RELEASE="-source 1.8 -target 1.8"
rm -rf build && mkdir -p build/classes
"$JAVAC" $RELEASE -nowarn -encoding UTF-8 -cp "$CP" -d build/classes "$SRC"
cp src/main/resources/plugins.config LICENSE build/classes/
"$JAR" cf SMLM_Localization.jar -C build/classes .
echo "Built SMLM_Localization.jar and SMLM_Localization_Fiji.groovy"
