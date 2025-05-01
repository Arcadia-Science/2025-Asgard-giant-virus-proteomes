#!/usr/bin/env bash

# Wrapper script to launch InterProScan 5, handling configuration,
# environment setup, and Java validation.

# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# Prevent errors in a pipeline from being masked.
set -euo pipefail

# --- Functions ---

log_info() {
    echo "INFO: $1"
}

log_error() {
    # Print error messages to stderr
    echo "ERROR: $1" >&2
}

# Extracts the value of 'bin.directory' from a Java properties file.
# Expects lines like 'bin.directory=/path/to/bin'
# Args:
#   $1: Path to the properties file.
# Outputs:
#   The extracted path, or empty string if not found.
get_bin_directory() {
    local properties_file="$1"
    if [ ! -f "$properties_file" ]; then
        log_error "Properties file not found: '$properties_file'"
        return 1 # Indicate failure
    fi
    # Use grep and cut for potentially more robust parsing than sed/awk
    grep '^bin.directory=' "$properties_file" | cut -d'=' -f2-
}

# --- Initialization & Path Handling ---

# Record the directory from which the script was invoked
USER_INVOKE_DIR=$PWD

# Determine the script's absolute directory, resolving symlinks
SCRIPT_SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT_SOURCE" ]; do # Resolve $SCRIPT_SOURCE until it's no longer a symlink
  SCRIPT_DIR="$( cd -P "$( dirname "$SCRIPT_SOURCE" )" &> /dev/null && pwd )"
  SCRIPT_SOURCE="$(readlink "$SCRIPT_SOURCE")"
  # If $SCRIPT_SOURCE was a relative symlink, resolve it relative to the symlink's directory
  [[ $SCRIPT_SOURCE != /* ]] && SCRIPT_SOURCE="$SCRIPT_DIR/$SCRIPT_SOURCE"
done
SCRIPT_DIR="$( cd -P "$( dirname "$SCRIPT_SOURCE" )" &> /dev/null && pwd )"
log_info "Script directory identified as: ${SCRIPT_DIR}"

# --- Configuration Loading ---

DEFAULT_PROPERTIES_FILE="${SCRIPT_DIR}/interproscan.properties"
log_info "Looking for default properties file: ${DEFAULT_PROPERTIES_FILE}"

BIN_DIR=$(get_bin_directory "${DEFAULT_PROPERTIES_FILE}")
if [ -z "$BIN_DIR" ]; then
    log_error "Could not determine 'bin.directory' from ${DEFAULT_PROPERTIES_FILE}"
    exit 1
fi
log_info "Default bin directory: ${BIN_DIR}"

# Check for custom configuration file via environment variable
if [ -n "${INTERPROSCAN_CONF:-}" ] && [ -f "$INTERPROSCAN_CONF" ]; then
    # Use :- to avoid error if INTERPROSCAN_CONF is unset
    log_info "Custom configuration file specified: ${INTERPROSCAN_CONF}"
    PROPERTIES_FILE="$(readlink -m "$INTERPROSCAN_CONF")" # Get canonical path
    log_info "Using properties file: ${PROPERTIES_FILE}"

    CUSTOM_BIN_DIR=$(get_bin_directory "${PROPERTIES_FILE}")
    if [ -z "$CUSTOM_BIN_DIR" ]; then
        log_error "Could not determine 'bin.directory' from custom config ${PROPERTIES_FILE}"
        exit 1
    fi
    BIN_DIR="$CUSTOM_BIN_DIR" # Override BIN_DIR
    log_info "Overridden bin directory from custom config: ${BIN_DIR}"
    JAVA_PROPERTY_ARG="-Dsystem.interproscan.properties=${PROPERTIES_FILE}"
else
    log_info "Using default properties file: ${DEFAULT_PROPERTIES_FILE}"
    PROPERTIES_FILE="${DEFAULT_PROPERTIES_FILE}" # Keep track of which properties file is actually used
    JAVA_PROPERTY_ARG=""
fi

# --- Environment Setup for Embedded Tools (e.g., getorf) ---

# Check if the determined binary directory and subdirectories exist
NUCLEOTIDE_DIR="${BIN_DIR}/nucleotide"
if [ ! -d "$NUCLEOTIDE_DIR" ]; then
     log_error "Required EMBOSS directory not found: ${NUCLEOTIDE_DIR}"
     log_error "Check the 'bin.directory' setting in '${PROPERTIES_FILE}'"
     exit 1
fi
export EMBOSS_ACDROOT="$NUCLEOTIDE_DIR"
export EMBOSS_DATA="$NUCLEOTIDE_DIR"
log_info "EMBOSS environment variables set using: ${NUCLEOTIDE_DIR}"

# --- Java Validation ---

# Find Java executable
JAVA_EXEC=$(type -p java)
if [ -z "$JAVA_EXEC" ]; then
    log_error "Java executable not found in PATH."
    log_error "Please install Java 11 (or later) and ensure it's in your PATH,"
    log_error "or specify the path in '${PROPERTIES_FILE}' (e.g., java.command=/path/to/java)."
    exit 1
fi
log_info "Using Java executable: ${JAVA_EXEC}"

# Check Java version
log_info "Checking Java version..."
# Execute java -version, redirect stderr to stdout for parsing
JAVA_VERSION_OUTPUT=$("$JAVA_EXEC" -Xms32M -Xmx32M -version 2>&1)
JAVA_VERSION_EXIT_CODE=$?
if [ $JAVA_VERSION_EXIT_CODE -ne 0 ]; then
    log_error "Failed to execute '$JAVA_EXEC -version' (Exit code: ${JAVA_VERSION_EXIT_CODE}). Output:"
    echo "$JAVA_VERSION_OUTPUT" >&2 # Show output on error
    exit 1
fi

# Parse version string (handles formats like "11.0.1", "17", etc.)
# Looks for 'version "X.Y.Z"' pattern
JAVA_VERSION_STRING=$(echo "$JAVA_VERSION_OUTPUT" | grep -oP 'version "?\K[0-9]+\.[0-9]+\.[0-9]+' | head -n 1)
# Fallback for simpler formats like "11" or "17" if the above fails
if [ -z "$JAVA_VERSION_STRING" ]; then
    JAVA_VERSION_STRING=$(echo "$JAVA_VERSION_OUTPUT" | grep -oP 'version "?\K[0-9]+' | head -n 1)
fi

if [ -z "$JAVA_VERSION_STRING" ]; then
    log_error "Could not parse Java version string from output:"
    echo "$JAVA_VERSION_OUTPUT" >&2
    exit 1
fi

# Extract major version number
JAVA_MAJOR_VERSION="${JAVA_VERSION_STRING%%.*}"

# Check if major version is at least 11
# Using arithmetic comparison
if [ "${JAVA_MAJOR_VERSION}" -lt 11 ]; then
    log_error "Java version 11 or later is required to run InterProScan."
    log_error "Detected version: ${JAVA_VERSION_STRING} (Major version: ${JAVA_MAJOR_VERSION})"
    log_error "Please install or configure the correct Java version."
    exit 1
fi
log_info "Java version check passed: ${JAVA_VERSION_STRING}"

# --- Execution ---

# Check if the JAR file exists
JAR_FILE="${SCRIPT_DIR}/interproscan-5.jar"
if [ ! -f "$JAR_FILE" ]; then
    log_error "InterProScan JAR file not found: ${JAR_FILE}"
    exit 1
fi

# Define memory settings - consider making these configurable
# Example: export IPRSCAN_JAVA_MEM_OPTS="-Xms2028M -Xmx9216M"
JAVA_MEM_OPTS="-Xms2028M -Xmx9216M"
# Define GC threads - consider making configurable
JAVA_GC_OPTS="-XX:ParallelGCThreads=8"

# The '-u $USER_INVOKE_DIR' argument tells InterProScan to use the original
# invocation directory for relative path resolution for input/output files.
log_info "Launching InterProScan 5..."
log_info "Command: ${JAVA_EXEC} ${JAVA_GC_OPTS} ${JAVA_MEM_OPTS} ${JAVA_PROPERTY_ARG} -jar ${JAR_FILE} $* -u ${USER_INVOKE_DIR}"

"$JAVA_EXEC" \
 ${JAVA_GC_OPTS} \
 ${JAVA_MEM_OPTS} \
 ${JAVA_PROPERTY_ARG} \
 -jar "${JAR_FILE}" "$@" -u "$USER_INVOKE_DIR" # Pass script arguments, ensure USER_INVOKE_DIR is quoted

# Capture and report the exit code of the java command
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    log_error "InterProScan exited with status ${EXIT_CODE}"
else
    log_info "InterProScan finished successfully (Exit code: ${EXIT_CODE})."
fi

exit $EXIT_CODE
