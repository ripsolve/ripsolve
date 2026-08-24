#!/usr/bin/env bash
# Build the Python extension module.
#
# PyO3 produces a cdylib that CPython can import once it is named for the module.
# maturin would do this (and build a wheel); this script exists so the extension
# can be built and tested with nothing but cargo and a Python with headers.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cargo build --release --manifest-path "$here/../Cargo.toml" -p ripsolve-py

# Linux names it lib<name>.so; Python wants <name>.so on the import path.
case "$(uname -s)" in
    Darwin) built="libripsolve.dylib" ;;
    *)      built="libripsolve.so" ;;
esac
cp "$here/../target/release/$built" "$here/ripsolve.so"
echo "built $here/ripsolve.so"
