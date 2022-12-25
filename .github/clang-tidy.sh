#!/usr/bin/env bash
clang-tidy --checks=-*,bugprone-*,performance-*,readability-*,-readability-magic-numbers,misc-*,clang-analyzer-*,-clang-analyzer-cplusplus* --warnings-as-errors=* "$@" --
