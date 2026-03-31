#!/bin/bash
# Pre-commit hook: blocks git commit if cargo fmt or clippy fail.
# Runs on PreToolUse when the Bash tool is about to execute a git commit.
set -euo pipefail

cd "$(git rev-parse --show-toplevel 2>/dev/null || echo "$CLAUDE_PROJECT_DIR")"

# Check formatting
if ! cargo fmt --all -- --check 2>/dev/null; then
    echo "BLOCKED: cargo fmt --check failed. Run 'cargo fmt --all' first." >&2
    exit 2
fi

# Check clippy
if ! cargo clippy --all-targets --all-features -- -D warnings 2>/dev/null; then
    echo "BLOCKED: cargo clippy found warnings. Fix them before committing." >&2
    exit 2
fi

exit 0
