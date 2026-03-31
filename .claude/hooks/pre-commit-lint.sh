#!/bin/bash
# PreToolUse hook: blocks git commit if cargo fmt or clippy fail.
# Only activates on git commit commands — skips everything else.

# Consume stdin (required — Claude Code pipes JSON on stdin)
INPUT=$(cat)

# Only run on git commit commands
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty')
case "$COMMAND" in
    git\ commit*) ;;
    *) exit 0 ;;
esac

cd "$CLAUDE_PROJECT_DIR"

# Check formatting
if ! cargo fmt --all -- --check >/dev/null 2>&1; then
    echo "cargo fmt --check failed. Run 'cargo fmt --all' first." >&2
    exit 2
fi

# Check clippy
if ! cargo clippy --all-targets --all-features -- -D warnings >/dev/null 2>&1; then
    echo "cargo clippy found warnings. Fix them before committing." >&2
    exit 2
fi

exit 0
