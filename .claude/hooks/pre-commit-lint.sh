#!/bin/bash
# PreToolUse hook: blocks git commit if cargo fmt or clippy fail.
# Only activates on git commit commands — skips everything else.

# Consume stdin (required — Claude Code pipes JSON on stdin)
INPUT=$(cat)

# Require jq for command filtering
if ! command -v jq >/dev/null 2>&1; then
    echo "pre-commit-lint.sh: jq is required but not installed" >&2
    exit 2
fi

# Only run on git commit commands
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty')
case "$COMMAND" in
    git\ commit*) ;;
    *) exit 0 ;;
esac

cd "$CLAUDE_PROJECT_DIR"

# Check formatting (forward diagnostics on failure)
FMT_OUTPUT=$(cargo fmt --all -- --check 2>&1)
if [ $? -ne 0 ]; then
    echo "cargo fmt --check failed:" >&2
    echo "$FMT_OUTPUT" >&2
    exit 2
fi

# Check clippy (forward diagnostics on failure)
CLIPPY_OUTPUT=$(cargo clippy --all-targets --all-features -- -D warnings 2>&1)
if [ $? -ne 0 ]; then
    echo "cargo clippy found warnings:" >&2
    echo "$CLIPPY_OUTPUT" >&2
    exit 2
fi

exit 0
