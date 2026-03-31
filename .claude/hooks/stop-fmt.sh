#!/bin/bash
# Stop hook: runs cargo fmt once per turn after Claude finishes responding.
# Follows the official Claude Code hooks pattern from code.claude.com/docs/en/hooks

# Consume stdin (required — Claude Code pipes JSON on stdin)
INPUT=$(cat)

# Prevent infinite loop: if this is a re-fire, allow stop
if echo "$INPUT" | jq -re '.stop_hook_active' >/dev/null 2>&1; then
    exit 0
fi

cd "$CLAUDE_PROJECT_DIR"
cargo fmt --all >/dev/null 2>&1 || true
exit 0
