#!/bin/bash
# Stop hook: runs cargo fmt once per turn after Claude finishes responding.

# Consume stdin (required — Claude Code pipes JSON on stdin)
cat > /dev/null

cd "$CLAUDE_PROJECT_DIR"
cargo fmt --all >/dev/null 2>&1 || true
exit 0
