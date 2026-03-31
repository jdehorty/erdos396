#!/bin/bash
# Stop hook: runs cargo fmt once per turn, after Claude finishes responding.
# This avoids the PostToolUse cache-staleness problem where formatting a file
# mid-turn causes Claude to see stale content on its next read.
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo "$CLAUDE_PROJECT_DIR")"
cargo fmt --all 2>/dev/null || true
exit 0
