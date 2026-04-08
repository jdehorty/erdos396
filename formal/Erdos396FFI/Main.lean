import Erdos396FFI.Basic

/-!
# Erdős 396 — Lean Oracle CLI

Batch-mode executable for differential random testing (DRT).

## Protocol

One command per line on stdin, one result per line on stdout.

| Command    | Output            | Description                         |
|------------|-------------------|-------------------------------------|
| `W k n`    | `0` or `1`        | Is `n` a `k`-witness?              |
| `G n`      | `0` or `1`        | Is `n` in the Governor Set?         |
| `VC n p`   | integer           | `v_p(C(2n, n))` via Legendre       |
| `VF n p`   | integer           | `v_p(n!)` via Legendre             |
| `VK n p`   | integer           | `v_p(C(2n, n))` via Kummer         |
| `VP m p`   | integer           | `v_p(m)`                           |
| `Q`        | *(exits)*         | Quit                                |

All other lines produce `ERR`.
-/

open Erdos396FFI

def processLine (primes : Array Nat) (line : String) : String :=
  let parts := line.trim.splitOn " "
  match parts with
  | ["W", kStr, nStr] =>
    match kStr.toNat?, nStr.toNat? with
    | some k, some n => if checkWitness k n primes then "1" else "0"
    | _, _ => "ERR"
  | ["G", nStr] =>
    match nStr.toNat? with
    | some n => if isGovernor n primes then "1" else "0"
    | none => "ERR"
  | ["VC", nStr, pStr] =>
    match nStr.toNat?, pStr.toNat? with
    | some n, some p => toString (vpCentralBinom n p)
    | _, _ => "ERR"
  | ["VF", nStr, pStr] =>
    match nStr.toNat?, pStr.toNat? with
    | some n, some p => toString (vpFactorial n p)
    | _, _ => "ERR"
  | ["VK", nStr, pStr] =>
    match nStr.toNat?, pStr.toNat? with
    | some n, some p => toString (vpCentralBinomKummer n p)
    | _, _ => "ERR"
  | ["VP", mStr, pStr] =>
    match mStr.toNat?, pStr.toNat? with
    | some m, some p => toString (vp m p)
    | _, _ => "ERR"
  | _ => "ERR"

def parseSieveLimit (args : List String) : Option Nat :=
  match args with
  | ["--batch"] => some 5000000
  | ["--batch", limitStr] => some (limitStr.toNat?.getD 5000000)
  | _ => none

def main (args : List String) : IO UInt32 := do
  match parseSieveLimit args with
  | none =>
    IO.eprintln "Usage: erdos396ffi --batch [sieve_limit]"
    IO.eprintln "Reads commands from stdin, writes results to stdout."
    IO.eprintln "Commands: W k n | G n | VC n p | VF n p | VK n p | VP m p | Q"
    return 1
  | some sieveLimit =>
    let primes := sievePrimes sieveLimit
    let stderr ← IO.getStderr
    stderr.putStrLn s!"READY {primes.size}"
    stderr.flush

    let stdin ← IO.getStdin
    let stdout ← IO.getStdout

    let mut done := false
    while !done do
      let line ← stdin.getLine
      let trimmed := line.trim
      if trimmed == "Q" || trimmed.isEmpty then
        done := true
      else
        let result := processLine primes trimmed
        stdout.putStrLn result
        stdout.flush

    return 0
