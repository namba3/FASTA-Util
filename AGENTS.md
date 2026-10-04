# Repository guidance

## Project

- This repository contains the `fasta-util` Rust CLI and the separate `generate_random_data` helper crate.
- The main crate uses nightly-only Rust features. Use the nightly toolchain for Cargo commands.
- Read `README.md` for the user-facing commands and examples before changing CLI behavior.

## Changes

- Keep FASTA headers (lines beginning with `>`) distinct from sequence data. Sequence validation currently accepts the uppercase symbols in `src/nucleic_acid.rs`, including `-`; coordinate or alphabet changes should be deliberate and covered by tests.
- The `slice` range is zero-based and follows Rust range syntax: `a..b` excludes `b`, while `a..=b` includes it. Preserve this behavior unless the requested change explicitly changes it.
- Keep output streaming for large FASTA files. The main file reader memory-maps input, so review lifetime and safety assumptions carefully when changing `src/lib.rs`.
- Keep the random-data helper independently runnable through its own manifest at `generate_random_data/Cargo.toml`.
- Avoid adding dependencies when the standard library or existing crates are sufficient.

## Useful commands

```sh
cargo +nightly fmt --check
cargo +nightly test
cargo +nightly build --release
cargo +nightly run --manifest-path=generate_random_data/Cargo.toml -- 10000
```

The root crate contains nightly benchmark code; run benchmarks with `cargo +nightly bench` when performance work requires it. Do not treat benchmark timings as stable across machines or runs.
