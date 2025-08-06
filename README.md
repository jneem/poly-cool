# poly-cool

[![GitHub Actions CI status.](https://img.shields.io/github/actions/workflow/status/jneem/linesweeper/ci.yml?logo=github&label=CI)](https://github.com/jneem/poly-cool/actions)
[![Latest published version.](https://img.shields.io/crates/v/poly-cool.svg)](https://crates.io/crates/poly-cool)
[![Documentation build status.](https://img.shields.io/docsrs/poly-cool.svg)](https://docs.rs/poly-cool)

A pre-alpha crate for finding roots of low-degree polynomials.

This crate's goal is to provide well-tested and decently-optimized
methods for finding roots of polynomials. Currently, we have
implemented Blinn's method for cubics and Yuksel's algorithm
for polynomials of arbitrary degree. The "well-tested"
and "decently-optimized" parts are still aspirational.
