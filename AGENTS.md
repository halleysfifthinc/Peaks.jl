# AGENTS

This file provides guidance to LLM agents when working with code in this repository.

Add under a ## Julia Development section at the top level\n\n

## Development style

- When debugging or reviewing code, always re-read the current state of files before making suggestions or changes to make sure that diffs/etc apply cleanly to the current state of things - previous context may be outdated.
- Never repeat approaches that failed in previous attempts within the same problem domain - ask user about prior attempts first if this is a new session

### Tool use

- Always avoid changing the directory with `cd`, instead prefer using tool options to handle different directories, e.g. cmake has `-B`/`--build` options, and julia has `--project[={<dir>|@temp|@.|@script[<rel>]}]`
- Never give an absolute path to `--project`

## Development Commands

### Package Management
```bash
# Activate package environment in project root
julia --project=@.

# activate environment of subdirectory
julia --project=subdir

# Install/update dependencies
julia --project=@. -e "using Pkg; Pkg.instantiate()"

# Check package status
julia --project=@. -e "using Pkg; Pkg.status()"
```

Whenever changing (adding, removing, etc) project dependencies, try to use the Pkg.jl in julia first, before resorting to editing the Project.toml file directly.

### Testing
```bash
# Run all tests (includes Aqua, JET, ExplicitImports quality checks)
julia --project=@. -e "using Pkg; Pkg.test()" # generate coverage statistics with the `coverage` keyword argument
```

### Documentation
```bash
# Build documentation locally
julia --project=docs docs/make.jl

# Documentation is built using Documenter.jl with auto-deployment to gh-pages
```

### Benchmarking
```bash
# Run benchmarks (benchmark environment available)
julia --project=benchmark benchmark/benchmarks.jl
```

## Architecture Overview

### Core Module Structure
- **Main module**: `src/Peaks.jl` - Exports all public functions and includes all submodules
- **Peak finding algorithms**: `src/minmax.jl` - Core extrema detection
- **Simple peak detection**: `src/simple.jl` - Fast, basic peak finding functions with SIMD optimization
- **Peak characteristics**:
  - `src/peakprom.jl` - Peak prominence calculations
  - `src/peakwidth.jl` - Peak width measurements (including FWHM)
  - `src/peakheight.jl` - Peak height calculations
- **Utilities**: `src/utils.jl` - Helper functions and utilities
- **Visualization**: `src/plot.jl` - RecipesBase integration for plotting

### Key Dependencies
- **SIMD.jl**: Used for vectorized operations in core algorithms for performance
- **RecipesBase.jl**: Provides plotting integration with Plots.jl ecosystem

### Algorithm Design
The package uses a window-based approach for peak detection with support for:
- Configurable window sizes for peak isolation
- Strict/non-strict handling of NaN/missing values via `strict` parameter
- Plateau detection for flat peaks
- SIMD-optimized computations for performance

### Function Naming Convention
- `arg*` functions return indices (e.g., `argmaxima`, `argminima`)
- `simple*` functions provide fast basic detection (e.g., `simplemaxima`, `simpleminima`)
- `find*` functions offer comprehensive peak finding with filtering options
- `peak*` functions calculate peak characteristics (prominence, width, height)
- `is*` functions provide boolean tests (e.g., `ismaxima`, `isminima`, `isplateau`)

### Error Handling
- Custom error hints registered in `__init__()` to guide users when `simplemaxima`/`simpleminima` are used with missing values
- Domain errors for invalid window sizes
- Graceful handling of edge cases with NaN/missing data

### Testing Strategy
- Comprehensive test suite covering all major functions
- Code quality checks using Aqua.jl, JET.jl, and ExplicitImports.jl
- SIMD-specific tests for performance-critical code paths
- Plotting integration tests
- Manual test files for interactive development

### Documentation Structure
- Comprehensive docs using Documenter.jl with examples
- Interactive examples with plots for visualization
- Glossary of peak-finding terminology
- Benchmark comparisons and performance analysis
