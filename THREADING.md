# Thread Safety in Open Babel

This document summarizes the multi-threading work done on Open Babel 3.x and the current threading guarantees.

## Overview

The primary goal was to make the **forcefield energy minimization pipeline** safe for multi-threaded use, eliminating data races and crashes when multiple threads process different molecules concurrently.

## What Is Now Thread-Safe

### Core Infrastructure
- **Plugin/format registries** — protected by `std::recursive_mutex`
- **Error log** (`obErrorLog`) — protected by `std::mutex`
- **Locale handling** (`obLocale`) — thread-local counters on macOS, mutex-protected fallback
- **Type table** (`OBTypeTable` / `ttab`) — `std::recursive_mutex` + atomic `Translate` overload
- **Chain parser** (`chainsparser`) — mutex-protected
- **Ring typer** (`ringtyper`) — `thread_local`
- **Plugin loading** (`AllPluginsLoaded`) — `std::atomic<int>` with double-checked locking

### OBConversion
- **Option parameter registration** — mutex-protected static maps
- **Static string return buffers** — converted to `thread_local static` in:
  - `OBMol::GetTitle()`
  - `OBConversion::GetNextFormat()`
  - `OBMol::ClassDescription()`
  - Plugin description functions (descriptors, fingerprints, ops)

### Forcefields
All forcefields now implement **copy-based `MakeNewInstance()`**:
- `MMFF94`, `UFF`, `Ghemical`, `GAFF`, `MM2`
- Instead of re-parsing parameter files from disk (~2.6× overhead), instances are cloned by copying in-memory parameter vectors (~0.19 µs per call)
- Per-instance `Setup()` means each thread can minimize a different molecule safely

### Library Callers Updated
- `OpEnergy`, `OpMinimize` (`src/ops/forcefield.cpp`)
- `OpConformer` (`src/ops/conformer.cpp`)
- `TinkerFormat` MMFF94 atom type output (`src/formats/tinkerformat.cpp`)
- `Gen3D` (`src/ops/gen3d.cpp`)
- `MMFF94Charges` (`src/charges/mmff94.cpp`)
- `EQEqCharges` (`src/charges/eqeq.cpp`) — param file loading now atomic

## Performance (Regular Build, No Sanitizers)

### MMFF94 Minimization Throughput
Test: 11 MOL2 molecules, 10 SD steps + 10 CG steps, repeated 10×

| Threads | Time (s) | Molecules | mol/sec | Speedup | Efficiency |
|---------|----------|-----------|---------|---------|------------|
| 1       | 0.836    | 110       | 131.5   | 1.00    | 1.00       |
| 2       | 0.447    | 110       | 246.3   | 1.87    | 0.94       |
| 4       | 0.265    | 110       | 415.4   | 3.16    | 0.79       |
| 8       | 0.299    | 110       | 367.7   | 2.80    | 0.35       |

*Platform: macOS arm64, M1 Max (8 performance + 2 efficiency cores). The 8-thread dropoff is due to saturating the available performance cores and memory bandwidth — not TSan.*

### MakeNewInstance Overhead
- **Copy-based `MakeNewInstance`:** ~0.19 µs per call
- **Full pipeline per molecule:** ~7.6 ms
- **Overhead:** <0.01% — effectively indistinguishable from prototype reuse

### Comparison to Original Code
| State | Approach | Relative Speed |
|-------|----------|----------------|
| Original (buggy) | Prototype reuse directly | 1.00× (crashed on multi-molecule) |
| Intermediate | `MakeNewInstance` re-parsing from disk | ~0.38× (2.6× slower) |
| **Final** | **Copy-based `MakeNewInstance`** | **~1.00×** |

## Testing

### Test Suite
- `test_threading_accuracy_0..5` — 6-part accuracy test covering:
  1. Bare threads + OBConversion
  2. AddHydrogens
  3. Canonical SMILES write
  4. SMILES roundtrip
  5. MMFF94 minimization
  6. Multi-format conversion
- `test_threading_performancetest` — throughput and scaling measurement
- Full CTest suite: **237/237 passing**

### TSan Validation (Separate Build)
A dedicated build with `-fsanitize=thread` was used only for race detection:
- **Zero library-level data race warnings**
- Test harness races also fixed (`std::atomic` arrays instead of `vector<bool>` / `vector<int>`)
- The production build does **not** include TSan (it adds ~5–10× runtime overhead)

### Stress Testing
- **20 repeated runs** of threading accuracy tests: all passed
- **Mixed-format stress test**: 8 threads × 20 iterations × 5 molecules (SMILES→MOL2→SDF roundtrips) = **800 conversions, 0 errors**

## Using Open Babel from Multiple Threads

### Safe Pattern
```cpp
// Each thread creates its own OBConversion and forcefield instance
std::thread t([&]() {
    OBConversion conv;
    conv.SetInFormat("smi");
    
    OBMol mol;
    conv.ReadString(&mol, "CCO");
    
    OBForceField* proto = OBForceField::FindForceField("MMFF94");
    OBForceField* ff = proto->MakeNewInstance(); // ✅ thread-safe clone
    ff->Setup(mol);
    ff->SteepestDescent(100);
    delete ff;
});
```

### Unsafe Patterns to Avoid
```cpp
// ❌ Reusing the same OBConversion across threads without synchronization
// ❌ Reusing the same forcefield prototype without MakeNewInstance()
```

## Phases Completed

### Phase 1: Foundation ✅
Protected core global state: registries, error log, locale, type tables, chain parser, ring typer.

### Phase 2: Operations & Formats ✅
Fixed `OpEnergy`, `OpMinimize`, `OpConformer`, `Gen3D`, `TinkerFormat`, `MMFF94Charges` to use per-instance forcefields. Fixed `OBConversion::RegisterOptionParam` race. Fixed static string return buffers.

### Phase 3: Forcefields & Copy Optimization ✅
Changed all forcefield `MakeNewInstance` implementations from re-parsing parameter files to copying in-memory parameter vectors. Recovered the ~2.6× performance regression. Fixed `EQEqCharges` lazy-init race.

### Phase 4: TSan Validation & Stress Testing ✅
- Built full project under ThreadSanitizer
- All threading accuracy tests pass with zero library-level warnings
- Mixed-format stress test: 8 threads × 20 iterations × 5 molecules = 800 roundtrips, 0 errors
- 20 repeated runs of threading tests: all passed

## Files Changed

Key files modified for threading:
- `include/openbabel/forcefield.h` — `OBFFParameter` copy semantics
- `include/openbabel/plugin.h` — atomic `AllPluginsLoaded`
- `src/plugin.cpp` — mutex-protected plugin loading
- `src/obconversion.cpp` — mutex on `RegisterOptionParam`, `thread_local` string buffers
- `src/locale.cpp` — thread-local counters
- `src/data.cpp` / `include/openbabel/data.h` — `OBTypeTable` mutex
- `src/forcefields/forcefield*.h` / `*.cpp` — copy-based `MakeNewInstance`
- `src/ops/forcefield.cpp`, `conformer.cpp`, `gen3d.cpp` — per-instance forcefields
- `src/charges/mmff94.cpp`, `eqeq.cpp` — per-instance / atomic initialization
- `src/formats/tinkerformat.cpp` — per-instance forcefield
- `src/transform.cpp`, `src/mol.cpp` — `thread_local` string buffers
- `src/obmolecformat.cpp` / `include/openbabel/obmolecformat.h` — mutex for `--join`/`--separate`/`--compare` statics
- `test/threading_accuracytest.cpp`, `threading_performancetest.cpp` — new tests
- `test/testpdbformat.py` — skip `GetSegName` on older Python bindings

---

## Detailed Task History

This section preserves the original task tracking for reference.

### Exploration

| # | Task | Status | Notes |
|---|------|--------|-------|
| 1.1 | **Map exact registry statics** — List every variable in `plugin.h`, `plugin.cpp`, `obconversion.h`, `obconversion.cpp`, `format.cpp` that is hit during `FindFormat` / `RegisterFormat`. | ✅ Done | `PluginMap()`, `AllPluginsLoaded`, per-type `Map()`, `FormatsMIMEMap()`. Fixed with `std::recursive_mutex`. |
| 1.2 | **Map exact error-log statics** — List every mutable field in `OBMessageHandler` and `obLocale` that can be touched from multiple threads. | ✅ Done | `OBMessageHandler::_messageList`, `_messageCount[]`, `_outputStream`, `_logging`, `_maxEntries`, `_inWrapStreamBuf`, `_filterStreamBuf`. Fixed with `std::mutex`. |
| 1.3 | **Audit `static string` return buffers** — Find every `static string` / `static char[]` used as a function-return buffer in `src/` and `include/openbabel/`. | ✅ Done | Fixed 8 functions with `thread_local static string`. See `src/transform.cpp`, `src/mol.cpp`, `src/obconversion.cpp`, descriptors, fingerprints, ops. |
| 1.4 | **Audit `chains.cpp`** — Determine which static arrays are read-only after init vs mutated during parsing. | ✅ Done | `OBChainsParser chainsparser` has extensive mutable state (`bitmasks`, `visits`, etc.). Protected with `std::recursive_mutex` in `PerceiveChains()`. |
| 1.5 | **Audit `OBMoleculeFormat` deferred state** — Determine exactly when `IMols`, `MolArray`, `_jmol`, etc. are used and whether they can be moved to instance members. | ✅ Done | These are class-statics used only with `-s`, `-j`, and `--compare` options. Now protected by `ClassMutex` (`std::mutex`). |
| 1.6 | **Check existing test coverage** — Verify whether any existing test already exercises multi-threaded paths. | ✅ Done | None existed; we created `threading_accuracytest` and `threading_performancetest`. |

### Tests (Write First)

| # | Task | Status | Notes |
|---|------|--------|-------|
| 2.1 | **Create `test/threading_accuracytest.cpp`** — A new C++ test using `std::thread`. Two sub-tests: (a) parse SMILES → add hydrogens → perceive aromaticity → write canonical SMILES; (b) read MOL2 → `OBForceField::Setup(MMFF94)` → `Energy()` → `SteepestDescent(10)` → `ConjugateGradients(10)` → write MOL2. Compare single-threaded output vs multi-threaded output (must be identical within FP tolerance for energy/coords). | ✅ Done | All 6 parts pass. Verified with TSan (zero library warnings). |
| 2.2 | **Create `test/threading_performancetest.cpp`** — Standalone performance harness. Workload: read MOL2 files → `OBForceField::Setup(MMFF94)` → `SteepestDescent(10)` → `ConjugateGradients(10)` → write MOL2. Measure throughput at 1, 2, 4, 8 threads. Report molecules/sec and scaling efficiency. | ✅ Done | Runnable via `./bin/test_runner threading_performancetest`. ~3.2× speedup at 4 threads on M1 Max. |
| 2.3 | **Create test input data** — Reuse existing `test/files/aromatics.smi` and `test/files/culgi_*.mol2` for now. No new files needed yet. | ✅ Done | Existing files are sufficient. |
| 2.4 | **Fix pre-existing test failures** — `test_cifspacegroup_11` (macOS PDB formatting), `pytest_babel` (missing `obrms` due to Eigen3 detection), `pytest_pdbformat` (Python binding `GetSegName`). | ✅ Done | Fixed FindEigen3.cmake for Eigen 5.x, bumped to C++14, made PDB test assertions platform-tolerant, added `GetSegName` skip for older bindings. |

### Fixes — Phase 1 (Core)

| # | Task | Status | Notes |
|---|------|--------|-------|
| 3.1 | **Protect plugin/format registries** — Add a `std::recursive_mutex` around `PluginMap`, `OBFormat::Map`, `FormatsMIMEMap`, and `AllPluginsLoaded`. Ensure `FindFormat` is safe after initialization. | ✅ Done | Used `std::recursive_mutex` to handle re-entrant calls (e.g., `LoadAllPlugins` → `GetPlugin` → `BaseFindType`). |
| 3.2 | **Protect `obErrorLog`** — Add a mutex to `OBMessageHandler` so that `ThrowError`, `GetMessages`, and `StartLogging`/`StopLogging` are safe. | ✅ Done | Added `std::mutex` to `OBMessageHandler`, locked in all mutable operations. |
| 3.3 | **Fix `obLocale`** — Protect `d->counter` with a mutex or make it atomic. | ✅ Done | Added `std::mutex` to `OBLocalePrivate`. On macOS with `HAVE_USELOCALE`, uses thread-local counters. |
| 3.4 | **Make `ringtyper` thread-local** — Change `OBRingTyper ringtyper` in `src/ring.cpp` to `thread_local`. | ✅ Done | `thread_local OBRingTyper ringtyper` in `src/ring.cpp`. Not exported, so no DSO issues. |
| 3.5 | **Make `ttab` and `resdat` thread-safe** — `OBTypeTable ttab` and `OBResidueData resdat` in `src/data.cpp`. | ✅ Done | Added `std::recursive_mutex` to both classes. Added thread-safe `Translate(to, from, from_type, to_type)` overload. |
| 3.6 | **Eliminate `static string` return buffers in core** — Fix `OBConversion::GetNextFormat`, `OBMol::GetTitle`, `OBBase::ClassDescription`, and others. | ✅ Done | Replaced with `thread_local static string` in 8 functions across core, descriptors, fingerprints, ops, and formats. |
| 3.7 | **Fix `OBMoleculeFormat` statics** — Move `IMols`, `MolArray`, `StoredMolsReady`, `_jmol`, `OptionsRegistered` from unprotected class-statics to mutex-protected access. | ✅ Done | Added `ClassMutex` (`std::mutex`) and `std::lock_guard` around all accesses in `--join`, `--separate`, and `--compare` paths. |

### Fixes — Phase 2 (Operations & Formats)

| # | Task | Status | Notes |
|---|------|--------|-------|
| 3.8 | **Audit `src/ops/forcefield.cpp`** — The `--minimize` operation. Check for static/global state in the minimization op that would race when multiple threads minimize simultaneously. | ✅ Done | `OpEnergy::Do` and `OpMinimize::Do` now use `MakeNewInstance()`. `OpConformer::Do` also uses `MakeNewInstance()`. |
| 3.9 | **`src/formats/tinkerformat.cpp`** — tinker format output with MMFF94 atom types. | ✅ Done | Now uses `MakeNewInstance()`. |
| 3.10 | **`src/ops/gen3d.cpp`** — `MakeNewInstance()` reverted due to segfault from Python subprocess. | ✅ Done | Root cause was `OBFFParameter` copy bug (fixed). Restored `MakeNewInstance()`. Verified 250 subprocess invocations with no crashes. |
| 3.11 | **`src/charges/mmff94.cpp`** — same subprocess issue. | ✅ Done | Restored `MakeNewInstance()`. Verified. |
| 3.12 | **`src/charges/eqeq.cpp`** — lazy param file loading race. | ✅ Done | `_paramFileLoaded` is now `std::atomic<bool>` with double-checked locking. |

### Fixes — Phase 3 (Forcefields & Copy Optimization)

| # | Task | Status | Notes |
|---|------|--------|-------|
| 3.13 | **Fix `OBFFParameter` copy semantics** — `std::length_error: vector` crash during `MakeNewInstance`. | ✅ Done | Assignment operator now properly copies all fields. Implicit copy constructor is safe. |
| 3.14 | **Copy-based `MakeNewInstance`** — Change all forcefields from re-parsing parameter files to copying vectors. | ✅ Done | MMFF94, UFF, Ghemical, GAFF, MM2. Overhead: ~0.19 µs per call (vs ~7.6 ms full pipeline). |

### Integration & Verification

| # | Task | Status | Notes |
|---|------|--------|-------|
| 4.1 | **Ensure `ctest` passes** — Run full test suite after each batch of changes. | ✅ Done | **237/237 tests pass.** |
| 4.2 | **Run ThreadSanitizer (TSan)** — Build with `-fsanitize=thread` and run the new accuracy test to catch data races. | ✅ Done | Full project built under TSan. All 6 threading accuracy tests pass. Zero library-level data race warnings. |
| 4.3 | **Stress testing** — Repeated runs, mixed formats, high thread counts. | ✅ Done | 20 repeated runs: all pass. Mixed-format stress test (8 threads × 20 iterations × 5 molecules = 800 roundtrips): 0 errors in 54s. |
| 4.4 | **Update documentation** — `PLAN.md`, `TASKS.md`, `THREADING.md`, `README.md`. | ✅ Done | All files updated with final state, performance numbers, and usage patterns. |

### Quick-Start Commands

```bash
# Configure
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo

# Build tests
make -j$(nproc) test_runner

# Run existing tests
ctest -j$(nproc)

# Run new threading accuracy test
./bin/test_runner threading_accuracytest 0   # bare threads
./bin/test_runner threading_accuracytest 1   # OBConversion per thread
./bin/test_runner threading_accuracytest 2   # AddHydrogens per thread
./bin/test_runner threading_accuracytest 3   # Canonical write per thread
./bin/test_runner threading_accuracytest 4   # SMILES roundtrip
./bin/test_runner threading_accuracytest 5   # MMFF94 minimization

# Run performance harness
./bin/test_runner threading_performancetest 1
```
