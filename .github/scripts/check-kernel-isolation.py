"""Check that the SIMD kernel objects share no code with the library.

GCC and Clang objects must define no shared functions. MSVC objects may
define them with baseline instructions only.

The LLVM tools read the ELF, Mach-O and COFF objects alike, the CI
environment installs them with llvm-tools.
"""

import re
import subprocess
import sys
from pathlib import Path

OBJECT_NAME = re.compile(
    r"CalculationKernel(SSE2|AVX2|AVX512|NEON)(?:\.cpp)?\.(o|obj)$"
)
NM_LINE = re.compile(r"(.+) ([A-Za-z?]) [0-9a-fA-F]*(?: .*)?")
FUNCTION_TYPES = {"T", "W", "i", "I"}
LABEL = re.compile(r"^[0-9a-f]+ <(.+)>:$")
INSTRUCTION = re.compile(r"^\s+[0-9a-f]+:\s+((?:[0-9a-f]{2} )+)\s*([a-z]\S*)")
# Instructions above the baseline that have no VEX prefix.
EXTRA_ISA = {"tzcnt", "lzcnt", "popcnt", "crc32"}


def run(*command):
    result = subprocess.run(
        command, capture_output=True, text=True, errors="replace"
    )
    if result.returncode != 0:
        raise SystemExit(
            f"{' '.join(command)} failed with exit status "
            f"{result.returncode}:\n{result.stderr}"
        )
    return result.stdout


def private_symbol(name, isa):
    return (
        f"xsimd_{isa.lower()}::" in name
        or "(anonymous namespace)" in name
        or "`anonymous namespace'" in name
    )


def kernel_entry(name, isa):
    return f"diffpy::srreal::kernels::{isa.lower()}CalculationKernel(" in name


def defined_functions(path):
    output = run(
        "llvm-nm",
        "--defined-only",
        "--extern-only",
        "--demangle",
        "--format=posix",
        str(path),
    )
    functions = []
    for line in output.splitlines():
        match = NM_LINE.fullmatch(line)
        if not match:
            raise ValueError(f"Unrecognized llvm-nm output: {line}")
        name, kind = match.groups()
        if kind in FUNCTION_TYPES:
            functions.append(name)
    return functions


def disassembly(path):
    """Instructions of the functions in path as (bytes, mnemonic) by
    name."""
    output = run(
        "llvm-objdump",
        "--disassemble",
        "--demangle",
        "--x86-asm-syntax=intel",
        str(path),
    )
    bodies = {}
    current = None
    for line in output.splitlines():
        label = LABEL.match(line)
        if label:
            current = bodies.setdefault(label[1], [])
            continue
        match = INSTRUCTION.match(line)
        if match and current is not None:
            current.append(match.groups())
    return bodies


def non_baseline(body):
    for raw, mnemonic in body:
        # c4 and c5 start VEX instructions, 62 starts EVEX ones.
        if (
            raw.split()[0] in {"c4", "c5", "62"}
            or mnemonic.startswith("v")
            or mnemonic in EXTRA_ISA
        ):
            return mnemonic
    return None


def check(path, isa, msvc):
    functions = defined_functions(path)
    if not any(kernel_entry(name, isa) for name in functions):
        raise ValueError(f"Missing {isa} entry point in {path} (LTO object?)")
    shared = sorted(
        {
            name
            for name in functions
            if not private_symbol(name, isa) and not kernel_entry(name, isa)
        }
    )
    if not msvc:
        return [f"shared code symbol: {name}" for name in shared]
    bodies = disassembly(path)
    errors = []
    for name in shared:
        if not bodies.get(name):
            raise ValueError(
                f"No disassembly for shared function {name} in {path}, "
                f"{len(bodies)} functions disassembled"
            )
        mnemonic = non_baseline(bodies[name])
        if mnemonic:
            errors.append(f"non-baseline instruction {mnemonic}: {name}")
    return errors


def main():
    build_dir = Path(sys.argv[1])
    failures = []
    found = False
    for path in sorted(build_dir.rglob("CalculationKernel*")):
        match = OBJECT_NAME.fullmatch(path.name)
        if not match:
            continue
        found = True
        errors = check(path, match[1], msvc=match[2] == "obj")
        if errors:
            failures.extend(f"{path}: {error}" for error in errors)
        else:
            print(f"PASS {path}")
    if not found:
        failures.append(f"No SIMD kernel objects in {build_dir}")
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
