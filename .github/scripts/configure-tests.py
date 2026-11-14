"""Write the conda environment file for the CI tests."""

import json
from pathlib import Path

# llvm-tools for check-kernel-isolation.py
requirements = ["python", "llvm-tools"]
for filename in ("build.txt", "conda.txt", "tests.txt"):
    for line in (Path("requirements") / filename).read_text().splitlines():
        spec = line.split("#", 1)[0].strip()
        if spec:
            requirements.append(spec)

Path("build").mkdir(exist_ok=True)
Path("build/test-environment.yml").write_text(
    json.dumps(
        {
            "name": "test",
            "channels": ["conda-forge"],
            "dependencies": requirements,
        },
        indent=2,
    )
    + "\n"
)
