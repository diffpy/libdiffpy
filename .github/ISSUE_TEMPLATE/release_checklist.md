---
name: Release
about: Checklist and communication channel for GitHub and conda-forge release
title: "Ready for <version-number> GitHub/conda-forge release"
labels: "release"
assignees: ""
---

### GitHub rc-release preparation checklist:

- [ ] All PRs/issues attached to the release are merged.
- [ ] All the badges on the README are passing.
- [ ] License information is verified as correct. If you are unsure, please comment below.
- [ ] Locally rendered documentation contains all appropriate pages and other human-written text is up-to-date with any changes in the code.
- [ ] Installation instructions in the README, documentation, and the website are updated.
- [ ] Grammar and writing quality are checked (no typos).
- [ ] `FALLBACK_VERSION` in `site_scons/fallback_version.py` is updated to `<version-number>.post0`. Without this, the SCons build fails after tagging.
- [ ] `scons test` passes locally with the latest conda-forge compilers and dependencies.
- [ ] Dispatch matrix testing to test the release on all systems. If you do not have permission to run this workflow, tag the maintainer and say `@maintainer, please dispatch matrix testing workflow`. Matrix testing also runs automatically after a GitHub release is published, but it does not block the release, so check that it passes.

Please tag the maintainer (e.g., @username) in the comment here when you are ready for the GitHub release. Include any additional comments necessary, such as version information and details about the pre-release here:

### GitHub full-release preparation checklist:

- [ ] The rc GitHub release and its source tarball appear correctly on the releases page.
- [ ] Docs are deployed successfully to `https://www.diffpy.org/libdiffpy`.
- [ ] Successfully run all tests and do functional testing, e.g., with `diffpy.srreal` built against the rc.

Please let the maintainer know that all checks are done and the package is ready for full release.

### conda-forge release preparation checklist:

<!-- After the maintainer releases on GitHub, please check the following when creating a PR for conda-forge release.-->

- [ ] Ensure that the full release has appeared on GitHub successfully.
- [ ] `version` and `sha256` of the GitHub source tarball are updated in `meta.yaml` in the feedstock.
- [ ] New package dependencies listed in `build.txt`, `conda.txt` and `tests.txt` are added to `meta.yaml` in the feedstock.
- [ ] `license_file` in `meta.yaml` points to `LICENSE.rst` (renamed from `LICENSE.txt`).
- [ ] Close any open issues on the feedstock. Reach out to the maintainer if you have questions.
- [ ] Tag the maintainer for conda-forge release.

### Post-release checklist

<!-- Before closing this issue, please complete the following: -->

- [ ] Install the released package from conda-forge and conduct functional testing using the installation guide in the README. Attach screenshots/results as comments.
- [ ] Documentation (README, API references, and websites) is deployed without broken links or missing figures.
