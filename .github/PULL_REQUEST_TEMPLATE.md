## Description

<!-- What does this PR do and why? -->

## Type of change

- [ ] Bug fix
- [ ] New feature / module
- [ ] Configuration change
- [ ] Refactor
- [ ] Documentation only

## Checklist

**All PRs**
- [ ] Pre-commit checks pass locally (`pre-commit run --all-files`)
- [ ] Pipeline stub tests pass (`nf-test test tests/default.nf.test --filter "stub"`)
- [ ] No large files (>5 MB) accidentally staged

**New or modified modules**
- [ ] Module tests added/updated (`nf-test test modules/local/<name>/tests/`)
- [ ] Stub block present and correct
- [ ] `versions.yml` emitted
- [ ] `publishDir` configured in `conf/modules.config`

**Workflow / config changes**
- [ ] Full pipeline test passes (`nf-test test tests/default.nf.test --filter "three B. uniformis genomes$"`)
- [ ] Stub task counts updated if processes were added or removed
- [ ] `README.md` updated if parameters or outputs changed
- [ ] `CLAUDE.md` updated if architecture or conventions changed

**Notebook changes**
- [ ] Notebook tests pass (`task test:notebooks`)
- [ ] `README.md` updated if render parameters changed
