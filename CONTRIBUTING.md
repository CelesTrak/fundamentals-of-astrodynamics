# Contributing

Thank you for your interest in contributing to the Fundamentals of Astrodynamics code repository! This guide explains how to get involved effectively.

---

## Before You Start

1. **Check the [existing issues](https://github.com/CelesTrak/fundamentals-of-astrodynamics/issues)** to see if the bug or feature has already been reported or is actively being worked on.
2. **Read the relevant README** for the language you are contributing to (e.g., [`software/python/README.md`](software/python/README.md)).
3. **Open an issue** (see below) and wait for maintainer acknowledgment before writing code. This prevents wasted effort if the change is out of scope or already in progress.

---

## Opening an Issue

Use the GitHub issue templates to report a bug or propose a feature. Choose the appropriate template when creating a new issue:

- **Bug report** – Something is broken or producing incorrect results.
- **Feature request** – You want a new function, option, or enhancement.

After opening an issue, a maintainer will triage it and confirm whether a pull request would be welcome. Once you have a green light, you can proceed with the steps below.

---

## Setting Up Your Environment

1. **Fork the repository** by clicking the "Fork" button on the GitHub page.

2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/<your-username>/fundamentals-of-astrodynamics.git
   cd fundamentals-of-astrodynamics
   ```

3. **Add the upstream remote** so you can keep your fork in sync:
   ```bash
   git remote add upstream https://github.com/CelesTrak/fundamentals-of-astrodynamics.git
   ```

4. **Install Python dependencies** (if contributing to the Python package):
   ```bash
   cd software/python
   pip install -e .
   pip install black==24.10.0 flake8==7.0.0 pytest==8.1.1 pytest-cov==5.0.0
   ```

---

## Making Changes

1. **Sync your fork** with the latest upstream changes before starting:
   ```bash
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

2. **Create a branch** using one of the prefixes below and a short, descriptive name:

   | Prefix | Use for |
   |--------|---------|
   | `feature/` | New functions or capabilities |
   | `fix/` | Bug fixes |
   | `docs/` | Documentation-only changes |
   | `refactor/` | Code improvements with no behavior change |
   | `test/` | Adding or updating tests |

   ```bash
   git checkout -b feature/add-mu-param
   git checkout -b fix/rv2coe-edge-case
   ```

3. **Follow the coding standards** outlined in the [Coding Standards](#coding-standards) section below.

4. **Add or update tests** for any new or changed functionality.

5. **Update docstrings and documentation** as needed.

---

## Testing and Linting (Python)

Before submitting, make sure all checks pass. From the `software/python` directory:

**Run tests:**
```bash
python -m pytest -v --disable-warnings tests/
```

**Check formatting:**
```bash
black --check --skip-magic-trailing-comma .
```
If formatting issues are found, auto-fix with:
```bash
black --skip-magic-trailing-comma .
```

**Check linting:**
```bash
flake8 --statistics .
```

These same checks run automatically in CI on every push and pull request.

---

## Submitting a Pull Request

1. **Commit your changes** with a clear message:
   ```bash
   git add <files>
   git commit -m "Add gravitational parameter support to rv2coe and coe2rv"
   ```

2. **Push your branch** to your fork:
   ```bash
   git push origin feature/add-mu-param
   ```

3. **Open a pull request** against `main` on the original repository. Use the PR template and include:
   - A summary of what changed and why.
   - A reference to the issue: `Closes #123`
   - The language(s) affected.

4. **Add labels** for the type of change (`bug`, `enhancement`, `documentation`, etc.) and the language (`python`, `matlab`, `c#`, etc.).

5. **Assign yourself** to the PR.

Reviewers are automatically assigned based on the files you changed (see [CODEOWNERS](.github/CODEOWNERS)). All CI checks must pass before a PR can be merged.

---

## Review Process

A maintainer will review your PR and may request changes. Please:

- Respond to feedback promptly and push updates to the same branch.
- Keep the PR focused — one logical change per PR makes review easier.
- Do not force-push after review has started, as it makes it harder to see what changed.

---

## Coding Standards

- **Style**: Follow [PEP 8](https://peps.python.org/pep-0008/) for Python. Use `black` and `flake8` to enforce this automatically.
- **Docstrings**: Use [Google-style docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) with `Args`, `Returns`, and `References` sections. Reference the relevant section of Vallado (2022) where applicable.
- **Testing**: Write `pytest` unit tests for any new function or behavior change. Validate against examples from Vallado (2022) where possible.
- **Consistency**: Match the conventions of the surrounding code in whichever language you are contributing to.

---

## Reporting Issues

If you encounter a bug or have a feature request but are not planning to fix it yourself:
1. Check the [existing issues](https://github.com/CelesTrak/fundamentals-of-astrodynamics/issues) to avoid duplicates.
2. Open a new issue using the appropriate template.

---

## License

By contributing, you agree that your contributions will be licensed under the [GNU Affero General Public License v3.0](LICENSE), the same license as this project.