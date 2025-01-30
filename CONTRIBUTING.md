
# Contributing

Thank you for your interest in contributing to the project! Contributions are welcome and greatly appreciated. This guide will help you get started.

## Project Status: Work in Progress 🚧

This project is in an early stage, and we are still finalizing processes for contributions, testing, and documentation. While we appreciate interest in contributing, we ask for patience as we establish guidelines and ensure a stable foundation!

If you have ideas, feedback, or bug reports, we encourage you to **open an issue** rather than directly submitting a pull request for now.

---

## Getting Started

Before making any changes, please **open an issue** to discuss your proposed contribution. This will help ensure your work aligns with the project goals and prevent duplicate efforts.

1. **Fork the Repository**  
   If you are new to the project, start by forking the repository to your GitHub account. You can do this by clicking the "Fork" button at the top right of the repository page.

2. **Clone Your Fork**  
   Clone the forked repository to your local machine:
   ```bash
   git clone https://github.com/<your-username>/fundamentals-of-astrodynamics.git
   cd fundamentals-of-astrodynamics
   ```

3. **Set Up the Upstream Remote**  
   Link the original repository as the upstream to keep your fork up to date:
   ```bash
   git remote add upstream https://github.com/CelesTrak/fundamentals-of-astrodynamics.git
   ```

---

## Making Changes

1. **Create a New Branch**  
   Create a branch for your changes. Use a meaningful name that describes your work:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Your Changes**  
   - Follow the coding standards outlined below and keep consistency with existing code.
   - Add or update documentation as needed (headers, docstrings, etc.).
   - Write tests for any new functionality.


3. **Run Tests and Linting (Python Only)**  
   Before submitting your changes, ensure all tests pass and the code adheres to formatting standards.

   **Run Tests:**
   ```bash
   python -m pytest -v --disable-warnings tests/
   ```
   
   **Check for Linting Issues (flake8):**
   Navigate to the `software/python` directory and run:
   ```bash
   flake8 .
   ```
   
4.   **Check Formatting (black):**
In the `software/python` directory, run:
   ```bash
   black --check --skip-magic-trailing-comma .
   ```
   If any formatting issues are found, run:
   ```bash
   black --skip-magic-trailing-comma .
   ```

4. **Commit Your Changes**  
   Commit your changes with a clear and descriptive message:
   ```bash
   git add .
   git commit -m "Add feature: Description of your feature"
   ```

5. **Push Your Branch**  
   Push your changes to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

---

## Submitting Your Contribution

1. **Open a Pull Request**  
   - Navigate to the original repository and click **"New Pull Request."**
   - Choose your fork and branch as the source and the main repository as the target.
   - Provide a descriptive title and explain your changes in the description.
   - Add appropriate labels*: `bug`, `documentation`, `enhancement`, etc. We also have labels for each programming language, so please tag the relevant language(s) you modified.  
   - Assign yourself to the pull request.
   - Add reviewers: Dave for now, and Samira for Python changes.


2. **Participate in the Review Process**  
   - Respond to feedback and make changes as needed.
   - Ensure all pipeline checks pass before resubmitting.

---

## Coding Standards

- **Code Style**: Follow PEP 8 guidelines for Python code, and maintain consistency with existing code.
- **Testing**: Include unit tests for new features and ensure existing tests pass.
- **Documentation**: Update relevant sections in the documentation.

---

## Reporting Issues

If you encounter a bug or have a feature request:
1. Check the existing [Issues](https://github.com/CelesTrak/fundamentals-of-astrodynamics/issues) to see if it’s already reported.
2. If not, open a new issue with a clear and detailed description.

---

## License

By contributing, you agree that your contributions will be licensed under the [same license](LICENSE) as the project.
