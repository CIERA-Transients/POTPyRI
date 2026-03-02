# Contributing to POTPyRI

Thank you for your interest in contributing to POTPyRI. We welcome issues, suggestions, and pull requests.

## Reporting issues

If you encounter a bug, have a question about usage, or have a feature request:

- **GitHub Issues:** Open an issue at [github.com/CIERA-Transients/POTPyRI/issues](https://github.com/CIERA-Transients/POTPyRI/issues). Please include a clear description, steps to reproduce (for bugs), and your environment (Python version, instrument, OS) where relevant.
- **Email:** You can also contact the developers at `ckilpatrick@northwestern.edu` for matters you prefer not to discuss in the issue tracker.

## Contributing code or documentation

1. **Fork the repository** and create a branch from the default branch.
2. **Make your changes** and add or update tests if applicable. Run the test suite with `pytest tests` (use `pytest tests -m "not integration"` for offline runs).
3. **Open a pull request** against the main repository. Describe your changes clearly and reference any related issues.
4. **Code style:** Follow the existing style in the codebase. The project uses standard Python packaging and type hints where appropriate.

Instrument-specific changes (e.g. new instruments or header/sorting logic) should include a brief justification and, if possible, a note in the PR description on how the change was tested.

## Support

- For **usage questions, bugs, or feature requests,** please open a [GitHub issue](https://github.com/CIERA-Transients/POTPyRI/issues).
- For **other inquiries** (e.g. collaboration, adding an instrument), you can contact the developers at `ckilpatrick@northwestern.edu`.

We aim to respond to issues and pull requests in a timely manner; please allow at least a week for non-urgent responses.
