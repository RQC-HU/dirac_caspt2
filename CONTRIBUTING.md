## Contributing

First, thank you for considering contributing to this project!
You can contribute by opening a pull request or an issue.
Please make sure to read the [Developer's wiki](https://github.com/RQC-HU/dirac_caspt2/wiki/developers-wiki) before you start contributing.
If you have any questions, please feel free to contact us via [GitHub discussions](https://github.com/RQC-HU/dirac_caspt2/discussions)

### How to contribute

- If you find a bug or have a feature request, please [open an issue](https://github.com/RQC-HU/dirac_caspt2/issues/new/choose) or [open a pull request](https://github.com/RQC-HU/dirac_caspt2/compare)

- You can contribute by opening a pull request the following steps below
1. Fork the repository
2. Add your changes
3. Open a pull request
4. Wait for the review

### Setting up the development environment

We set up the development environment using Docker, VSCode devcontainer and GitHub Codespaces.
So even if you don't have prerequisite softwares, you can start developing the project right away.

- If you have Docker, docker compose, VSCode and VSCode remote container extension installed, you can start developing the project by opening the project in VSCode and selecting "Reopen in Container" from the command palette.
- If you don't have Docker, you can use GitHub Codespaces to start developing the project.
  - Please refer to the [GitHub Codespaces documentation](https://docs.github.com/en/codespaces) for more details.
- If you have Docker, docker compose, but don't have VSCode, you can start developing the project by running the following command in the project .devcontainer directory.
  ```sh
  docker compose up -d --build
  docker exec -it dirac-caspt2 /bin/bash
  ```

You can find the details of the development environment in the [Developer's wiki](https://github.com/RQC-HU/dirac_caspt2/wiki/configure-dev-environment)

### Building the project

Please refer to the build instructions at [README.md](README.md#how-to-build)

### How to run tests

- We use pytest to run tests
Please refer to the test instructions at [Developer's wiki](https://github.com/RQC-HU/dirac_caspt2/wiki/about-test)
