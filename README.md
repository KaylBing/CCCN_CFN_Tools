This is the supporting package for the projects in the papers listed below.

Please feel free to reach out with errors, critiques, or questions.

## Installation

Before proceeding, ensure you have R and RStudio installed. You will also need the devtools package, which can be installed with:

install.packages("devtools")

Once devtools is installed, you can install this package directly from GitHub using:

devtools::install_github("your_github_username/repository_name")

## Usage

After installation, load the package in R with:

library(packageName)

You can then use the available functions as described in the package documentation.

## Development & Contribution

If you wish to modify or contribute to the package, clone the repository locally using:

git clone https://github.com/your_github_username/repository_name.git

Then, in R, navigate to the package directory and use devtools to load the package for development:

devtools::load_all()

To check for issues before committing changes, run:

devtools::check()

