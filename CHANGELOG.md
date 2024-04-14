# Changelog

All notable changes to the zigzag project are documented in this file.

## [0.1.0] - 2024-03-10

### Added
- Initial release of the package. Used in publication: [Thompson et al., 2020, PNAS](https://doi.org/10.1073/pnas.1919748117)

## [1.0.0] - 2024-04-01

### Added
- Model with all mixture components share variance param (whether active or inactive).
- Enhanced posterior predictive simulation plotting features
-- Posterior predictive plots show inactive and active component distributions


### Changed
- By default mixture component share a single variance parameter
-- Need to set shared_variance = FALSE for original default model in v. 0.1.0


### Fixed
- minor bugs associated with pdf devices not closing properly and number of generations run being two cycles too high.



