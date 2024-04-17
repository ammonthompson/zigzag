# Changelog

All notable changes to the zigzag project are documented in this file.

## [0.1.0] - 2024-03-10

### Added
- Initial release of the package. Used in publication: [Thompson et al., 2020, PNAS](https://doi.org/10.1073/pnas.1919748117)



## [1.0.0] - 2024-04-16

### Added
- Model with all mixture components share variance param (whether active or inactive).
- Enhanced posterior predictive simulation plotting features
-- Posterior predictive plots show inactive and active component distributions


### Changed
- By default mixture component now share a single variance parameter instead of one for the inactive component and one shared by all active components.
-- Need to set shared_variance = FALSE to set model to original default model in v0.1.0.
- Reorganized tutorial and simulation directories.


### Fixed
- minor bugs associated with pdf devices not closing properly and number of generations run being two cycles too high.
- fixed bug in append for mcmc and burnin
- Fixed problem from specifying prior thresholds on component means without specifying number of components. Zigzag now infers the number of active components from the number of prior thresholds in threshold_a. Previously, zigzag set threshold_a to 'auto' and disregarded the manual setting if num_active_components was 'auto'.



