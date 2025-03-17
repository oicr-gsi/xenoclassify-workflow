# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.6.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to medata only).

## [1.5.1] - 2024-04-08
### Changed
- Changed default filterSupplemental to tru for sorting task.

### Fixed
- Changed how the above task handled flag (it appears that previous version did not have the right implementation for filterSupplemental flag check).

## [1.5.0] - 2024-03-08
### Changed
- Updated imported bwaMem (to bwamem2 2.2.1) and star (to 2.3.0).
- Updated samtools from 1.9 to 1.14.

### Added
- Flag for removing supplemental alignments (ON by default to prevent failures when using updated STAR).

## [1.4.0] - 2022-06-15
### Added
- Assembly-specific modules are specified inside the workflow.

## [1.3.0] - 2022-03-08
### Added
- Multi-lane STAR support for WT/MR data, some re-design of the workflow.

## [1.2.1] - 2022-02-08
### Added
- Added support for TS, EX in the workflow.

## [1.2] - 2022-01-18
### Added
- Adding ability to run on WT data (using imported STAR).

## [1.1.1] - 2021-05-31
### Changed
- Migrating to vidarr.

## [1.1] - 2020-10-21
### Added
- Added an ability to remove tags on a list, bwaMem is now a sub-workflow, .json report with totals for different tags is also provisioned.

## [1.0] - 2020-05-07
### Added
- Initial Release.
