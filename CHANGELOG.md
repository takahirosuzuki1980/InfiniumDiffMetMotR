# Changelog
## [1.3] - 2018-08-07
### Added
- IMAGE motidf database (https://genome.cshlp.org/content/28/2/243.long)

### Fixed
- probeID2position_EPIC.R was integrated to ProbeID2positon.R and deleted.
- Conversion of duplicated probeIDs to position is available. (ProbeID2positon.R)
- dependecy of packages was fixed.
- Date was removed from output file.

## [1.4] - 2018-11-09
### Added
- sampling option
- Output of aignificantly enriched motif plots to a directory

### Fixed
- Numer of core to be used has been changed from 16 to 4.

## [1.5.0] - 2019-10-05
### Update
- Format of scripts was hanged to Roxygen format.

### Added
- When using the sampling option, indicate it on standerd out
- Output DMPs
- lumiMethyNorm support idat files as input
- lumiMethyNorm returns a selDataMatrix object
- MotScr support data.frame as input
- getIdat2M function: get, normalize, and M-value conputation from GEO ID

### Fixed
- Error with the case of DMR = 0/1

## [1.5.1] - 2019-12-27
### Fixed
- debug of probeID2position.R
- debug of seqExtract.R

## [1.5.2] - 2020-05-11
### Fixed
- add "force" option to lumiMethyNorm.R, which is inherited option from readEPIC

## [1.5.3] - 2020-05-23
### Updated
- add "set.seed" to stratSampling.R to fix the randomization
