# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).



##  version [1.0.0](https://github.com/i2bc/orfmine/releases) - 2022-11-25
### Added
- softwares.ini file to inform on absolute local path to tango and/or iupred2a


### Changed
- use of docker image through a handler process - commands to use docker image are the same as without docker (breaking)
- output directory options for all orfmine callable scripts
- rename output of orfplot (replace "output.png" by "outdir/filename.png")
- place -o argument option in orfget with -outdir & -outname in orfget (not required by default)


##  version [0.8.7](https://github.com/i2bc/orfmine/releases) - 2022-11-12
### Added
- docker image of ORFmine (breaking)

### Changed
- whole new structure dependent on docker container filesystem  (breaking)


### Removed
