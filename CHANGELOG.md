# Change Log

Notable changes to the project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and the
project adheres to the [Haskell Package Versioning
Policy (PVP)](https://pvp.haskell.org)

## [1.2.0.0] - 2018-04-03
### Changed
  * update for AoS representation of complex numbers
  * improve pure FFT implementation

### Contributors

Special thanks to those who contributed patches as part of this release:

  * Trevor L. McDonell (@tmcdonell)
  * Rinat Striungis (@Haskell-mouse)


## [1.1.0.0] - 2017-09-21
### Changed
  * build against FFTW and cuFFT foreign implementations by default

### Fixed
  * fix to ignore `sh` parameter in inverse mode ([#5])

### Removed
  * Drop support for (deprecated) `accelerate-cuda` backend


## 1.0.0.0 - 2017-03-31

[1.2.0.0]:          https://github.com/AccelerateHS/accelerate-fft/compare/1.1.0.0...1.2.0.0
[1.1.0.0]:          https://github.com/AccelerateHS/accelerate-fft/compare/1.0.0.0...1.1.0.0

[#5]:               https://github.com/AccelerateHS/accelerate-fft/pull/5

