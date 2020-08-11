<div align="center">
<img width="450" src="https://github.com/AccelerateHS/accelerate/raw/master/images/accelerate-logo-text-v.png?raw=true" alt="henlo, my name is Theia"/>

# FFT components for the Accelerate language

[![GitHub CI](https://github.com/tmcdonell/accelerate-fft/workflows/CI/badge.svg)](https://github.com/tmcdonell/accelerate-fft/actions)
[![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/AccelerateHS/Lobby)
<br>
[![Stackage LTS](https://stackage.org/package/accelerate-fft/badge/lts)](https://stackage.org/lts/package/accelerate-fft)
[![Stackage Nightly](https://stackage.org/package/accelerate-fft/badge/nightly)](https://stackage.org/nightly/package/accelerate-fft)
[![Hackage](https://img.shields.io/hackage/v/accelerate-fft.svg)](https://hackage.haskell.org/package/accelerate-fft)

</div>

FFT library for the embedded array language Accelerate. For details on
Accelerate, refer to the [main repository][GitHub].

The following build flags control whether optimised implementations are used.
Note that enabling these (which is the default) will require the corresponding
Accelerate backend as a dependency:

  - `llvm-ptx`: For NVIDIA GPUs
  - `llvm-cpu`: For multicore CPUs

Contributions and bug reports are welcome!<br>
Please feel free to contact me through GitHub or [gitter.im][gitter.im].

  [GitHub]:       https://github.com/AccelerateHS/accelerate
  [gitter.im]:    https://gitter.im/AccelerateHS/Lobby

