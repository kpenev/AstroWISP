A list of possible improvements {#possible_improvements_page}
===============================

- Introduce non-parametric PSF
    - Piecewise constant
    - Bi-cubic, requiring \f$C^1\f$ continuity
- Correct external source positions for sub-pixel map effects
    - After fitting for PSF one can analytically calculate the offset
      introduced by the sub-pixel map for each source and re-fit the PSF.
      Iterater to convergence.
    - Re-derive the astrometric solution simultaneously with the PSF,
      assuming the source match is correct.
- Implement non-piecewise constant sub-pixel maps
    - The whole map can be some finite order polynomial
    - Each sub-pixel can be a polynomial
