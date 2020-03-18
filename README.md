# Deconvolution

Fourier and wavelet analysis tools in 1 dimension for deconvolution

# Requirements

OCTAVE

# Usage

### Wavelet related functions:
* `[u, v] = filt(type, p)` returns the parent wavelets `u`, `v` used in wavelet transform used of type `type`.
* `w = wtrans(z, type, p)` returns the `p`th resolution wavelet transform of the vector `z` using the parent wavelets of type `type`.
* `z = iwtrans(w, type, p)` returns the `p`th resolution _inverse_ wavelet transform of the vector `w` using the parent wavelets of type `type`. If the `type` and `p` are the same, `iwtrans` and `wtrans` should be inverses of each other upto some _small_ error.
* `coeffs(w, p, q)` returns the `q`-th level wavelet coefficients of the `p`-th resolution wavelet transform. `q` can be from `1` to `p+1`. q=p+1 represents the coarsest wavelet level.
* `B = getbasismat(type, p, N)` generates a matrix of dimension `(p+1)xN` whose j-th row contains the basis element which is translated to generate the vectors used to compute the j-th level wavelet coefficients of a `p`-th resolution wavelet transform of a signal of length `N` (an integer power of 2).

### Parameters:
* `type` can be `shan` for Shannon's wavelets, `d`n for Daubechies wavelets where n can be `2`, `4`, ... , `20`.
* If `type`= `shan`, `p` must divide the length of the signal `z` or `w`.

### Helper functions for wavelet related tools (codes are self-explanatory):
* `filt(type, N)` contains the parent wavelets of various wavelet bases. `N` should be an integer power of 2.
* `v = getother(u)`: given one parent wavelet `u`, returns the other `v` using the conjugate mirror filter lemma.
* `wrec`, `iwrec` are recursive implementation of wavelet transform and inverse wavelet transforms, also called Fast Wavelet Transform. Works best (fastest) on vectors of length of type 2^n for some natural number n.
* `up`, `down` upsamples and downsamples a signal by zeros.
* `fold(z)` folds a vector in half.
* `realconv(a,b)` is the circular convolution of two real  signals (i.e. it returns the real part of the inverse Fourier transform of the product of Fourier transforms of the participating signals)

### Plotting tools:
* `plotcoeffs(w,p)` plots the coefficients of wavelet transform `w`, the `p`-th resolution wavelet transform
* `plotbasis(type, p, N)` (better) plots the basis elements by first computing the full basis matrix using `getbaismat`

<!-- Compression and error-->
<!--* `keeplarge` zeros out smaller values-->
<!--* `w = compress(z, type, p, K)` returns -->
<!--* `relerr(z, typelist, kMax, p, q)` returns the relative error matrix computed in `q`-norm ...-->
<!--* `compareErr(z, typelist, kMax, p, normlist)`-->


### Fourier based Deconvolution
* `[fw, mult] = fdecwien(fsig, fimp, fori, noise, scaling)` returns the Fourier transform `fw` of the Wiener deconvolution with the (usually unknown) signal`mult` returns the Fourier shrinkage parameter used in the deconvolution.
```
 Input: 	fobs: Fourier transform of the noisy observed signal
		fimp: Fourier transform of the impulse response	
		fori: Fourier transform of the guess for the original signal
 			A possible guess will be fobs./sqrt(var(fimp))
		noise: a scalar for noise standard deviation, or the noise data (time domain)
		scaling: A regularized parameter to reduce the ringing, default =1
 Output:	fw:	Fourier transform of the estimate
		mult:	Fourier shrinkage multiplier
```

### Wavelet based deconvolution
* `alpha = getoptsc(z, K, type, p, sigma, rootmethod)` compute the optimal  scaling parameter vector `alpha` of length `p+1` that minimizes the error in ForWaRD algorithm. Here `rootmethod` can be either `search` (for searching through uniformly located numbers between 0 and 1) or `bisec` (for a bisection method, much faster and accurate)

### ANITA related tools
* `getantenna(n)` outputs the filename of the impulse response of the antenna, given the absolute index of the antenna
* `[ax, aximp] = prepsig(num)` outputs matrix `ax` with noisy blurred ANITA signals and `aximp` contains the impulse responses used in the observation
* `[snrval, M, m, a, b] = getsnr(z)` computes the signal-to-noise ratio of a temporally localized signal `z`. Here, `M` and `m` are the max and the min of the signal, the time interval `[a,b]` has 10 times the length of the peak region and is a neighborhood of the peak region.
<!--* `wpxsc = getwpxsc(type, p, sigma)` returns the optimal scaling values `alpha` required for the deconvolving the ANITA signal `wpx` (see `getdata` for the test signal), using the bisection method -->

* `plotfanita(z)` plots the absolute value of the Fourier transform of an ANITA signal in the frequency domain in MHz unit
* `plotfanitaS(z)` (with a different sample rate) plots the absolute value of the Fourier transform of an ANITA signal in the frequency domain in MHz unit

### Other tools
* `z = padfreq(w, N, eps)` increases the resolution of the signal (i.e. upsamples) `w` by zero-padding followed by smoothing. If `N > length(w)` it adds zeros in the frequency domain. However, the process can add Gibbs phenomenon in the time domain if the frequency does not decay to zero. To avoid this, we apply a smooth Fourier cut-off function (a Plank-taper window in the frequency domain) with `eps` (between 0 and 1) being the proportion of the _actual_ frequency of the signal `w` to get scaled down. For example, `eps = 1` scales the entire frequency range of the original signal. 
To apply the Plank-taper window without padding, use `N = length(w)` and any value of `eps`.
* `planktaper(N,eps)` generates a Plank-taper window of length `N` and tapering window of length `eps*N`. E.g. `eps = 0.5` the function obtains value 1 on a set of measure zero
* `croscor(f,q)` computes the circular cross-correlation between two vectors of same length
<!--* `cpuow(z)` computes the cumulative power distribution of a signal `z`-->
* `deriv(z)` computes the difference of consecutive terms (a crude derivative)

### Demo
* `set_params` specifies the parameters to be used. This also produces the basis matrix for the wavelet transform.
* `set_data` prepares a set of test data. The original simulates signal is called `testyori`. It produces the following matrices, rows (15 of those) of which are individual signals.
 1. `ax` contains the ANITA observations (blurred and noisy)
 1. `aximp` contains the corresponding impulse responses
 1. `noiseax` contains the noise data
 1. `wax` contains the simulated noisy blurred signals (`conv(testyori, aximp(i)) + noiseax(i)`)


# Reference
* M. Frazier, _An Introduction to Wavelets through Linear Algebra_
* Thanks to Dr. Peter Gorham for providing support and ANITA data and Manuel Olmedo for making the data matlab-friendly
