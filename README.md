# FPW: A Fast Periodogram for Unevenly Spaced Data

FPW is a Python package for computing periodograms from unevenly spaced time series data. It provides a fast and efficient method for identifying periodic signals in astronomical and other time-domain datasets by fitting a series of piecewise-continuous boxes.

## Installation

You can install FPW via pip:

```sh
pip install fpwperiodic
```

## Usage

Import FPW and run the periodogram analysis using:

```python
import fpw

# Example usage
periodogram = fpw.run_fpw(times, flux, fluxerrs, freq_grid, N_bins)
```

### Parameters
- `times`: Array of time observations.
- `flux`: Array of observed flux values.
- `fluxerrs`: Array of uncertainties in flux.
- `freq_grid`: Array of frequencies to evaluate the periodogram.
- `N_bins`: Number of bins for phase folding. More bins is better for eclipses and complex waveforms, less for sinusoids.

### Returns
The function returns the computed periodogram values corresponding to the provided `freq_grid`.

## Contributing
If you'd like to contribute to FPW, feel free to submit a pull request or open an issue on GitHub.

## License
This project is licensed under the MIT License.

## Related Paper
For more details, please refer to our published paper: [Finkbeiner et al. 2025](https://doi.org/10.48550/arXiv.2502.00243).

## Contact
For questions or issues, please open a GitHub issue or reach out to the maintainers.

