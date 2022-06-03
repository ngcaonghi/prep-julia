# prep-julia

A Julia implementation of the PREP EEG preprocessing pipeline designed by Shamlo *et al.* (2015) [[1]](#1).

## Features and checklist
- [x] Choose to preprocess either a single file or a folder of .csv files
- [x] Denoise data by bandpass and linenoise filter 
- [x] Detect bad channels by NaN, flatness, deviation, high-frequency noise, and correlation
- [ ] Rereference channels by robust referencing (in progress)
- [ ] Interpolate bad channels by spherical spline interpolation (in progress)
- [ ] Detect bad channels by random sample consensus (start and finish dates TBD)
- [ ] Return preprocessed data in .csv format (start and finish dates TBD)
- [ ] Plot before-and-after raw data and power density spectra (start and finish dates TBD)

## Getting started
- Run `TerminalInterface.jl` by typing `julia` in the terminal and typing the file name. In VSCode, the file can be run by Ctrl+F5.
- Choose "1" to detect the bad channels of only 1 file in the sample folder, or "2" for the whole folder.
- Type the local directory of the file or the folder, *e.g.*, `/Users/username/Documents/prep-julia/samples`.
- Wait for the result! Right now only bad channel detection is fully implemented. In the future when the whole pipeline has been completed, the resulted data will be saved in the .csv format.

## Data source
The data provided in the samples folder are extracted from the BED dataset by Arnau-González *et al.* (2021) [[2]](#2). Find the full dataset [here](https://zenodo.org/record/4309472?token=eyJhbGciOiJIUzUxMiIsImV4cCI6MTY1ODUyNzE5OSwiaWF0IjoxNjMwMzk5NjcyfQ.eyJkYXRhIjp7InJlY2lkIjo0MzA5NDcyfSwiaWQiOjE2ODExLCJybmQiOiIyZmU5Nzk0ZiJ9.mCdQaX9123h0Cm37l2qPq9FFrC_g0D5YRW1R5ztilrRd_TI9ssvpw-hUl17sN4wU8DI6E7C0LqzZ-diYOaZDGg#.YXimL3Uzb1E).

## Reference
<a id="1">[1]</a>
Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for large-scale EEG analysis. *Frontiers in neuroinformatics*, 9, 16.

<a id="2">[2]</a>
Arnau-González, P., Katsigiannis, S., Arevalillo-Herráez, M., & Ramzan, N. (2021). BED: A New Data Set for EEG-Based Biometrics. *IEEE Internet of Things Journal*, 8(15), 12219-12230.
