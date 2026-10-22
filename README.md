## Polymer Pyrolysis Search

## Description
An app to search pyrolysis GC/MS data and identify polymer(s) contained therein.  

## Installation
Download and unzip or clone the repository.  

## Usage
To use the NIST Pyrolysis Polymer Search, simply double click "run-NPPS.bat".  Once in the app select retention index calibration file and click "Create RI calibration" button then load query pyrogram, and click "Run AMDIS" button and when finished click "Polymer Search" button .  The APP will call NIST AMDIS in the background after both the "Create RI calibration" and the "Run AMDIS" buttons.  You can change the search parameters and click "Polymer Search" and it will not call AMDIS again.  A *.MSPEC file (combined file of all individual spectra AMDIS extracted from the data file in the ./NIST_PPS/data/ folder) and *.csv files of results (raw pyrolyzate hitlist, sorted pyrolyzate hitlist, and polymer list in the ./NIST_PPS/results/ folder) are written.  Data files (RI or polymer pyrolysis) cannot have spaces in the name as AMDIS called not run.

## Support
edward.erisman@nist.gov

## Roadmap
-Polymer library increased in size (pyrolysis and reactive pyrolysis)
-possible modification of retention index handling

## Authors and acknowledgment
Edward P Erisman

## License
See [LICENSE.md](LICENSE.md)
