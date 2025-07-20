# FITSelector
A python utility for selecting FITS files from a collection, based on FITS Header keywords and values, 
along with some statistics and map plotting capabilities.  This is especially useful for selecting specfic
FITS files from a large collection prior to stacking, such as when working with Smart Telescope data from 
diverse sources.  The original motivation for writing this program came from the my work with 
The Seestar Collective Discord community. 

## Installation

1. Click on the Green "CODE" button and either download a zip file, or copy the HTTPS address and run:

   ```bash
   cd /directory/where/you/want/the/code
   git clone https://github.com/Astrobirder/FITSelector
   ```
   
2. Create a virtual environment and activate it:

   ```bash
   python3 -m venv .venv
   .venv/scripts/activate
   ```

3. Install the required dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Launch FITSelector:

   ```bash
   python FITSelector.py   # Shows help message

   Try these 2 versions on Seestar Files:
   python3 FITSelector.py -i /path/to/your/fits/files -t "TELESCOP" -p
   python3 FITSelector -i /path/to/your/fits/files -t "ALLHEADERS"
   ```

## Usage

   ```bash
   usage: FITSelector.py [-h] [-i INPUTDIR] [-o OUTPUTDIR] [-k KEY=VALUE [KEY=VALUE ...]]
             [-t, --tablekey TABLEKEY] [-p, --plot] [-m, --move | -c, --copy] [-l, --log]

This program lets you select FITS files based on specific FITS Header keywords and their values.
You can then generate statistics, plot data on a map and move or copy the files that match
the optional keywords.
HINT: On Seestar FITS files, TRY:
python3 FITSelector.py -i /path/to/your/fits/files -t "TELESCOP" -p

Also try the special tableky "ALLHEADERS" to get a list of all headers in all files:
python3 FITSelector -i /path/to/your/fits/files -t "ALLHEADERS"

options:
  -h, --help            show this help message and exit
  -i, --inputdir INPUTDIR
                        The directory containing your FITS files
  -o, --outputdir OUTPUTDIR
                        The directory where matching FITS files will be moved/copied.
                        If unspecified, inputdir/output will be used
  -k, --keywords KEY=VALUE [KEY=VALUE ...]
                        FITS Header Keyword and value in the format -k "EQMODE=1".
                        You can use multiple -k "FITSKEY1=VALUE1 FITSKEY2=VALUE2" keys
                        if you want all to be true to select a file. You can also use
                        -k "FITSKEY1>=VALUE1 FITSKEY2<=VALUE2" to use comparison operators.
  -t, --tablekey TABLEKEY
                        Create a summary table with a count of all the files for each unique
                        value of the "tablekey". NOTE: if used with -k, all the files must
                        also match all keywords. If you use ALLHEADERS as the tablekey,
                        a list of all the headers in all of the files that are in inputdir
                        that match all (optional) keywords, with a count of the files that
                        contain each keyword.
  -p, --plot            Plot the count of matching files on a map. 
  -m, --move            Move files that match all keywords into the outputdir.
  -c, --copy            Copy files that match all keywords into the outputdir.
  -l, --log             Log debug info to inputdir/FITSelector.log
  ```
