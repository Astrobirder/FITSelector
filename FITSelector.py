#
#  This python script sorts through FITS files (Flexible Image Transport System) that are used in 
#  astronomical imaging and related fields.  It selects files based on FITS Header Keywords and their values.
#  You may then either print a table with counts of the matching files, and/or move/copy the selected files to a destination directory.
#
#  The script takes the following arguments:
#
#  -i, --inputdir : The directory containing the fits files to be scanned/summarized/moved or copied
#  -k

import numpy as np
import pandas as pd
import argparse
import os
from astropy.io import fits
import signal
import sys
import glob
import pprint
import re
import operator
from collections import defaultdict
from argparse import ArgumentParser


inputdir='.'                # the directory where your FITS files are stored
outputdir='./output'        # destination directory where FITS files will be optionally moved or copied
keyword_ndarray=[]          # a 2D numpy ndarray for specified matching Keyword/value combinations \
                            # each row contains a Keyword, Comparison operator, and Value
                            # note, Keywords and Comparison operators are stored as strings, Values are either strings or float.
tablekey=''                 # FITS Header variable used to create a summary of all values of said key across a set of FITS files.
hdu=0                       # Header Data Unit setting to pull the FITS Header from each file.
args=''                     # ArgumentParser object containing the command line arguments
total_fits_files=0          # Counter for total fits files found in inputdir
movefiles=False             # Flag for moving files to outputdir
copyfiles=False             # Flag for copying files to outputdir
operator_map = {            # Map for converting keyword comparison values into operators for math operations.
    "==": operator.eq,
    "!=": operator.ne,
    "<": operator.lt,
    ">": operator.gt,
    "<=": operator.le,
    ">=": operator.ge
}
keyword_match_count = 0     # Counter to keep track of potentially multiple keyword matches.
all_keywords_match = False  # Flag for when all keywords are matched.


def main():

    # Make sure the argument variables used here are all using the globals
    global inputdir
    global keyword_ndarray
    global outputdir
    global tablekey
    global args
    global movefiles
    global copyfiles
    global keyword_match_count
    global all_keywords_match
    global process_counter
    global matched_files_counter

    signal.signal(signal.SIGINT, signal_handler) #Handle CTRL+C
    args=get_arguments()
    # print('Args returned are: \n')
    # pprint.pprint(args)
    # print('\n')
    #print('Args are inputdir: ',args.inputdir, ' Keywords: ', args.keywords, ' outputdir: ',  args.outputdir, ' table: ', args.table)
    #print('ON return, inputdir: ',{inputdir}, ' Keywords: ', keyword_ndarray, ' outputdir: ', {outputdir})

    # Verify that arguments are acceptable
    if not os.path.isdir(inputdir):
        print(f'Directory ',{inputdir},' does not exist.')
        exit(1)
    
    os.makedirs(outputdir, exist_ok=True) # if outputdir doesn't exist, make it.

    # setup have_keywords for easy checking below.
    # print('Args.keywords is: ')
    # pprint.pprint(args.keywords)

    if (args.keywords) is None:
        #print ('No Keywords!')
        have_keywords = False
    else:
        #print('Have Keywords!')
        have_keywords = True


    # Count how many fits files we have:
    global total_fits_files
    total_fits_files=count_files_with_glob(inputdir, '*.fit*')

    print('Processing FITS files in ',inputdir)
    print('Total FITS files: ', total_fits_files, '\n')

    if (total_fits_files<1):
        print('Directory doesn\'t contain any valid files')
        print('Please use a directory with either *.fit or *.fits files')
        exit(0)

    if tablekey:
    # Main option Tablekey processing--generating summary table.
        #print('We have a tablekey! and it is: ', tablekey)

        # use a defaultdict for the results table, with tablekey values and counts
        results=defaultdict(int)
        # We have a tablekey, so generate a summary table based on differing tablekeys across all the files.
        # Iterate across all the files, then iterate across all the keywords, testing if the headers match all keywords
        #print("Starting Loop!")
        # Loop through all FITS files in the source directory


        #print ('\nWe have ', keyword_rows, ' rows and ', keyword_cols, ' columns in our keywords.')
        #pprint.pprint(keyword_ndarray)
        #print('\n')

        
        process_counter = 0
        matched_files_counter = 0
        for filename in os.listdir(inputdir):
        # loop through files looking for keyword matches
            # print(' \n ' )
            #print('Filename is: ',{filename})
            #user_input = input("Hit <enter> to process")

            # If group below filters for only .fit or .fits files.
            if filename.lower().endswith('.fits') or filename.lower().endswith('.fit'):
            # Make sure we only use filenames/filepaths that are fits files!
                filepath = os.path.join(inputdir, filename)
                #print({filepath})
                try:
                    header = fits.getheader(filepath, hdu)
                    #print('Got Header from: ', filepath)
                    #pprint.pprint(header)
                    #print('\n\n')
                    # for row in keyword_ndarray:
                    #     if keyword_ndarray[row,0] in header:
                    #         keyword_match_count += 1
                    #         print(filename, ' has keyword ', keyword_ndarray[row,0])
                    #     else:
                    #         print(filename, ' MISSING keyword ', keyword_ndarray[row,0])
                     
                    # Iterate across the keywords now:
                except Exception as e:
                    print(f"Error reading {filename}: {e}, skipping it")
                    continue


            # If group checks if we have keywords  
            # and for loop goes through all keywords and sets all_keywords_match if they do
            if (have_keywords):    
            #Verify that all keywords in the file under test match the requirements
                keyword_match_count = 0 # zero our counter
                all_keywords_match = False
                # Counting required keywords:
                keyword_rows, keyword_cols = keyword_ndarray.shape

                # print('Here''s our keyword_ndarray: ')
                # pprint.pprint(keyword_ndarray)

                # for loops through all keywords in the keyword_ndarray looking for matches
                for row in keyword_ndarray: # remember that row here is a 1-D slice 
                # Search through all the keywords and see if they all match the file in question.
                    #print('Checking keyword ',row)
                    # pprint.pprint(row)
                    if (row[0]) in header:  # if the current keyword (column 0 is keyword text) exists in the header of the file under test
                        # Convert the operator in the row[1] from a string to something usable
                        #print ('Filename: ', filename, ' has ', row[0], ' in header')
                        op_function = operator_map.get(row[1])
                        # Compare keyword value (row[2]) to the value in the file to see if it passes
                        if op_function(row[2],header[row[0]]):
                            keyword_match_count += 1

                if (keyword_match_count==keyword_rows):
                    all_keywords_match = True
                    #print('We did!')
                    # All the keywords have matched!
                    # Check to see if the tablekey exists in the header of the file under test
        
            # If every keyword matches or we aren't using keywords, add the file's tablekey to results
            if(all_keywords_match or not have_keywords):
            # everythings good with this file, add it to the table
                matched_files_counter += 1
                if (tablekey) not in header:
                    #print (tablekey, ' is not in file ', filename, 'Setting it to None' )
                    tablekey_value='None'
                if (tablekey) in header:  # if the key exists, put the value out of the header
                    tablekey_value=header[tablekey]
           
                # print('\n results before modification: ')
                # pprint.pprint(results)
                # Now we have the "row" we need to increment or add to the results
                # See if the entry already exists, and if so, increment the count, otherwise, create it and set count to 1
                if (tablekey_value) in results:
                    results[tablekey_value] += 1
                    # print('Incrementing ',tablekey_value)
                    # pprint.pprint(results)
                else:
                    # print('New table entry added for ', tablekey_value)
                    results[tablekey_value] = 1

            process_counter += 1
            # Update print-in-place
            print(f'Processing: {process_counter} of {total_fits_files}   Matched: {matched_files_counter}', end='\r' )
            #pprint.pprint(results)
                
        # Print the final results
        print('\n')
        total_filtered_files=0
        print(tablekey,'\t Count')
        print('======================')
        for firstcol, count in sorted(results.items(), key=lambda items:items[1], reverse=True):
            print(f"{firstcol:<14} {count:>7}")
            total_filtered_files += count
        print('======================')
        total_label='Total'
        print(f"{total_label:<14} {total_filtered_files:>7}")
        #print('\n')
        if (have_keywords):
            print('* that match these keywords: ')
            for row in keyword_ndarray:
                print(row[0],row[1],row[2])
        print('\n')


    else:
        # iterate across all the files, then iterate across all the keywords, testing if the headers match all keywords
        # if so, do the appropriate move or copy of the matching files to the outputdir.
        # What if no Move or Copy is set????
        print('No tablekey, doing move/copy')
    

    exit(0)


    # EXAMPLE LOOP ONLY!!!
    print("Starting Loop!")
    # Loop through all FITS files in the source directory
    for filename in os.listdir(inputdir):
        #print('Filename is: ',{filename})
        if filename.lower().endswith('.fits') or filename.lower().endswith('.fit'):
            filepath = os.path.join(inputdir, filename)
            #print({filepath})
            try:
                header = fits.getheader(filepath, hdu)
                #print(header)
                eqmode = header.get('EQMODE', None)
                thisscope = header['TELESCOP']
                print('Scope is: ',{thisscope}, ' EQ Mode is: ', {eqmode})
                if eqmode == 1:
                    #shutil.copy(filepath, os.path.join(destination_dir, filename))
                    print(f"Found {filename} (EQMODE=1)")
                #else:
                    #print(f"Skipped {filename} (EQMODE={eqmode})")
                #hdul.close()
            except Exception as e:
                print(f"Error reading {filename}: {e}")

    exit(0)


    # Loop through all FITS files in the source directory
    # for filename in os.listdir(inputdir):
    #     if filename.lower().endswith('.fits') or filename.lower().endswith('.fit'):
    #         filepath = os.path.join(inputdir, filename)
    #         try:
    #             with fits.open(filepath, mode='readonly') as hdul:
    #                 header = hdul[0].header
    #                 #print(header)
    #                 eqmode = header.get('EQMODE', None)
    #                 thisscope = header['TELESCOP']
    #                 #print('Scope is: ',{thisscope}, ' EQ Mode is: ', {eqmode})
    #                 if eqmode == 1:
    #                     #shutil.copy(filepath, os.path.join(destination_dir, filename))
    #                     print(f"Found {filename} (EQMODE=1)")
    #                 #else:
    #                     #print(f"Skipped {filename} (EQMODE={eqmode})")
    #                 hdul.close()
    #         except Exception as e:
    #             print(f"Error reading {filename}: {e}")

def get_arguments():
    global inputdir
    global keyword_ndarray
    global outputdir
    global tablekey
    global args
    global movefiles
    global copyfiles

    parser=ArgumentParser(
        description='This program lets you select FITS files based on specific Header keywords and their values. You can then do statistics, move or copy the files'
        )

    parser.add_argument("-i", "--inputdir", help='The directory containing your FITS files', dest='inputdir', default='.')
    parser.add_argument("-o", "--outputdir", help='The directory where matching FITS files will be moved/copied', dest='outputdir', default="./output")
    parser.add_argument("-k", "--keywords",  metavar='KEY=VALUE', nargs='+', 
                        help='FITS Header Keyword and value in the format -k EQMODE=1. You can use multiple \
                        -k FITSKEY=VALUE keys if you want all to be true to select a file. \
                        You may also ')
    parser.add_argument("-t, --tablekey", dest='tablekey', help='Create a summary table with a count of all the files for each unique value of the "tablekey". \
                        NOTE: if used with -k, all the files must also match all keywords.')
    MoveOrCopy=parser.add_mutually_exclusive_group()
    MoveOrCopy.add_argument("-m, --move", action='store_true', dest='movefiles', help='Move files that match all keywords into the outputdir.')
    MoveOrCopy.add_argument("-c, --copy", action='store_true', dest='copyfiles', help='Copy files that match all keywords into the outputdir.')
    parser.add_argument("-d, --debug", dest="DEBUG", action='store_true', help='Turn on DEBUG printing to console')

    args = parser.parse_args()

 

    inputdir=args.inputdir
    # print('Keywords is:')
    # pprint.pprint(args.keywords)
    if args.keywords is not None:
        # print('Keywords is not None')
        keyword_ndarray=parse_keywords(args.keywords)
    #print('Keyword_ndarray is: \n', keyword_ndarray)
    outputdir=args.outputdir
    tablekey=args.tablekey

    #print('Source directory: ',{inputdir},'- Keyword/Value Pairs: ',keyword_ndarray)
    return(args)



def parse_keywords(entries):

   # Regular expression to capture keyword, operator and value
    pattern = r'(\w+)(<=|>=|!=|=|<|>)(.+)'

    processed = []
    for entry in entries:
        match = re.match(pattern, entry)
        if match:
            keyword, operator, value = match.groups()
            value = value.strip() # clean up any leading or trailing spaces

            # Convert bare '=' to '=='
            if operator == '=':
                operator = '=='

            # do any type conversion inferred from operator type    
            if operator in ('<', '>', '<=', '>='):
                value = float(value)
            elif operator in ('==', '!='):
                # could be numeric or string; try to parse number
                try:
                    if '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass  # keep as string
            processed.append([keyword, operator, value])

    # Convert to numpy array with dtype=object
    ndarray = np.array(processed, dtype=object)
    return ndarray

#OLD CODE BELOW--PROBABLY DELETE.

# def parse_single_key_value_pair(s):
#     # Parse a single key/value pair which are separated by an =
    
#     items = s.split('=')
#     #print(items)
#     key = items[0].strip() # remove blanks on keys
#     if len(items) > 1: #i.e. there is a value included
#         # rejoin the value:
#         value = '='.join(items[1:])
#     else: # Set value to wildcard if there's a Keyword listed without an = or if value is ''
#         value='*'
#     if value=='': 
#         value='*'
#     return (key, value)

# def parse_keywordsX(items):
#     """
#     Parse a series of key-value pairs and return a dictionary
#     """
#     d = {}

#     if items:
#         for item in items:
#             key, value = parse_single_key_value_pair(item)
#             d[key] = value
#     return d

# OLD CODE ABOVE, PROBABLY DELETE # 

def signal_handler(sig, frame): # handle CTRL+C
    print('\nSIGINT or CTRL+C detected.  Exiting.')
    sys.exit(0)


def count_files_with_glob(directory_path, pattern):
    """
    Counts files in a directory matching a pattern using glob.

    Args:
        directory_path (str): The path to the directory to search.
        pattern (str): The pattern to match (e.g., "*.txt", "image_*.png").

    Returns:
        int: The number of matching files.
    """
    full_pattern = os.path.join(directory_path, pattern)
    matching_files = glob.glob(full_pattern)
    return len(matching_files)

if __name__ == '__main__':
    raise SystemExit(main())