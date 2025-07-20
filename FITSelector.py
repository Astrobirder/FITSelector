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
import shutil
import plotly.express as px
import plotly.graph_objects as go


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
log=False
logfile='FITSelector.log'   # Log file for debug output
results=defaultdict(int)    # use a defaultdict for the results table, with tablekey values and counts
plotmap=False               # Flag to indicate if we should plot the results on a map
column_names=['Latitude', 'Longitude', 'Count']  # Column names for the DataFrame used for plotting results on a map
mymap=pd.DataFrame(columns=column_names)        # DataFrame for plotting results on a map


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
    global total_fits_files
    global results
    global log
    global logfile
    global plotmap
    global mymap

    signal.signal(signal.SIGINT, signal_handler) #Handle CTRL+C
    args=get_arguments()
    # print('Args returned are: \n')
    # pprint.pprint(args)
    # print('\n')
    log=args.log
    copyfiles=args.copyfiles
    movefiles=args.movefiles
    plotmap=args.plot

    if log:
        print('Logging debug info to ', logfile)
        logfile=os.path.join(inputdir, logfile)  # set the logfile to be in the output directory
        # Set up logging to a file in the inputdir
        import logging
        logging.basicConfig(filename=os.path.join(inputdir, logfile), filemode='a', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        logging.info('\n\n Starting FITSelector with args:\n %s', args)

    # Verify that arguments are acceptable
    if (inputdir is None):
        inputdir = "."
        print('No input directory specified, using current directory: ', inputdir)
        logging.info('No input directory specified, using current directory: %s', inputdir)  if log else None
    if not os.path.isdir(inputdir):
        print(f'Requested Input Directory ',{inputdir},' does not exist.')
        logging.error('Requested Input Directory %s does not exist.', inputdir) if log else None
        exit(1)

    if ( (copyfiles or movefiles) and (outputdir is None) ):
        outputdir = os.path.join(inputdir, 'output')  # if no outputdir specified, use inputdir/output
        print('No output directory specified with move/copy!')
        logging.info('No output directory specified with move/copy! Using %s', outputdir) if log else None
        print (f'Creating & using {outputdir}')
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

    total_fits_files=count_files_with_glob(inputdir, '*.fit*')

    print('Processing FITS files in ',inputdir)
    logging.info('Processing FITS files in %s', inputdir) if log else None
    print('Total FITS files: ', total_fits_files, '\n')
    logging.info('Total FITS files: %d', total_fits_files) if log else None

    if (total_fits_files<1):
        print('Directory doesn\'t contain any valid files')
        print('Please use a directory with either *.fit or *.fits files')
        exit(0)

# Zero out counters before we start looping through the files.
    process_counter = 0
    matched_files_counter = 0

# Setup the results table based on whether we have a ALLHEADERS tablekey or not.
    if (tablekey != 'ALLHEADERS'):
        results= defaultdict(lambda: {'count': 1}) # use a single column defaultdict to store file counts for each unique tablekey value  
    else:
        results=  defaultdict (lambda: {'count': 0,  'example':  None,  'comment':  None})
        # use a 3 column defaultdict where key is each unique FITS Header, the first column has file counts 
        # for each unique header, and then columns for the example value and comment from the header 
        # of the file that first contained that particular header.

# Main file loop
    for filename in os.listdir(inputdir):
        #print('Filename is: ',{filename})
        #user_input = input("Hit <enter> to process")

# Filter for only .fit or .fits files.
        if filename.lower().endswith('.fits') or filename.lower().endswith('.fit'):

            filepath = os.path.join(inputdir, filename)
            # print({filepath})
            logging.info('Examining file: %s', filepath) if log else None
            try:
                header = fits.getheader(filepath, hdu)
                #print('Got Header from: ', filepath)
                #pprint.pprint(header)
                #print('\n\n')                   
            except Exception as e:
                print(f"Error reading {filename}: {e}, skipping it")
                logging.error('Error reading %s: %s', {filename}, 'Skipping') if log else None
                continue
        else:
            print('Skipping file ', filename, ' as it is not a FITS file.')
            logging.info('Skipping file %s as it is not a FITS file.', filename) if log else None
            continue

        #print('Processing file: ',{filepath})
        # Check if we have keywords and loop goes through them.
        # Set all_keywords_match if they do
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
                    if op_function(header[row[0]],row[2]):
                        keyword_match_count += 1

            if (keyword_match_count==keyword_rows):
                all_keywords_match = True
                #print('We did!')
                # All the keywords have matched!
                # Check to see if the tablekey exists in the header of the file under test
    
        # If every keyword matches or we aren't using keywords, increment the matched files counter
        if(all_keywords_match or not have_keywords):
        # everythings good with this file, add it to the table, then 
            matched_files_counter += 1

            # If we have a tableky, add the file's info to the results table.
            if tablekey:
            # Main option Tablekey processing--generating summary table.
                #print('We have a tablekey! and it is: ', tablekey)
                # We have a tablekey, so generate a summary table based on differing tablekeys across all the files.

                if tablekey == 'ALLHEADERS':
                    #print('We have ALLHEADERS as the tablekey, so we will create a list of all headers in all files')
                    # If the tablekey is ALLHEADERS, we want to create a list of all the headers in all of the files
                    # that are in inputdir that match all (optional) keywords, with a count of the files that contain each keyword.
                    for key in header:
                        if key not in results:
                            this_example = header[key] if key in header else 'None'
                            this_comment = header.comments[key] if key in header.comments else 'No comment'

                            results[key] = {'count': 1, 'example': header[key], 'comment': header.comments[key]}
                        else:
                            results[key]['count'] += 1

                else: # regular tablekey processing              
                    # Handle any file that doesn't have a tablekey in the header
                    if (tablekey) not in header:
                        #print (tablekey, ' is not in file ', filename, 'Setting it to None' )
                        tablekey_value='None'
                    if (tablekey) in header:  # if the key exists, pull the value out of the header
                        tablekey_value=header[tablekey]

                    # See if the entry already exists, and if so, increment the count, otherwise, create it and set count to 1
                    if (tablekey_value) in results:
                        results[tablekey_value] += 1
                        # print('Incrementing ',tablekey_value)
                        # pprint.pprint(results)
                    else:
                        # print('New table entry added for ', tablekey_value)
                        results[tablekey_value] = 1

            if plotmap:
                # If we are plotting the results on a map, we need to extract the RA and DEC from the header
                if 'SITELAT' in header and 'SITELONG' in header:
                    latitude = round(header['SITELAT'],2)
                    longitude = round(header['SITELONG'],2)
                    # update the mymap DataFrame appropriately.
                    mask = (mymap['Latitude'] == latitude) & (mymap['Longitude'] == longitude)
                    if mask.any(): # If the latitude and longitude already exist in the DataFrame, just increment the count
                        mymap.loc[mask, 'Count'] += 1
                    else: # New latitude and longitude, add a new row with count=1
                        new_row = {'Latitude': latitude, 'Longitude': longitude, 'Count': 1}
                        mymap = pd.concat([mymap, pd.DataFrame([new_row])], ignore_index=True)
                    
                else:
                    print(f"Skipping {filename} as it does not have SITELAT and SITELONG keywords.")
                    logging.warning('Skipping %s as it does not have SITELAT and SITELONG keywords.', filename) if log else None
            
            if copyfiles:
                #print('Copying file ', filename, ' to ', outputdir)
                # Copy the file to the output directory
                try:
                    shutil.copy(filepath, os.path.join(outputdir, filename))
                except Exception as e:
                    print(f"Error copying {filename}: {e}")
            
            if movefiles:
                #print('Moving file ', filename, ' to ', outputdir)
                # Move the file to the output directory
                try:
                    shutil.move(filepath, os.path.join(outputdir, filename))
                except Exception as e:
                    print(f"Error moving {filename}: {e}")

        # Completed processing of this file, increment the process counter, print results.
        process_counter += 1
        # Update print-in-place
        print(f'Processing: {process_counter} of {total_fits_files} FITS files - Matched: {matched_files_counter}', end='\r' )
        #pprint.pprint(results)

    # End of main loop through files.

    # Print the final results table if we have a tablekey
    if (tablekey):
        print('\n')
        total_filtered_files=0
        if (tablekey != 'ALLHEADERS'):
            print(tablekey,'\t Count')
            print('======================')
            for firstcol, count in sorted(results.items(), key=lambda items:items[1], reverse=True):
                print(f"{firstcol:<14} {count:>7}")
                total_filtered_files += count
            print('======================')
            total_label='Total'
            print(f"{total_label:<14} {total_filtered_files:>7}")
            print('\n')
            if (have_keywords):
                print('* that match these keywords: ')
                for row in keyword_ndarray:
                    print(row[0],row[1],row[2])
            print('\n')
        else:
            print('All Headers in all files that match the keywords:')
            print('\n')
            print('Header     Count                 Typ. Value   Comment')
            print('======================================++=============')
            for header_name, data in sorted(results.items(), key=lambda items:items[1]['count'], reverse=True):
                count = data['count']
                example = str(data['example'])
                if header_name == 'COMMENT':
                    comment = data['example']  
                    example = ""
                    # COMMENT key is a special case, use example as comment
                else:
                    comment = data['comment']
                print(f"{header_name:<10} {count:>5} {example:>26}   {comment}")
                total_filtered_files += count
            print('======================================++=============')
            total_label='Total'
            print(f"{total_label:<10} {total_filtered_files:>5}")
            print('\n')
    
    if plotmap:
        mymap = mymap.astype({'Latitude': float, 'Longitude': float, 'Count': int})  # Ensure correct data types
        # pprint.pprint(mymap)
        # print('\n\n')
        # print('Here\'s the unique pairs')
        # print(mymap[['Latitude', 'Longitude']].drop_duplicates())
        # print('\n Datatypes are:')
        # print(mymap.dtypes)
        # mymap['Latitude'] = mymap['Latitude'].astype(float)
        # mymap['Longitude'] = mymap['Longitude'].astype(float)
        # mymap['Count'] = mymap['Count'].astype(int)
        # print(mymap.dtypes)

        fig = go.Figure()
        fig.add_trace(go.Scattergeo(lat=mymap['Latitude'], lon=mymap['Longitude'], text=mymap['Count'],
                                   marker=dict( color=mymap['Count'], showscale=True, colorbar=dict(title='Count')),
                                   mode='markers', name='FITS Frames'))
        
   
        fig.update_layout(
            title='FITS Frames Count by Location',
            geo=dict(
                scope='world',
                showland=True,
                projection_type='natural earth',
                landcolor='rgb(217, 217, 217)',
                showlakes=True,
                lakecolor='rgb(255, 255, 255)',
            )
        )

        fig.show()
        #print(mymap.to_dict())

    exit(0)


def get_arguments():
    global inputdir
    global keyword_ndarray
    global outputdir
    global tablekey
    global args
    global movefiles
    global copyfiles
    global log
    global logfile


    parser=ArgumentParser(
        description='This program lets you select FITS files based on specific FITS Header keywords and their values. You can then do statistics, move or copy the files'
        )

    parser.add_argument("-i", "--inputdir", help='The directory containing your FITS files', dest='inputdir')
    parser.add_argument("-o", "--outputdir", help='The directory where matching FITS files will be moved/copied. If unspecified, inputdir/output will be used', dest='outputdir')
    parser.add_argument("-k", "--keywords",  metavar='KEY=VALUE', nargs='+', 
                        help='FITS Header Keyword and value in the format -k \"EQMODE=1\". You can use multiple \
                        -k \"FITSKEY1=VALUE1 FITSKEY2=VALUE2\" keys if you want all to be true to select a file.')
    parser.add_argument("-t, --tablekey", dest='tablekey', help='Create a summary table with a count of all the files for each unique value of the "tablekey". \n \
                        NOTE: if used with -k, all the files must also match all keywords. \n \
                        If you use ALLHEADERS as the tablekey, a list of all the headers in all of the files \n \
                        that are in inputdir that match all (optional) keywords, with a count of the files that \n \
                        contain each keyword, will be printed.')
    parser.add_argument("-p, --plot", action='store_true', dest='plot', help='Plot the count of matching files on a map.')
    MoveOrCopy=parser.add_mutually_exclusive_group()
    MoveOrCopy.add_argument("-m, --move", action='store_true', dest='movefiles', help='Move files that match all keywords into the outputdir.')
    MoveOrCopy.add_argument("-c, --copy", action='store_true', dest='copyfiles', help='Copy files that match all keywords into the outputdir.')
    parser.add_argument("-l, --log", dest="log", action='store_true', help='Log debug info to inputdir/FITSelector.log')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

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
    log

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