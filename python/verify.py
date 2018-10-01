import sys
import os.path

# Check TT
try:
    import tt
except ImportError:
    pass
    print('ttpy toolbox is not found')
    print('Please install it using pip/pip3')
    print('See https://github.com/oseledets/ttpy for more details')
    exit()

# Check tt_irt
try:
    import tt_irt
except ImportError:
    pass
    print('tt_irt module is not found')
    print('Please change to directory tt_irt_py and run')
    print('  python setup.py install [--user]')
    exit()

# Check UWerr
if (not os.path.exists('puwr.py')):
    print('puwr.py is not found. Downloading...')
    if sys.version_info>(3,0):
        import urllib.request
        url = 'https://raw.githubusercontent.com/dhesse/py-uwerr/master/puwr.py'
        urllib.request.urlretrieve(url, 'puwr.py')
    else:
        import urllib2
        filedata = urllib2.urlopen('https://raw.githubusercontent.com/dhesse/py-uwerr/master/puwr.py')
        datatowrite = filedata.read()
        with open('puwr.py', 'wb') as f:
            f.write(datatowrite)


print('Ready to go!')
