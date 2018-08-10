'''
Created on Dec 12, 2016

@author: uvilla
'''
import json
import urllib
import ssl

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

#key = 'AIzaSyCnCjlxS5LUvP1jjXansnk1ymt-Liga2aQ'
#key = 'AIzaSyBLuHKVMU8kxU5R2ykwEgGTdhuINg-3HUE'
#key = ' AIzaSyBCcGtouSRX4uLUyzujx_DMRcLDmcKWGT0'
key = 'AIzaSyD-SX9PQc3F7B88JGyJr9SzDjgjwOUgxGo'
versions = [('v2.1.0', 'https://goo.gl/N7g6wU'),
            ('v2.0.0', 'https://goo.gl/hYazvu'),
            ('v1.6.0', 'https://goo.gl/FsoZsG'),
            ('v1.5.0', 'https://goo.gl/doJZRB'),
            ('v1.4.0', 'https://goo.gl/37bskk'),
            ('v1.3.0', 'https://goo.gl/NgJ887'),
            ('v1.2.0', 'https://goo.gl/OcvROZ'),
            ('v1.1.0', 'https://goo.gl/pDb10B'),
            ('v1.0.2', 'https://goo.gl/wMb18C'),
            ('v1.0.1', 'https://goo.gl/srnvsc'),
            ('v1.0.0', 'https://goo.gl/h7HUFJ')
           ]

request = 'https://www.googleapis.com/urlshortener/v1/url?shortUrl={0}&projection=FULL&key={1}'

print '{0:7} {1:9}'.format('version', 'downloads')
total_downloads = 0
for v in versions:
    filehandle = urllib.urlopen(request.format(v[1], key), context=ctx)
    data = json.load(filehandle)
    try:
    	total_downloads += int(data['analytics']['allTime']['shortUrlClicks'])
    	print '{0:7} {1:9}'.format(v[0], int(data['analytics']['allTime']['shortUrlClicks']))
    except KeyError:
        print data

print '{0:7} {1:9}'.format('all', total_downloads)
