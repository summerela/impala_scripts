#!/usr/bin/env python
import sys
import cStringIO
import xml.etree.ElementTree as xml

def cleanResult(element):
    result = ""
    if element is not None:
        result = str(element.text).strip()
    return result


def process(val):
    root = xml.fromstring(val)
    sceneID = cleanResult(root.find('sceneID'))
    cc = cleanResult(root.find('cloudCover'))
    returnval = ("%s,%s") % (sceneID,cc)
    return returnval.strip()

if __name__ == '__main__':
    buff = None
    intext = False
    for line in sys.stdin:
        line = line.strip()
        if '<metaData>' in line:
            intext = True
            buff = cStringIO.StringIO()
            buff.write(line)
        elif '</metaData>' in line:
            intext = False
            buff.write(line)
            val = buff.getvalue()
            buff.close()
            buff = None
            print process(val)
        else:
            if intext:
                buff.write(line)