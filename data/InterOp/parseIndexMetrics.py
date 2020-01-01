# BCL to FASTQ file converter
# Copyright (c) 2007-2017 Illumina, Inc.
#
# This software is covered by the accompanying EULA
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein.
#
# \file parseIndexMetrics.py
#
# \brief Parsing of binary SAV-compatible metrics.
#
# \author Peter Saffrey
# \author Mauricio Varea
#

import sys
import struct
import csv

row_order = [ "sample_name", "project_name", "index_name", "cluster_count", "cluster_proportion" ]

class Lane(object):
    def __init__(self, lanenum):
        self.name = lanenum
        self._indexdata = {}

    def addClusterCount(self, index, value):
        self._indexdata[index]["cluster_count"] += value

    def findIndex(self, index):
        return index in self._indexdata

    def createIndex(self, index):
        self._indexdata[index] = {}
        self._indexdata[index]["index_name"] = index
        self._indexdata[index]["cluster_count"] = 0

    def addIndexValue(self, index, key, value):
        self._indexdata[index][key] = value

    def computeClusterPercentage(self):
        clustercounts = [ self._indexdata[index]["cluster_count"] for index in self._indexdata ]
        totalclusters = float(sum(clustercounts))
        for index in self._indexdata:
            cluster_proportion = float(self._indexdata[index]["cluster_count"]) / totalclusters
            self.addIndexValue(index, "cluster_proportion", cluster_proportion)

    def dump(self, outfh):
        writer = csv.writer(outfh, delimiter="\t")
        for index in self._indexdata:
            row = [ self.name ]
            for entry in row_order:
                row.append(self._indexdata[index][entry])
            writer.writerow(row)

def dumpHeader(outfh):
    writer = csv.writer(outfh, delimiter="\t")
    writer.writerow([ "lane_number" ] + row_order)

def getInt(fin, format, length):
    return struct.unpack(format, f.read(length))[0]

def getStr(fin, length):
    fmt = "{}s".format(length)
    return struct.unpack(fmt, f.read(length))[0]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "need an IndexMetricsOut.bin file"
        sys.exit(1)
    lanes = {}
    binfile = sys.argv[1]
    with open(binfile, "rb") as f:
        version = f.read(1)
        #print "version: %s" % ord(version)
        while f:
            # kind of hacky way of dropping out of the loop
            try:
                lanenum = getInt(f, 'H', 2)
                if lanenum not in lanes:
                    lanes[lanenum] = Lane(lanenum)
            except:
                break
            thislane = lanes[lanenum]
            tilenum = getInt(f, 'H', 2)
            readnum = getInt(f, 'H', 2)
            #print "...processing: lane %d, tile %d, read %d" % (lanenum,tilenum,readnum)
            index_length = getInt(f, 'H', 2)
            index_name = getStr(f, index_length)
            if not thislane.findIndex(index_name):
                thislane.createIndex(index_name)
            thislane.addClusterCount(index_name, getInt(f, 'I', 4))
            name_length = getInt(f, 'H', 2)
            thislane.addIndexValue(index_name, "sample_name", getStr(f, name_length))
            project_length = getInt(f, 'H', 2)
            thislane.addIndexValue(index_name, "project_name", getStr(f, project_length))
    dumpHeader(sys.stdout)
    for lanenum in lanes:
        lane = lanes[lanenum]
        lane.computeClusterPercentage()
        lane.dump(sys.stdout)
