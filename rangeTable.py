# Range Table module

import subprocess
import os
import time
import datetime
import re
from operator import mul

class RangeTable:
    """The RangeTable"""

    """Statics"""
    __therangetable = object()
    __table = dict()

    """Constants"""
    __PATH_OF_TABLE = "rangeTable.txt"

    def __init__(self, token):
        if token is not self.__therangetable:
            raise ValueError("don't construct directly, use make_foo")
        self.updateTable()

    @classmethod
    def getRangeTable(cls):
        return cls(cls.__therangetable)
    
    ##########################################################################
    ## Real Class Functions ##################################################
    ##########################################################################

    def isInTable(self, p, n):
        return ( (p,n) in self.__table )


    def updateTable(self):
        """
        Updates the table, i.e. reads the file again
        """
        with open(self.__PATH_OF_TABLE) as f:
            for l in f:
                split = l.split("\t")
                self.__table[eval(split[0])] = split[1].rstrip()
        f.close()


    def toTable(self, p,n, poly ):
        """
        Adds (p,n) poly to the table and updates the file
        """
        self.__table[ (p,n) ] = poly

        with open(self.__PATH_OF_TABLE,"a") as f:
            f.write("(%d,%d)\t%s\n"%(p,n,str(poly)))
        f.close()


theRangeTable = RangeTable.getRangeTable()
