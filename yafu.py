# Yafu Factoring module

import subprocess
import os
import time
import re
from operator import mul

class Yafu:
    """Connection to Yafu"""

    """Statics"""
    __theyafu = object()

    """Constants"""
    __PATH_OF_YAFU = "./yafu/yafu-1.34.3/yafu"
    __PATH_OF_YAFU_JOB_FOLDER = "./yafu_tmp/"

    def __init__(self, token):
        if token is not self.__theyafu:
            raise ValueError("don't construct directly, use make_foo")

    @classmethod
    def getYafu(cls):
        return cls(cls.__theyafu)
    
    ##########################################################################
    ## Real Class Functions ##################################################
    ##########################################################################

    def factor_to_lib(self,n,num_threads=1, yafu_additional_args=[]):
        """
        Factors n with Yafu and writes result to library.
        Does nothing if n is already in library.
        
        n: long to be factored
        """
        n = long(n)

        if thelib.from_lib(n, update_lib=True) != None:
            return n,thelib.from_lib(n)

        if not os.path.exists(self.__PATH_OF_YAFU_JOB_FOLDER):
            os.mkdir(self.__PATH_OF_YAFU_JOB_FOLDER)
        tstmp = str(long(time.time()*1E6))

        os.mkdir(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp)
        
        with open(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+"/job.bat","w") as f:
            f.write("factor("+str(n)+")\r\n")
        f.close()

        os.symlink(os.path.abspath(self.__PATH_OF_YAFU+".ini"),\
                self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+"/yafu.ini")

        #try:
        p = subprocess.Popen([os.path.abspath(self.__PATH_OF_YAFU), \
            "-batchfile", "job.bat",\
            #"-logfile"\
            #, os.path.abspath(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+".log"),\
            #"-session"\
            #, os.path.abspath(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+".session.log"),\
            #"-qssave", os.path.abspath(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+".dat"),\
            "-of", "out.txt",\
            "-threads",str(num_threads)]+yafu_additional_args,\
            cwd=self.__PATH_OF_YAFU_JOB_FOLDER+tstmp)
        p.wait()
        if p.returncode != 0:
            # if yafu fails try deep ecm
            return factor_to_lib(self,n,num_threads,\
                    yafu_additional_args=["-plan","deep"])

        factorization = []
        with open(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp+"/out.txt") as f:
            for l in f:
                match = re.search("\((\d+)\)(/.*)$",l)
                if match:
                    file_n = long(match.groups()[0])
                    if file_n != n:
                        raise ValueError("Not the same n!")
                    for pr in re.findall("/[\d\^]+",match.groups()[1]):
                        split = pr[1:].split("^")
                        if len(split) > 1:
                            factorization += \
                                [(long(split[0]),long(split[1]))]
                        else:
                            factorization += [(long(split[0]),1)]
        f.close()
        if n != reduce(mul,map(lambda pr: pr[0]**pr[1], factorization)):
            raise ValueError("factorization "+str(factorization)\
                    +" is not a fac of n = "+str(n))

        self.__purge(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp,".*")
        os.rmdir(self.__PATH_OF_YAFU_JOB_FOLDER+tstmp)
        thelib.to_lib(n,factorization)

        return n,factorization


    def batch_to_library(self,batchfile,num_threads=1,yafu_additional_args=[]):
        """factors every line in given batchfile

           batchfile: file containing a number to factor in every line
        """

        with open(batchfile) as fin:
            for phi in fin:
                self.factor_to_lib(phi,num_threads,yafu_additional_args)
        fin.close()



    def __purge(self,dirpath, pattern):
        """Purges a dir with given pattern"""
        for f in os.listdir(dirpath):
            if re.search(pattern, f):
                os.remove(os.path.join(dirpath, f))



class FactoringLibrary:
    """Library that stores already factored numbers"""

    """static objects"""
    _THE_MAGIC_WORD = object()
    __facdict = dict()

    """constants"""
    __MAX_GOOD_LEN = 20
    __PATH_OF_LIB_FILE = "./factor_lib.txt"

    # Make private init function
    def __init__(self, token):
        if token is not self._THE_MAGIC_WORD:
            raise ValueError("don't construct directly, use make_foo")

    @classmethod
    def getLib(cls):
        return cls(cls._THE_MAGIC_WORD)


    ##########################################################################
    ## Real Class Functions ##################################################
    ##########################################################################

    def is_fac_easy(self,n):
        """
        Checks if n can be factored easily

        n: Number to be factored
                    if p == 1: continue
                    if p in prime_factors:
                        idx = prime_factors.index(p)
                        factorization[idx] = \
                                (p,factorization[idx][1]+1)
                    else:
                        prime_factors += [p]
                        factorization += [(p,1)]
        """
        if len(str(long(n))) <= self.__MAX_GOOD_LEN:
            return True
        return False

    def to_lib(self,n,factorization):
        """
        Adds n to the factorization dictionary

        n: Number
        factorzation: the complete factorization of n
        """
        n = long(n)
        # first check factorization
        if n != reduce(mul,map(lambda pr: pr[0]**pr[1], factorization)):
            raise ValueError("factorization "+str(factorization)\
                    +" is not a fac of n = "+str(n))

        self.__facdict[n] = factorization

        with open(self.__PATH_OF_LIB_FILE,"a") as f:
            f.write(str(n)+"\t"+str(factorization)+"\r\n")
        f.close()

    def from_lib(self,n,update_lib=False):
        """
        Reads factorization of n from Library if there

        n: Number
        update_lib: updates library if n first not found

        Returns a factorization if there, None if no was found
        """
        n = long(n)
        if update_lib and not self.__facdict.has_key(n):
            self.update_lib()
        if self.__facdict.has_key(n):
            return self.__facdict[n]
        return None

    def update_lib(self):
        """
        Updates the FactoringLibrary, i.e. reads the file again
        """
        with open(self.__PATH_OF_LIB_FILE) as f:
            for l in f:
                split = l.split("\t")
                self.__facdict[long(split[0])] = eval(split[1].rstrip())
        f.close()


    def clean_lib(self):
        """
        Cleans the Factoringlibrary, i.e. writes the file again
        """
        self.update_lib()
        print("number of entries: "+str(len(self.__facdict.keys())))
        with open(self.__PATH_OF_LIB_FILE,'w') as f:
            for n in self.__facdict.keys():
                f.write(str(n)+"\t"+str(self.__facdict[n])+"\r\n")
        f.close()

    def merge_lib(self,other_file):
        """
        Merges other_file to the FactoringLibrary
        !!! other_file must have same format as specified here !!!
        """
        with open(other_file) as fin_other:
            with open(self.__PATH_OF_LIB_FILE,'a') as fout:
                fout.write(fin_other.read())
            fout.close()
        fin_other.close()
        self.clean_lib()




# use faclib as Factoringlibrary from the outside
thelib = FactoringLibrary.getLib()
yafu = Yafu.getYafu()
