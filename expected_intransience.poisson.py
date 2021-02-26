#! /usr/bin/python3

__author__ = "Patrick Denis Browne"
__date__ = "02/2021"
__email__ = "pdbr@plen.ku.dk"
__license__ = "GPLv3"

import re
import sys
import copy
import random
import warnings
import numpy as np
from numpy import average
from scipy.stats import poisson


class otuTable:
    def __init__(self, otufile, taxonomy=None):
        '''
        Parses a text otu table. Comment lines (beginning with '#') are
        allowed in the first lines only. The last of the first
        consecutive lines beginning with '#' will be assumed to be the
        header of the table. No further comment lines will be parsed
        appropriately. The header is assumed to contain sample names. If
        taxonomy information is included, it must be specified at
        instantiation with the relevant column's header (eg.
        taxonomy='Taxonomy'). The taxonomy assignments are assumed to
        conform to the following formats:
        k__<K>;p__<P>;c__<C>;... etc
        partial assignments like:
        K__<K>;P__<P>;c__;o__;f__,g__;s__
        All numbers in the table can be integer counts or relative
        frequencies
        '''
        if isinstance(otufile, str):
            self.matrix = [line.rstrip('\n').split('\t') for line in \
            open(otufile, 'rt').readlines()]
        else:
            # The file is assumed to be a file-handle opened for
            # reading in text mode with the cursor at the appropriate
            # position (usually you want to do `otufile.seek(0)' before
            # passing it into here
            self.matrix = [line.rstrip('\n').split('\t') for line in \
                           otufile.readlines()]
        
        # Allow Comments above the start of the OTU table only
        if self.matrix[0][0].startswith('#'):
	        while self.matrix[1][0].startswith('#'):
	            a = self.matrix.pop(0)
        
        # assert self.matrix[0][0] == '#OTU ID',\
        # "Improperly formatted OTU table encountered. The first line should begin with '#OTU ID'"
        
        assert len(list(set([row[0] for row in self.matrix]))) == \
               len([row[0] for row in self.matrix]), \
               "Duplicate OTUs or blank line(s) were/was found in the OTU table"
        
        self.sampleIDs = self.matrix.pop(0)
        _ = self.sampleIDs.pop(0)
        self.otuIDs = [row.pop(0) for row in self.matrix]
        if taxonomy != None:
            taxonomyColIndex = self.sampleIDs.index(taxonomy)
            self.sampleIDs.pop(taxonomyColIndex)
            self.lineages = [row.pop(taxonomyColIndex) for row in self.matrix]
        else:
            self.lineages = None
        self.matrix = [[float(col) for col in row] for row in self.matrix]
        self.updateCounts()
        self.backup()
    
    def numSamples(self):
        return len(self.sampleIDs)
    
    def numOtus(self):
        return len(self.otuIDs)
    
    def backup(self):
        '''
        Creates a backup of the otuTable object that can be restored
        to undo some modificantions. This method is called upon
        instantiation and this original backup can be over-written at
        any point in the future.
        '''
        self.__originalmatrix__ = copy.deepcopy(self.matrix)
        self.__originalsampleIDs__ = copy.copy(self.sampleIDs)
        self.__originalotuIDs__ = copy.copy(self.otuIDs)
        self.__originallineages__ = copy.copy(self.lineages)
    
    def restore(self):
        '''
        Restores a backup of the otuTable object. A backup is created
        upon instantiation, so this method can always be called.
        '''
        self.matrix = copy.copy(self.__originalmatrix__)
        self.sampleIDs = copy.copy(self.__originalsampleIDs__)
        self.otuIDs = copy.copy(self.__originalotuIDs__)
        self.lineages = copy.copy(self.__originallineages__)
        self.updateCounts()
    
    def updateCounts(self):
        '''
        Sums up values for OTUs and values for Samples. Called upon
        instantiation and for the restore method.
        '''
        self.otuCounts = [sum(row) for row in self.matrix]
        self.sampleCounts = [sum([row[i] for row in self.matrix]) for i in
                             range(len(self.matrix[0]))]
    
    def dropOTUbyIndex(self, index):
        '''Deletes an OTU from the object by index. This method is
        primarily intended for calling by other methods of this object
        '''
        assert isinstance(index, int), "'index' must be an integer"
        cut = self.otuIDs.pop(index)
        cut = self.matrix.pop(index)
        if self.lineages:
            cut = self.lineages.pop(index)
        self.updateCounts()
    
    def dropSamplebyIndex(self, index):
        assert isinstance(index, int), "'index' must be an integer"
        cut = self.sampleIDs.pop(index)
        for i in range(len(self.matrix)):
            cut = self.matrix[i].pop(index)
        self.updateCounts()
    
    def rarifyOTUtable(self, depth, keepLow=False):
        '''
        Rarifies the OTU table by random subsampling. This will only run
        on integer counts of OTUs.
        'depth': the integer number of OTUs to sample for each sampleID
        Assigning any value to 'keepLow' will override the default
        behaviour of dropping samples that are less represented than the
        'depth' value
        '''
        ################################################################
        # There's a bug where the otu table values are processed as    #
        # floats while integer counts are required for this process.   #
        ################################################################
        if max(self.getSampleTotalAbundances().values()) <=1:
            raise Exception("The values appear to be relative counts.")
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                self.matrix[i][j] = int(self.matrix[i][j])
        ################################################################
        assert isinstance(depth, int), "'depth' must be an integer"
        # Drop the samples which are too small
        if keepLow == False:
            dropSampleIndices = []
            for i in range(len(self.sampleIDs)):
                if self.sampleCounts[i] < depth:
                    dropSampleIndices.append(i)
            dropSampleIndices.sort() # Should already be sorted?
            dropSampleIndices.reverse()
            for dIndex in dropSampleIndices:
                self.dropSamplebyIndex(dIndex)
        
        # This will not modify the samples that are smaller than the 'depth'
        for i in range(len(self.sampleIDs)):
            if self.sampleCounts[i] >= depth:
                sampleOtus = []
                for j in range(len(self.otuIDs)):
                    otuCountInSample = self.matrix[j][i]
                    if not isinstance(otuCountInSample, int):
                        self.restore()
                        raise Exception("The values in the OTU table must be integers for this method")
                    for k in range(otuCountInSample):
                        sampleOtus.append(self.otuIDs[j])
                random.shuffle(sampleOtus)
                sampleOtus = sampleOtus[:depth]
                for j in range(len(self.otuIDs)):
                    self.matrix[j][i] = sampleOtus.count(self.otuIDs[j])
        self.dropZeroRows()
        self.updateCounts()
    
    def getLineages(self):
        return [lineage for lineage in self.lineages]
    
    def getOTUsLineage(self, otuID):
        '''
        When given an OTU's identifier as input, it returns the
        associated lineage term. If no lineage terms are present, it
        will return None
        '''
        if self.lineages == None:
            return None
        else:
            otuIndex = self.otuIDs.index(otuID)
            return self.lineages[otuIndex]
    
    def getOTUlineages(self, otuList):
        '''
        Given a (nested) list of otu identifiers as input, it returns
        a list of lineages associated with the otu identifiers. The
        returned list will have the same dimensions as the input. If
        no lineage terms were parsed in the instance of self, it will
        still return a list with the same dimensions as the input, where
        every element of the returned list will be None
        '''
        assert isinstance(otuList, list), "A list is required as input"
        output = []
        for elem in otuList:
            if isinstance(elem, str):
                output.append(self.getOTUsLineage(elem))
            elif isinstance(elem, list):
                output.append(self.getOTUlineages(elem))
            else:
                raise(ValueError, "Unexpected data type")
        return output
    
    def keepOnlySamples(self, keeplist):
        '''Drops samples from the otuTable if they are not in the 
        keeplist - a list of sampleIDs to be retained'''
        droplist = [sID for sID in self.sampleIDs if not sID in keeplist]
        dropindices = [self.sampleIDs.index(sampleID) for sampleID in droplist]
        dropindices.sort()
        dropindices.reverse()
        for dropindex in dropindices:
            cut = self.sampleIDs.pop(dropindex)
            for line in self.matrix:
                line.pop(dropindex)
        self.updateCounts()
    
    def dropSamples(self, sampleList):
        '''Drops samples from the OTU table if they are in the list'''
        dropindices = []
        for sampleID in sampleList:
            if sampleID in self.sampleIDs:
                dropindices.append(self.sampleIDs.index(sampleID))
        if len(dropindices) > 0:
            print("Dropping %i samples" % (len(dropindices)))
            dropindices.sort()
            dropindices.reverse()
            for dropindex in dropindices:
                cut = self.sampleIDs.pop(dropindex)
                for line in self.matrix:
                    line.pop(dropindex)
            self.updateCounts()
        else:
            print("Warning: Zero samples were dropped.\nAre the SampleIDs from the correct source?")
    
    def dropZeroRows(self):
        '''Deletes rows were all OTUs have no representation'''
        delIndices = []
        for otu in self.otuIDs:
            otuIndex = self.otuIDs.index(otu)
            abundance = sum(self.matrix[otuIndex])
            if abundance == 0:
                delIndices.append(otuIndex)
        delIndices.sort()
        delIndices.reverse()
        for delIndex in delIndices:
            self.dropOTUbyIndex(delIndex)
    
    def dropOtusByCount(self, mincount, minsamples):
        '''
        Drops any Otus that are not present at least 'mincount' times in
        at least 'minsamples'
        'mincount' = integer
        'minsamples' = integer
        '''
        dropindices = []
        for i in range(len(self.otuIDs)):
            otucount = 0
            for scount in self.matrix[i]:
                if scount >= mincount:
                    otucount += 1
            if otucount < minsamples:
                dropindices.append(i)
        dropindices.reverse()
        print("Dropping %i OTUs" % (len(dropindices)))
        for dropindex in dropindices:
            cut = self.otuIDs.pop(dropindex)
            cut = self.matrix.pop(dropindex)
            if self.lineages:
                cut = self.lineages.pop(dropindex)
        self.updateCounts()
    
    def dropOtusByProportion(self, proportion, minsamples, mincount=None):
        '''
        'proportion' must be a float between 0 and 1
        'minsamples' is the minimum number of samples for which an OTU
        must represent the minimum proportion.
        if 'mincount' is specified, then a sample with a count greater
        than or equal to mincount will not be counted as being below the
        'proportion' argument.
        '''
        assert proportion > 0 and proportion < 1, "'proportion' must be a number between 0 and 1"
        assert isinstance(minsamples,int), "'minsamples' must be an integer"
        self.updateCounts() # Just in case of what?
        dropindices = []
        for i in range(len(self.otuIDs)):
            otucount = 0
            for j in range(len(self.matrix[i])):
                scount = float(self.matrix[i][j])
                if float(scount)/self.sampleCounts[j] >= proportion:
                    otucount += 1
                elif mincount != None:
                    if float(scount) >= mincount:
                        otucount += 1
            if otucount < minsamples:
                dropindices.append(i)
        dropindices.sort()
        dropindices.reverse()
        print("Dropping %i OTUs" % (len(dropindices)))
        for dropindex in dropindices:
            cut = self.otuIDs.pop(dropindex)
            cut = self.matrix.pop(dropindex)
            if self.lineages:
                cut = self.lineages.pop(dropindex)
        self.updateCounts()
    
    def dropOTUs(self, otus, inverse=False):
        '''
        'otus' can be a string (dropping a single OTU) or a list (to
        drop one or more OTUs)
        '''
        if isinstance(otus, str):
            otus = [otus]
        elif isinstance(otus, list):
            otus = list(set(otus))
        else:
            raise Exception("'otus' must be a string or a list")
        if inverse:
            otuIndices = [self.otuIDs.index(otu) for otu in self.getOTUlist()
                          if not otu in otus]
            otuIndices.sort()
            otuIndices.reverse()
            for otuIndex in otuIndices:
                self.dropOTUbyIndex(otuIndex)
        else:
            for otu in otus:
                otuIndex = self.otuIDs.index(otu)
                self.dropOTUbyIndex(otuIndex)
    
    def getOTUlist(self):
        return [otuID for otuID in self.otuIDs]
    
    def getSampleIDs(self):
        return [sampleID for sampleID in self.sampleIDs]
    
    def getSampleAbundances(self, sampID, relative=False):
        sIDindex = self.sampleIDs.index(sampID)
        if relative == False:
            return [self.matrix[i][sIDindex] for i in range(len(self.matrix))]
        else:
            total = sum([self.matrix[i][sIDindex] for i in range(len(self.matrix))])
            return [float(self.matrix[i][sIDindex])/total for
                    i in range(len(self.matrix))]
    
    def getOTUabundance(self, otuID, sampID, relative=False, pseudocount=0):
        '''
        Returns the value for the OTU in the specified sample.
        If the otutable is count data and relative abundacnes are
        required, then set 'relatvie' to True.
        A pseudocount is only valid for relative abundances.
        '''
        if not sampID in self.sampleIDs:
            return np.nan
        else:
            otuIndex = self.otuIDs.index(otuID)
            sampIndex = self.sampleIDs.index(sampID)
            sample_abundance = self.matrix[otuIndex][sampIndex]
            if relative == False:
                return sample_abundance
            else:
                sample_abundance = float(sample_abundance) + float(pseudocount)
                total_abundance = float(self.getSampleCoverage(sampID)) + float(pseudocount)
                return sample_abundance/total_abundance
    
    def getOTUabundances(self, otuID, sampIDs, relative=False, pseudocount=0):
        '''
        otuID must be a string present in self.otuIDs
        sampIDs must be a list of sample identifiers present in
        self.sampleIDs
        returns a list with otu counts / abundances in the same order
        as the sampleIDs present in sampIDs
        set 'relative' to anything other than False to get relative
        abundances instead
        A pseudocount is only applied in the case of relative abundance lookups
        '''
        if isinstance(sampIDs, list):
            return [self.getOTUabundances(otuID, sampID, relative=relative,
                    pseudocount=pseudocount) for sampID in sampIDs]
        else:
            return self.getOTUabundance(otuID, sampIDs, relative=relative,
                                        pseudocount=pseudocount)
    
    def getOTUmatrix(self, otuID, sampIDs, relative=False,
                     ignore=['N/A','#N/A']):
        """
        otuID must be a valid string from self.otuIDs
        
        sampIDs must be a 2-dimensional list of sample identifiers or
        values from the ignore variable
        
        By default, the raw numbers in the original otutable will be
        returned. These are usually counts. If you want relative
        abundances instead, then change the relative flag to True.
        
        Returns a matrix of the OTUs for the otu(s) and sample(s)
        specified
        """
        otuMatrix = copy.deepcopy(sampIDs)
        for i in range(len(otuMatrix)):
            for j in range(len(otuMatrix[i])):
                if not otuMatrix[i][j] in ignore:
                    if relative == False:
                        otuMatrix[i][j] = self.getOTUabundances(otuID,
                            [otuMatrix[i][j]])[0]
                    else:
                        otuMatrix[i][j] = self.getOTUabundances(otuID,
                            [otuMatrix[i][j]], relative=True)[0]
        return otuMatrix
    
    def getSampleTotalAbundances(self, samples=None):
        '''
        'samples': if this is not specified, the abundances of all
                   sampleIDs will be returned
                   specifying a string (one sampleID) or a list of
                   stings (one more sampleIDs) will retun abundances of
                   only those sample IDS
        Returns: dictionary {str(sampleID):int/float(sampleCount)}
        '''
        self.updateCounts()
        if not samples:
            return {self.sampleIDs[i]:self.sampleCounts[i] for
                    i in range(len(self.sampleIDs))}
        else:
            if isinstance(samples, str):
                sampIndex = self.sampleIDs.index(samples)
                return {self.sampleIDs[sampIndex]:self.sampleCounts[sampIndex]}
            elif isinstance(samples, list):
                sampIndices = list(set([self.sampleIDs.index(sample) for
                                        sample in samples]))
                return {self.sampleIDs[sampIndex]:self.sampleCounts[sampIndex]
                        for sampIndex in sampIndices}
    
    def getSampleCoverage(self, sampID):
        """Returns the summed abundances of a sample"""
        sampIndex = self.sampleIDs.index(sampID)
        return sum([self.matrix[i][sampIndex] for i in range(len(self.matrix))])
    
    def getSamplesCoverages(self, sampIDs):
        """Returns the summed abundances for a list of samples"""
        return [self.getSampleCoverage(sampID) for sampID in sampIDs]
    
    def writeOTUtable(self,File):
        '''
        'File' is either a file name (string) or a filehandle. If the
        file already exists, it will be overwritten.
        '''
        if isinstance(File, str):
            File = open(File, 'wt')
        else:
            assert isinstance(File, file), \
            "'File' must either be a file name or a file handle"
        firstline = ['#OTU ID'] + self.sampleIDs
        if not self.lineages == None:
            firstline.append('taxonomy')
        firstline = '\t'.join(firstline) + '\n'
        File.write(firstline)
        lines = [[] for i in range(len(self.otuIDs))]
        for i in range(len(self.otuIDs)):
            lines[i].append(self.otuIDs[i])
            for col in self.matrix[i]:
                lines[i].append(str(col))
            if not self.lineages == None:
                lines[i].append(self.lineages[i])
            lines[i] = '\t'.join(lines[i])
        lines = '\n'.join(lines)
        File.write(lines)
    
    def dropSamplesNotInObj(self, mapFileObj):
        '''
        The mapFileObj (e.g. a QiimeMap object) must have an
        "object.sampleIDs" list. SampleIDs in the OTU table that are not
        in the map object will be removed from the OTU table. The map 
        object will remain unchanged.
        '''
        notInMap = [i for i in range(len(self.sampleIDs))
                    if not self.sampleIDs[i] in mapFileObj.getSampleIDs()]
        print("Dropping %i samples from the OTU table" % (len(notInMap)))
        if len(notInMap) == len(self.sampleIDs):
            raise Exception("The OTU table has no sample identifiers in common with the map file.\nThe OTU table was not altered")
        notInMap.sort()
        notInMap.reverse()
        for notIndex in notInMap:
            self.dropSamplebyIndex(notIndex)
    
    def taxonomyOfOTUs(self, otuIDs):
        '''
        Takes an otu as a string, or a list of one or more otus and
        returns the corresponding taxonomy associations as a list.
        If no taxonomy assignment is available, this will return (None)
        '''
        if isinstance(otuIDs, str):
            otuIDs = [otuIDs]
        if self.lineages == None:
            sys.stderr.write("There are no taxonomy assignments!\n")
            return None
        else:
            return [self.lineages[self.otuIDs.index(otuID)] for otuID in otuIDs]


class QiimeMap:
    def __init__(self, filename):
        '''
        A Qiime file begins with '#SampleID' in the first line.
        The first line is a load of headers.
        Every subsequent line is the relevant data under the header.
        Data may be omitted but a null string will still count as a
        category or datapoint or whatever.
        All lines must have the same number of columns as the header.
        All delimitation must be tabbed.
        '''
        self.matrix = [line.rstrip('\n').split('\t') for line in
            open(filename,'r').readlines() if not line.startswith('#')
            or line.startswith('#SampleID\t')]
        assert self.matrix[0][0] == '#SampleID', \
            "A valid Qiime map file must begin with '#SampleID'"
        for row in self.matrix[1:]:
            assert len(row) == len(self.matrix[0]), \
                "The row:\n%s\nhad a different number of columns to the header" % ('\t'.join(row))
        self.headers = self.matrix.pop(0)
        cut = self.headers.pop(0)
        self.sampleIDs = [self.matrix[i].pop(0) for i in range(len(self.matrix))]
        self.backup()
    
    def backup(self):
        '''Saves a copy of the current self.matrix and self.headers for
        convenient restoration. This method is called automatically upon
        instantiation and is user-callable thereafter (to save changes at
        an intermediate point without the possibility of reverting to
        the original state any more).
        '''
        self.__originalmatrix__ = copy.deepcopy(self.matrix)
        self.__originalheaders__ = copy.copy(self.headers)
        self.__originalsampleIDs__ = copy.copy(self.sampleIDs)
    
    def restore(self):
        self.matrix = copy.copy(self.__originalmatrix__)
        self.headers = copy.copy(self.__originalheaders__)
        self.sampleIDs = copy.copy(self.__originalsampleIDs__)
    
    def __len__(self):
        return(len(self.sampleIDs))
    
    def listColumns(self, indices=False):
        '''List the columns of the Map File. If any non-False value is
        given for indices, the columns will be printed in a numbered
        fashion (0-based counting), one column per line.
        Otherwise, the columns will be printed without indices on one
        tab-delimited line.
        '''
        if not indices:
            print('\n'.join(self.headers))
        else:
            for a in enumerate(self.headers):
                print("%i\t%s" % (a[0], a[1]))
    
    def getSampleIDs(self):
        return [sampleID for sampleID in self.sampleIDs]
    
    def getColumns(self):
        '''Returns a list of the column headers'''
        return self.headers
    
    def listCategories(self, column):
        '''List all unique values under a certain column. The column
        can be specified by an exact string or by an index (0-based)
        '''
        if isinstance(column, int):
            colindex = column
        elif isinstance(column, str):
            try:
                colindex = self.headers.index(column)
            except ValueError:
                print("Error: %s is not in the map file's headers" % (column))
        else:
            raise Exception("'column' must be an integer index or a string")
        categories = list(set([line[colindex] for line in self.matrix]))
        print("\n".join(categories))
    
    def getCategories(self, column):
        '''
        'column' may be an integer index of or a string match to a
        column in the map file's header.
        Returns a list of unique strings under the column.
        '''
        if isinstance(column, int):
            colindex = column
        elif isinstance(column, str):
            colindex = self.headers.index(column)
        else:
            raise Exception("'column' must be a string or an integer")
        categories = list(set([row[colindex] for row in self.matrix]))
        return categories
    
    def dropCategories(self, column, categories):
        '''Drops one or more rows from the mapfile object if they are 
        labeled with a category from 'categories' under the specified 
        column. The column can be an exact string match of a column 
        header or an integer index. The categories can be a single 
        category as a string, or a list of one or more categories (as 
        strings).
        '''
        if isinstance(column, int):
            colindex = column
        elif isinstance(column, str):
            colindex = self.headers.index(column)
        else:
            raise Exception("'column' must be a string or an integer")
        
        if isinstance(categories, str):
            cats = [categories]
        elif isinstance(categories, list):
            cats = categories
        else:
            raise Exception("'categories' must be a string or a list of strings")
        dropindices = [i for i in range(len(self.matrix)) if 
                       self.matrix[i][colindex] in cats]
        dropindices.reverse()
        for dropindex in dropindices:
            cut = self.matrix.pop(dropindex)
            cut = self.sampleIDs.pop(dropindex)
        del(dropindices) # Shouldn't this happen automatically once this namespace is out of context?
    
    def dropSamples(self, sampleIDs):
        '''
        'sampleIDs' can be a string (one sampleID) or a list of one or
        more sampleIDs (a list of strings)
        The specified sampleID(s) will be dropped from the map file.
        '''
        if isinstance(sampleIDs, str):
            sIDs = [sampleIDs]
        elif isinstance(sampleIDs, list):
            sIDs = sampleIDs
        else:
            raise Exception("'sampleIDs' must be a string or a list")
        dropindices = [i for i in range(len(self.sampleIDs)) if
                       self.sampleIDs[i] in sIDs]
        dropindices.reverse()
        for dropindex in dropindices:
            cut = self.matrix.pop(dropindex)
            cut = self.sampleIDs.pop(dropindex)
        del(dropindices)
    
    def dropSamplesNotInObj(self, OTUtabObj):
        '''
        This compares the sample IDs in the current MapFile instance 
        with using an 'Obj.getSampleIDs()' call. Any sampleIDs that are 
        in the MapFile, but not in the object will be removed from the 
        MapFile instance. The object instance provided in the argument 
        is not changed.
        '''
        droppedSamples = [self.sampleIDs[i] for i in range(len(self.sampleIDs))
                          if not self.sampleIDs[i] in OTUtabObj.getSampleIDs()]
        for droppedSample in droppedSamples:
            self.dropSamples(droppedSample)
            
    def getSampleIdsByCategory(self, column, category, inverse=False):
        '''
        'column' can be an integer index or a string
        'categories' must be a string
        Returns a list of the corresponding sampleIDs.
        If inverse is set to True, all sample IDs except those matching
        will be returned
        '''
        if isinstance(column, int):
            colindex = column
        elif isinstance(column, str):
            colindex = self.headers.index(column)
        else:
            raise Exception("'column' must be an integer or a string")
        catIndices = [i for i in range(len(self.matrix)) if
                      self.matrix[i][colindex] == category]
        if inverse == True:
            catIndices = [i for i in range(len(self.matrix))
                          if not i in catIndices]
        sampleIDs = [self.sampleIDs[i] for i in catIndices]
        return sampleIDs
    
    def getCategoryFromSampleID(self, column, sampleID):
        '''
        Given a sampleID and a header column from the mapfile, this will
        return the string value for the sampleID under the given column
        '''
        colIndex = self.headers.index(column)
        sampleIndex = self.sampleIDs.index(sampleID)
        return self.matrix[sampleIndex][colIndex]
    
    def getValueOfSample(self, sampleID, column, trytype=str, allowException='N/A'):
        cIndex = self.headers.index(column)
        try:
            sIndex = self.sampleIDs.index(sampleID)
            try:
                v = trytype(self.matrix[sIndex][cIndex])
            except ValueError:
                v = self.matrix[sIndex][cIndex]
            return v
        except ValueError:
            if sampleID == allowException or sampleID in allowException:
                return(sampleID)
            else:
                raise ValueError("ValueError: '%s' is not in the list" % sampleID)
    
    def getValuesOfSamples(self, sampleIDs, column, trytype=str, allowException='N/A'):
        """
        Input
        -----
        A sample identifier in string format OR a (nested) list of sample
        identifiers which will be traversed recursively.
        
        A column heading from which to retrieve values corresponding to
        sample identifiers.
        
        
        Returns
        -------
        Either a string or a (nested) list of values, corresponding to the
        input
        """
        if isinstance(sampleIDs, str):
            return self.getValueOfSample(sampleIDs, column, trytype=trytype,
                   allowException=allowException)
        elif isinstance(sampleIDs, list):
            return [self.getValuesOfSamples(elem, column, trytype=trytype,
                    allowException=allowException) for elem in sampleIDs]
        else:
            print("Error: Unexpected data type. Either a list or a string was expected")
            quit()
    
    
    def timeCourse(self, participantIdHeader, timeCourseHeader,
                   timeIsNumeric=False, timeCourseOrder=None,
                   participants=None, subset=None, dataNotFound='N/A'):
        """
        Each participant listed under participantIdHeader may have a
        sample at multiple time points.
        
        Each entry under timeCourseHeader will be sorted using python's
        built in sort function. These entries will be treated as strings
        unless they are explicitly stated as numeric using 
        timeIsNumeric=True. A specific order for the time course may be
        specified with a list given as a value for timeCourseOrder.
        
        The output can be filtered to include only a certain subset of
        the participants. This is done by passing a two element tuple to
        'subset'. This first element of the tuple must be a string
        exactly matching a header in the mapfile. The second element can
        be either a single string or a list of strings exactly matching
        values under the header in the mapfile.
        
        In cases where a sample is missing or it is excluded because it
        had too few reads etc. the value returned for this sample will
        be the value specified for 'dataNotFound'.
        
        Returns:
        --------
        A three element list:
        1. A two dimensonal list of sample Identifiers. The number of
           lists corresponds to the number of participants. The length
           of each list corresponds to the number of time-points.
        2. A list of participant identifiers corresponding to the order
           in 1.
        3. A list of the timepoints corresponding to the order in 1.
        """
        # It would be good to change this to ignore rows where the 
        # timeCourseHeader column are ignored when they can't be
        # converted into floats
        pIndex = self.headers.index(participantIdHeader)
        pIDs = list(set([self.matrix[i][pIndex] for i in range(len(self.matrix))]))
        pIDs.sort()
        if isinstance(participants,list):
            pIDs = participants
        if isinstance(subset,tuple):
            cIndex = self.headers.index(subset[0])
            if isinstance(subset[1],str):
                js = [subset[1]]
            else:
                js = subset[1]
            for ji in js:
                if not ji in self.getCategories(subset[0]):
                    warn_message = "Warning: " + ji + " was not found under " + subset[0] + "\n"
                    warnings.warn(warn_message)
            ks = list(set([self.matrix[j][pIndex] for j in range(len(self.matrix)) if self.matrix[j][cIndex] in js]))
            pIDs = [pID for pID in pIDs if pID in ks]
        tIndex = self.headers.index(timeCourseHeader)
        tps = list(set([self.matrix[i][tIndex] for i in range(len(self.matrix))]))
        tps.sort()
        if timeIsNumeric == True:
            tpz = [float(tp) for tp in tps]
            tps = [tp for _,tp in sorted(zip(tpz,tps))]
        if isinstance(timeCourseOrder,str):
            timeCourseOrder = [timeCourseOrder]
        if isinstance(timeCourseOrder,list):
            tps = [str(x) for x in timeCourseOrder]
        tc = [[dataNotFound for tp in tps] for pID in pIDs]
        for i in range(len(pIDs)):
            for j in range(len(tps)):
                for k in range(len(self.matrix)):
                    if self.matrix[k][pIndex] == pIDs[i] and self.matrix[k][tIndex] == tps[j]:
                        tc[i][j] = self.sampleIDs[k]
        return [tc,pIDs,tps]
        
    
    def writeQiimeMap(self,File):
        '''
        Write a qiime v1.x compatible map file to File. If File already
        exists it will be OVERWRITTEN!!!
        '''
        fh = open(File,'w')
        output = '#SampleID\t' + '\t'.join(self.headers) + '\n'
        fh.write(output)
        outputs = []
        for i in range(len(self.sampleIDs)):
            output = self.sampleIDs[i] + '\t' + '\t'.join([str(col) for col in self.matrix[i]]) # This list comprehension will support modifying some columns of the metadata to ints or floats etc.
            outputs.append(output)
        output = '\n'.join(outputs)
        fh.write(output)
        fh.close()


MAPFILE = 'mapfile.txt' # A qiime1.X format mapfile is compatible
# The mapfile used will be available as a supplementary file in a
# peer-reviewed publication
OTUTAB = 'zotutab.tax.tab' # A tab-delimited feature table
# The feature table used here will be available as a supplementary file
# in a peer-reviewed publication

def dropNA(listObj):
    return [elem for elem in listObj if not elem == 'N/A']

def countPositive(floatList):
    count = 0
    for i in range(len(floatList)):
        if floatList[i] > 0:
            count += 1
    return count

def ceiling(num, ceil):
    '''Applies an upper limit to how big a number can be
    
    If `num' is greater than `ceil' the return will be `ceil'.
    Otherwise the return will be `num'.
    '''
    if num > ceil:
        num = ceil
    return num

def cleanTC(tcObj):
    '''Removes subjects without enough samples
    
    tcObj is the return from a call of
    AmpliconDataStructures.QiimeMap.timeCourse
    
    If the first time-point is missing (equal to 'N/A') or any two
    time-points are missing for a single subject, that subject will be
    removed.
    
    There is no return. The timeCourse structure is modified in place
    '''
    dropIndices = []
    for i in range(len(tcObj[0])):
        if tcObj[0][i][0] == 'N/A' or tcObj[0][i].count('N/A') >= 2:
            dropIndices.append(i)
    dropIndices.sort(reverse=True)
    for dI in dropIndices:
        _ = tcObj[0].pop(dI)
        _ = tcObj[1].pop(dI)

def presence_sc(sid1, sid2, otuID, otutab): #, sigC=6):
    '''Calculates a consistency score between two consecutive samples
    
    sid1:
    
    sid2:
    
    otuID:
    
    otutab:
    
    sigC:
    '''
    total1 = int(otutab.getSampleCoverage(sid1))
    total2 = int(otutab.getSampleCoverage(sid2))
    abnd1 = int(otutab.getOTUabundance(otuID, sid1))
    abnd2 = int(otutab.getOTUabundance(otuID, sid2))
    exp1 = total1 * (float(abnd2) / total2)
    exp2 = total2 * (float(abnd1) / total1)
    if abnd1 == 0 and abnd2 == 0:
        prs = float(0)
    elif abnd1 > 0 and abnd2 == 0:
        prs = 1 - poisson.sf(0, exp2)
    elif abnd1 == 0 and abnd2 > 0:
        prs = 1 - poisson.sf(0, exp1)
    elif abnd1 > 0 and abnd2 > 0:
        prs = 1
    return prs

def loss_sc(obs, exp):
    '''Returns a loss score up to 1, adjusted based on an expected value
    
    obs: osbserved count
    
    exp: expected count (float)
    '''
    if obs > 0:
        lsc = 0
    else:
        lsc = poisson.sf(0, exp)
    return lsc

def add_block(listObj, decimalPlaces=2, blockSize=5):
    '''Adds a list of strings, integers or floats to a list for printing
    
    listObj: A list of integers, floats or strings. Cannot be longer
             than blockSize unless that's set to None
    
    decimalPlaces: An integer or None. Floats will be rounded off to
                   this precision (or printed to python's level of
                   precision if it's None)
    
    blockSize: This is used to keep lines lined up when each line may
               have lists of differing lengths.
    
    Returns a list of strings which can be printed with '\t'.join(). The
    first element of the list will always be an empty string (ie. '')
    because this is intended to add a block of related items to another
    block of related items and the empty element provides a column of
    padding. The returned list will have a length of blockSize + 1
    (or len(listObj) + 1).
    '''
    if blockSize == None:
        blockSize = len(listObj)
    assert len(listObj) <= blockSize
    outList = [elem for elem in listObj]
    for i in range(len(outList)):
        if isinstance(outList[i], int):
            outList[i] = str(outList[i])
        if isinstance(outList[i], float):
            if not decimalPlaces == None:
                prcsn = '{0:.' + str(decimalPlaces) + 'f}'
                outList[i] = prcsn.format(outList[i])
            else:
                outList[i] = str(outList[i])
    while len(outList) < blockSize:
        outList.append('')
    return [''] + outList
               

def glscore_V3(sids, otuID, otutab):
    '''Write something informative here'''
    sids = dropNA(sids)
    coverages = [int(otutab.getSampleCoverage(sid)) for sid in sids]
    abnds = [int(otutab.getOTUabundance(otuID, sid)) for sid in sids]
    t0expAbnd = average([(float(abnds[i]) / coverages[i]) * coverages[0] for
                  i in range(1, len(sids))]) # `max' instead of `average'?
    tTexpAbnds = [(float(abnds[0]) / coverages[0]) * coverages[i] for
                  i in range(1, len(sids))]
    if abnds[0] == 0:
        lsc = 0
        # If we don't expect to see much, if any, of the OTU at t0 due
        # to low sample coverage, then we should penalise the gain score
        penalty = poisson.sf(0, t0expAbnd) # ceiling(t0expAbnd, sigC) / sigC
        # If the OTU isn't consistently there, then it could be
        # an OTU that fluctuates naturally and wasn't gained
        prss = [presence_sc(sids[i], sids[i+1], otuID, otutab)
                for i in range(1, len(sids) - 1)]
        prss = [prss[i] * prss[i+1] for i in range(len(prss) - 1)]
        prs = average(prss)
        gsc = prs * penalty
    if abnds[0] > 0:
        gsc = 0
        # I want to penalise the loss scores of anything that could
        # have been lost because of uneven sample coverage
        lscs = [loss_sc(abnds[i+1], tTexpAbnds[i])
                for i in range(len(tTexpAbnds))]
        lscs = [lscs[i] * lscs[i+1] for i in range(len(lscs) - 1)]
        lsc = average(lscs)
    output = add_block(sids) + add_block(coverages) + add_block(abnds) \
             + add_block([t0expAbnd] + tTexpAbnds)
    output.append('')
    output.append('{0:.5f}'.format(gsc))
    output.append('{0:.5f}'.format(lsc))
    return output
    

if __name__ == '__main__':
    otutab = otuTable(OTUTAB, taxonomy='taxonomy')
    qmap = QiimeMap(MAPFILE)
    donor_sids = qmap.getSampleIdsByCategory('Treatment', 'Donor')
    qmap.dropCategories('Treatment', 'Donor')
    
    tcALL = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True)
    tcFMT = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                            subset=('Treatment', 'FMT'))
    tcPCB = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                            subset=('Treatment', 'Placebo'))
    tcFmtRes = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                               subset=('GroupOutcome', 'FmtResponder'))
    tcFmtNon = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                               subset=('GroupOutcome', 'FmtNonResponder'))
    tcPcbRes = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                               subset=('GroupOutcome', 'PlaceboResponder'))
    tcPcbNon = qmap.timeCourse('PatientNumber', 'TimepointInMonths', timeIsNumeric=True,
                               subset=('GroupOutcome', 'PlaceboNonResponder'))
    cleanTC(tcALL)
    cleanTC(tcFMT)
    cleanTC(tcPCB)
    cleanTC(tcFmtRes)
    cleanTC(tcFmtNon)
    cleanTC(tcPcbRes)
    cleanTC(tcPcbNon)
    
    # Comment out blocks below as appropriate to see different groupings
    # of the output.
    for otuID in otutab.getOTUlist():
        for i in range(len(tcALL[1])):
            pid = tcALL[1][i]
            sids = tcALL[0][i]
            gl = glscore_V3(sids, otuID, otutab)
            output = [pid, otuID] + gl
            print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcFMT[1])):
            # pid = tcFMT[1][i]
            # sids = tcFMT[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcPCB[1])):
            # pid = tcPCB[1][i]
            # sids = tcPCB[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcFmtRes[1])):
            # pid = tcFmtRes[1][i]
            # sids = tcFmtRes[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcFmtNon[1])):
            # pid = tcFmtNon[1][i]
            # sids = tcFmtNon[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcPcbRes[1])):
            # pid = tcPcbRes[1][i]
            # sids = tcPcbRes[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))
    
    # for otuID in otutab.getOTUlist():
        # for i in range(len(tcPcbNon[1])):
            # pid = tcPcbNon[1][i]
            # sids = tcPcbNon[0][i]
            # gl = glscore_V3(sids, otuID, otutab)
            # output = [pid, otuID] + gl
            # print('\t'.join(output))

