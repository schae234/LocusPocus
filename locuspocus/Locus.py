#!/usr/bin/python3
from collections import defaultdict
from dataclasses import dataclass,field,InitVar,FrozenInstanceError
from itertools import chain
from typing import Union, Any, List, Optional, cast, Callable, Iterable, ClassVar, Dict

from .Exceptions import StrandError, ChromosomeError

import re
import math
import hashlib
import locuspocus
import dataclasses
import uuid


import pandas as pd
import numpy as np

__all__ = ['Locus']

class SubLoci():
    # A restricted list interface to subloci
    def __init__(self,loci=None):
        if loci is None:
            loci = list()
        self._loci = loci

    def __iter__(self):
        return (x for x in self._loci)

    def add(self,locus):
        self._loci.append(locus)

    def __getitem__(self,index):
        return self._loci[index]

    def __len__(self):
        return len(self._loci)
    

class LocusAttrs():
    # a restricted dict interface to attributes
    def __init__(self,attrs=None):
        if attrs is None:
            attrs = dict()
        self._attrs = attrs

    def keys(self):
        return self._attrs.keys()

    def values(self):
        return self._attrs.values()

    def items(self):
        return self._attrs.items()

    def __getitem__(self,key):
        return self._attrs[key]

    def __setitem__(self,key,val):
        self._attrs[key] = val

    def __repr__(self):
        return repr(self._attrs)

@dataclass()
class Locus:
    chromosome: str
    start: int 
    end: int

    source: str = 'locuspocus'
    feature_type: str = 'locus'
    strand: str = '+'
    frame: int = None
    name: str = None

    # 
    uuid: uuid = field(init=False,repr=False)

    # Extra locus stuff
    #attrs: LocusAttrs = field(default_factory=LocusAttrs)
    attrs: ClassVar[Dict[str,str]] = {}
    # DAG info
    Tree: ClassVar[Dict[int,int]] = pd.DataFrame(columns=['parent','child'])
    parent: None = field(default=None,repr=False)

    def __post_init__(self):
        self.uuid = uuid.uuid4() 

    #-----------------------------------------
    #       Methods
    #-----------------------------------------


    def __len__(self):
        return abs(self.end - self.start) + 1

    def __eq__(self,locus):
        try:
            assert self.chromosome == locus.chromosome 
            assert self.start == locus.start 
            assert self.end == locus.end 
            assert self.source == locus.source 
            assert self.feature_type == locus.feature_type 
            assert self.strand == locus.strand 
            assert self.frame == locus.frame 
            assert self.name == locus.name 
            for key,val in self.attrs.items():
                assert val == locus.attrs[key]
            for x,y in zip(sorted(self.subloci),sorted(locus.subloci)):
                assert x == y
            return True
        except (AssertionError,KeyError):
            return False

    def __lt__(self,locus):
        if self.chromosome == locus.chromosome:
            return self.start < locus.start
        else:
            return self.chromosome < locus.chromosome

    def __le__(self,locus):
        if (self.chromosome,self.coor) == (locus.chromosome,locus.coor):
            return True
        else:
            return self < locus

    def __ge__(self,locus):
        if (self.chromosome,self.coor) == (locus.chromosome,locus.coor):
            return True
        else:
            return self > locus

    def __gt__(self,locus):
        if self.chromosome == locus.chromosome:
            return self.start > locus.start
        else:
            return self.chromosome > locus.chromosome

    def __hash__(self):
        '''
            Convert the locus to a hash, uses md5. The hash
            is computed using the primary

            Parameters
            ----------
            None

            Returns
            -------
            int : md5 hash of locus

        '''
        field_list = [str(x) for x in (
            self.chromosome,
            self.start,
            self.end,
            self.source,
            self.feature_type,
            self.strand,
            self.frame,
            self.name
        )]
        subloci_list = [str(hash(x)) for x in self.subloci]
        #attr_list = list(chain(*zip(map(str,self.attrs.keys()),map(str,self.attrs.values()))))
        # Create a full string
        loc_string  = "_".join(field_list + subloci_list) # attr_list + subloci_list)
        digest = hashlib.md5(str.encode(loc_string)).hexdigest()
        return int(digest, base=16)

    @property
    def stranded_start(self):
        if self.strand == '+':
            return min(self.coor)
        elif self.strand == '-':
            return max(self.coor)
        else:
            raise StrandError

    @property
    def stranded_end(self):
        if self.strand == '+':
            return max(self.coor)
        elif self.strand == '-':
            return min(self.coor)
        else:
            raise StrandError

    def __getitem__(self,item):
        return self.attrs[item]

    def __setitem__(self,key,val):
        self.attrs[key] = val

    def add_sublocus(self,locus):
        self.subloci.add(locus)

    def as_record(self):
        return ((
            self.chromosome,
            self.start,
            self.end,
            self.source,
            self.feature_type,
            self.strand,
            self.frame,
            self.name,
            hash(self)
        ),
            self.attrs
        )

    def default_getitem(self,key,default=None) -> Any:
        '''
        Returns the attr value of the Locus based on the key.
        similar to: Locus['Name'] -> "Locus1".
        If the attr key (e.g. 'Name') does not exist,
        a default can be specified.

        Parameters
        ----------
        key : str
            The key value of the attr
        default : Any
            If the key is not in the Locus attributes, return
            this value by default.
        '''
        try:
            val = self.attrs[key]
        except KeyError as e:
            val = default
        finally:
            return val

    @property
    def coor(self):
        '''
        Returns a tuple containing the start and end
        positions of the locus
        '''
        return (self.start,self.end)

    def upstream(self,distance: int) -> int:
        '''
        Calculates a base pair position 5' of the 
        locus.

        Parameters
        ----------
        distance : int
            The distance upstream of the locus
        '''
        if self.strand == '+':
            return max(0,self.start - distance)
        elif self.strand == '-':
            return self.end + distance

    def downstream(self,distance: int) -> int:
        '''
        Calculates a base pair position 3' of the 
        locus
        
        Parameters
        ----------
        distance : int
            The distance downstream of the locus
        '''
        if self.strand == '+':
            return self.end + distance
        elif self.strand == '-':
            return max(0,self.start) - distance
        #return self.end + distance

    @property
    def center(self):
        '''
        Calculates the center base pair position of
        the locus.

        NOTE: If the locus has an odd length, a 'half-bp'
              will be returned. 
              E.g: Locus('1',100,200).center == 

        Returns
        -------
        The center position
        '''
        return self.start + len(self)/2

    def distance(self,locus):
        '''
        Return the number of base pairs between two loci.
        NOTE: this excludes the start/end bases of the loci.

        Locus A             Locus B
        ==========---------=========
        1       10         20      30

        There are 9 bases between 10 and 20 (excluding positions
        10 and 20 themselves because they are included in the 
        loci).

        If the loci are on different chromosomes, return np.inf

        Parameters
        ----------
        locus: Locus Object
            A second locus object to calculate distance.

        Returns
        -------
        int: the number of bp between the loci
        np.inf: if on different chromosomes

        '''
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            x,y = sorted([self,locus])
            distance = y.start - x.end - 1
        return distance

    def center_distance(self, locus):
        '''
        Return the distance between the center of two loci. 
        If the loci are on different chromosomes, return np.inf.

        Parameters
        ----------
        locus : Locus Object
            A second locus object used to calculate a distance.

        Returns
        -------
        int : the distance between the center of two loci.
        np.inf: if on different chromosomes

        '''
        if self.chromosome != locus.chromosome:
            distance = np.inf
        else:
            distance = math.floor(abs(self.center - locus.center))
        return distance

    def combine(self, locus):
        '''
        Returns a new Locus with start and stop boundaries
        that contain both of the input loci. Both input loci are
        added to the subloci of the new Locus.

        NOTE: this ignores strand, the resultant Locus is just
              a container for the input loci.

        ___________Ascii Example__________________________
             
              Locus A          Locus B
        ------=========-------=============---------------
        
        A.combine(B)
        ------=============================---------------
              subloci=[A,B]
        '''
        if self.chromosome != locus.chromosome:
            raise ChromosomeError('Input Chromosomes do not match') 
        x,y = sorted([self,locus])
        start = x.start
        end = y.end
        return Locus(self.chromosome,start,end,subloci=[self,locus])

    def as_tree(self,parent=None): #pragma: no cover
        from anytree import Node, RenderTree
        root = Node(f'{self.feature_type}:{self.name}',parent=parent)
        for c in self.subloci:
            node = c.as_tree(parent=root)
        if parent is None:
            for pre, _, node in RenderTree(root):
                print("%s%s" % (pre, node.name))
        return root

    def __str__(self):
        return repr(self)
