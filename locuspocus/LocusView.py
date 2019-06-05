
from functools import wraps
from .Locus import Locus,LocusAttrs,SubLoci
from .Exceptions import *

__all__ = ['LocusView']


class AttrsView(LocusAttrs):
    def __init__(self,LID, refloci):
        self._LID = LID
        self._ref = refloci

    def keys(self):
        cur = self._ref._db.cursor()
        results = cur.execute('''
            SELECT key FROM loci_attrs 
            WHERE LID = ?
            ''',
            (self._LID,)
        )
        return [k[0] for k in results]

    def values(self):
        cur = self._ref._db.cursor()
        results = cur.execute('''
            SELECT val FROM loci_attrs 
            WHERE LID = ?
            ''',
            (self._LID,)
        ) 
        return [k[0] for k in results]

    def items(self):
        return zip(self.keys(),self.values())

    def __getitem__(self,key):
        cur = self._ref._db.cursor()
        try:
            val, = cur.execute('''
                SELECT val FROM loci_attrs
                WHERE LID = ? AND key = ?
            ''',(self._LID,key)).fetchone()
        except TypeError:
            raise KeyError(f'"{key}" in in attrs')
        return val 

    def __setitem__(self,key,val):
        cur = self._ref._db.cursor()
        cur.execute('''
            INSERT OR REPLACE INTO loci_attrs
            (LID,key,val)
            VALUES (?,?,?)
        ''',(self._LID,key,val))

    def __repr__(self):
        return '{'+','.join([':'.join([x,y]) for x,y in self.items()]) +'}'


class SubLociView(SubLoci):
    # A restricted list interface to subloci
    def __init__(self, LID, refloci):
        self._LID = LID
        self._ref = refloci
        # Grab Views for each sublocus
        self.LIDs = [x[0] for x in \
            self._ref._db.cursor().execute('''
            SELECT child FROM relationships
            WHERE parent = ?
        ''',(self._LID,)).fetchall()]

    def __iter__(self):
        return (LocusView(x,self._ref) for x in self.LIDs)

    def add(self,locus):
        pass

    def __getitem__(self,index):
        # return self.loci[index]
        return LocusView(self.LIDs[index],self._ref)

    def __len__(self):
        # return len(self.loci)
        return len(self.LIDs)

    def __repr__(self):
        return '['+','.join([repr(x) for x in self])+']'

class LocusView(Locus):
    '''
    A LocusView is an efficient way to access
    Locus objects 
    '''
    def __init__(self, LID, refloci):
        self._LID = LID
        self._ref = refloci

        self.attrs = AttrsView(self._LID,self._ref)
        self.subloci = SubLociView(self._LID,self._ref) 

    def _property(self,name):
        val, = self._ref._db.cursor().execute(
            f"SELECT {name} FROM loci WHERE LID = ?",
            (self._LID,)
        ).fetchone()
        return val

    @property
    def chromosome(self):
        return self._property('chromosome')
        
    @property
    def start(self):
        return self._property('start')

    @property
    def end(self):
        return self._property('end')

    @property
    def feature_type(self):
        return self._property('feature_type')

    @property
    def strand(self):
        return self._property('strand')

    @property
    def frame(self):
        return self._property('frame')

    @property
    def name(self):
        return self._property('name')

    def add_sublocus(self,locus):
        self.subloci.add(locus)

   #def __repr__(self):
   #    return (
   #        f'LocusView('
   #        f'chromosome={self.chromosome},'
   #        f'start={self.start},'
   #        f'end={self.end},'
   #        f'feature_type={self.feature_type},{self.strand},{self.frame},{self.name})'
