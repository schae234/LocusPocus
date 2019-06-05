import pytest
import numpy as np
import minus80 as m80

from itertools import chain
from locuspocus import Locus,RefLoci

from locuspocus.Exceptions import StrandError,ChromosomeError

@pytest.fixture(scope='module')
def SimpleRefLoci():
    a = Locus('1',10,20)
    b = Locus('1',20,30)
    c = Locus('2',30,40)
    d = Locus('2',40,50)

    x = Locus(
        '1', 100, 200,
        attrs={'foo':'bar'},
        name='x',
        subloci=[a,b]
    )
    y = Locus(
        '2', 100, 200,
        attrs={'foo':'bar'},
        name='y',
        subloci=[c,d]
    )
    if m80.Tools.available('RefLoci','test'):
        m80.Tools.delete('RefLoci','test',force=True)
    ref = RefLoci('test')   
    ref.add_locus(x)
    ref.add_locus(y)
    return ref

@pytest.fixture
def simpleLocusView(SimpleRefLoci):
    return SimpleRefLoci['x']

def test_initialization(simpleLocusView):
    # numeric chromosomes
    assert simpleLocusView.chromosome == '1'
    assert simpleLocusView.start == 100
    assert simpleLocusView.end == 200
    assert len(simpleLocusView) == 101

def test_getitem(simpleLocusView):
    assert simpleLocusView['foo'] == 'bar'

def test_default_getitem(simpleLocusView):
    assert simpleLocusView.default_getitem('name', 'default') == 'default'

def test_start(simpleLocusView):
    assert simpleLocusView.start == 100

def test_plus_stranded_start():
    l = Locus('1',1,100,strand='+')
    assert l.stranded_start == 1 

def test_minus_stranded_start():
    l = Locus('1',1,100,strand='-')
    assert l.stranded_start == 100 

def test_end(simpleLocusView):
    assert simpleLocusView.end == 200

def test_plus_stranded_end():
    l = Locus('1',1,100,strand='+')
    assert l.stranded_end == 100 

def test_minus_stranded_end():
    l = Locus('1',1,100,strand='-')
    assert l.stranded_end == 1 

def test_coor(simpleLocusView):
    assert simpleLocusView.coor == (100, 200)

def test_upstream(simpleLocusView):
    assert simpleLocusView.upstream(50) == 50

def test_upstream_minus_strand(simpleLocusView):
    l = Locus('1',1,100,strand='-')
    assert l.upstream(50) == 150

def test_downstream(simpleLocusView):
    assert simpleLocusView.downstream(50) == 250

def test_downstream_minus_strand(simpleLocusView):
    l = Locus('1',100,200,strand='-')
    assert l.downstream(50) == 50

def test_center():
    l = Locus('1',100,200,strand='-')
    assert l.center == 150.5 

def test_name(simpleLocusView):
    assert simpleLocusView.name == 'x'

def test_eq(simpleLocusView):
    another_Locus = Locus(1, 110, 220)
    assert simpleLocusView == simpleLocusView
    assert simpleLocusView != another_Locus

def test_loci_lt_by_chrom():
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert x < y

def test_loci_gt_by_chrom():
    x = Locus('1',1,1)
    y = Locus('2',1,1)
    assert y > x

def test_loci_lt_by_pos():
    x = Locus('1',1,100)
    y = Locus('1',2,100)
    assert x < y

def test_loci_gt_by_pos():
    x = Locus('1',1,100)
    y = Locus('1',2,200)
    assert y > x

def test_len(simpleLocusView):
    assert len(simpleLocusView) == 101
    assert len(Locus(1, 100, 100)) == 1

def test_lt(simpleLocusView):
    same_chrom_Locus = Locus('1', 110, 220)
    diff_chrom_Locus = Locus('2', 90, 150)
    assert simpleLocusView < same_chrom_Locus
    assert simpleLocusView < diff_chrom_Locus

def test_gt(simpleLocusView):
    same_chrom_Locus = Locus('1', 90, 150)
    diff_chrom_Locus = Locus('2', 90, 150)
    assert simpleLocusView > same_chrom_Locus
    assert diff_chrom_Locus > simpleLocusView

def test_repr(simpleLocusView):
    assert repr(simpleLocusView) 

def test_hash(simpleLocusView):
    assert hash(simpleLocusView) == 108424482858392584

def test_subloci_getitem():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus')
    x.add_sublocus(y)
    assert x.subloci[0].name == 'sublocus'

def test_subloci_iter():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.add_sublocus(y)
    x.add_sublocus(z)
    for sub in x.subloci:
        assert sub.chromosome == '1' 

def test_subloci_len():
    x = Locus('1',1,2)
    y = Locus('1',3,4,name='sublocus1')
    z = Locus('1',3,4,name='sublocus2')

    x.add_sublocus(y)
    x.add_sublocus(z)
    assert len(x.subloci) == 2

def test_attrs_keys_method():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert sorted(x.attrs.keys()) == ['bar','foo']

def test_attrs_keys_method_empty():
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.keys())) == 0

def test_attrs_vals_method():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert len(sorted(x.attrs.values())) == 2

def test_attrs_vals_method_empty():
    x = Locus('1',3,4,attrs={})
    assert len(list(x.attrs.values())) == 0

def test_attrs_items():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert len(sorted(x.attrs.items())) == 2

def test_attrs_contains():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert 'foo' in x.attrs

def test_attrs_repr():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    assert repr(x.attrs)

def test_le_equals():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x <= y

def test_le_less():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',30,40)
    assert x <= y

def test_ge_equals():
    x = Locus('1',3,4,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x >= y

def test_ge_greater():
    x = Locus('1',30,40,attrs={'foo':'locus1','bar':'baz'})
    y = Locus('1',3,4)
    assert x >= y

def test_stranded_start_invalid():
    # Strand cannot be '='
    x = Locus('1',3,4,strand='=')
    with pytest.raises(StrandError):
        x.stranded_start

def test_stranded_stop_invalid():
    # Strand cannot be '='
    x = Locus('1',3,4,strand='=')
    with pytest.raises(StrandError):
        x.stranded_end

def test_as_record():
    x = Locus('1',3,4,strand='+')
    # This doesn't compare the dictionaries of each ...
    assert x.as_record()[0] == ('1', 3, 4, 'locuspocus', 'locus', '+', None, None, 1030788771219373542)

def test_center_distance():
    x = Locus('1',1,100,strand='+')
    # This needs to be 201 since x starts at 1
    y = Locus('1',201,300,strand='=')
    assert x.center_distance(y) == 200

def test_center_distance_different_chroms():
    x = Locus('1',1,100,strand='+')
    # This needs to be 201 since x starts at 1
    y = Locus('2',201,300,strand='+')
    assert x.center_distance(y) == np.inf

def test_str():
    x = Locus('1',1,100,strand='+')
    assert str(x) == repr(x)

def test_combine():
    x = Locus('1',1,2)
    y = Locus('1',3,4)
    z = x.combine(y)

    assert z.start == 1
    assert z.end == 4
    assert x in z.subloci
    assert y in z.subloci

def test_combine_chromosome_mismatch():
    x = Locus('1',1,2)
    y = Locus('2',3,4)
    with pytest.raises(ChromosomeError):
        z = x.combine(y)

def test_distance():
    x = Locus('1',1,100)
    y = Locus('1',150,250)
    assert x.distance(y) == 49

def test_distance_diff_chroms():
    x = Locus('1',1,100)
    y = Locus('2',150,250)

def test_get_subloic_by_index(SimpleRefLoci):
    x = SimpleRefLoci['x']
    assert x.subloci[0]

def test_set_attr(SimpleRefLoci):
    x = SimpleRefLoci['x']
    assert x['foo'] == 'bar'
    x['foo'] = 'baz'
    assert x['foo'] == 'baz'
    # change it back
    x['foo'] = 'bar'
    # check it
    assert x['foo'] == 'bar'


def test_add_sublocus(SimpleRefLoci):
    x = SimpleRefLoci['x']
