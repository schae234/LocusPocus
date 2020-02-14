import pytest
import os

from locuspocus import Locus
from locuspocus import Loci
from locuspocus import Fasta

import minus80.Tools as m80tools

from locuspocus import Chromosome



@pytest.fixture(scope='module')
def realLocus():
    lines = [
        '9	ensembl	gene	831398	834611	.	+	.	ID=GRMZM2G158729;Name=GRMZM2G158729;biotype=protein_coding',
        '9	ensembl	mRNA	831398	834611	.	+	.	ID=GRMZM2G158729_T01;Parent=GRMZM2G158729;Name=GRMZM2G158729_T01;biotype=protein_coding',
        '9	ensembl	intron	831529	832121	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.132',
        '9	ensembl	intron	832202	832292	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.133',
        '9	ensembl	intron	832369	832441	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.134',
        '9	ensembl	intron	832653	832884	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.135',
        '9	ensembl	intron	833079	833156	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.136',
        '9	ensembl	intron	833263	833347	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.137',
        '9	ensembl	intron	833438	833654	.	+	.	Parent=GRMZM2G158729_T01;Name=intron.138',
        '9	ensembl	exon	831398	831528	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E02',
        '9	ensembl	exon	832122	832201	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E04',
        '9	ensembl	exon	832293	832368	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E08',
        '9	ensembl	exon	832442	832652	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E05',
        '9	ensembl	exon	832885	833078	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E11',
        '9	ensembl	exon	833157	833262	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E09',
        '9	ensembl	exon	833348	833437	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E03',
        '9	ensembl	exon	833655	834611	.	+	.	Parent=GRMZM2G158729_T01;Name=GRMZM2G158729_E07',
        '9	ensembl	CDS	831493	831528	.	+	.	Parent=GRMZM2G158729_T01;Name=CDS.147',
        '9	ensembl	CDS	832122	832201	.	+	0	Parent=GRMZM2G158729_T01;Name=CDS.148',
        '9	ensembl	CDS	832293	832368	.	+	2	Parent=GRMZM2G158729_T01;Name=CDS.149',
        '9	ensembl	CDS	832442	832652	.	+	0	Parent=GRMZM2G158729_T01;Name=CDS.150',
        '9	ensembl	CDS	832885	833078	.	+	1	Parent=GRMZM2G158729_T01;Name=CDS.151',
        '9	ensembl	CDS	833157	833262	.	+	0	Parent=GRMZM2G158729_T01;Name=CDS.152',
        '9	ensembl	CDS	833348	833437	.	+	1	Parent=GRMZM2G158729_T01;Name=CDS.153',
        '9	ensembl	CDS	833655	834394	.	+	1	Parent=GRMZM2G158729_T01;Name=CDS.154',
        '9	ensembl	mRNA	831418	834611	.	+	.	ID=GRMZM2G158729_T02;Parent=GRMZM2G158729;Name=GRMZM2G158729_T02;biotype=protein_coding',
        '9	ensembl	intron	831529	832121	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.155',
        '9	ensembl	intron	832202	832292	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.156',
        '9	ensembl	intron	832369	832441	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.157',
        '9	ensembl	intron	832653	832884	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.158',
        '9	ensembl	intron	833079	833156	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.159',
        '9	ensembl	intron	833438	833654	.	+	.	Parent=GRMZM2G158729_T02;Name=intron.160',
        '9	ensembl	exon	831418	831528	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E06',
        '9	ensembl	exon	832122	832201	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E04',
        '9	ensembl	exon	832293	832368	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E08',
        '9	ensembl	exon	832442	832652	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E05',
        '9	ensembl	exon	832885	833078	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E11',
        '9	ensembl	exon	833157	833437	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E01',
        '9	ensembl	exon	833655	834611	.	+	.	Parent=GRMZM2G158729_T02;Name=GRMZM2G158729_E10',
        '9	ensembl	CDS	831493	831528	.	+	.	Parent=GRMZM2G158729_T02;Name=CDS.168',
        '9	ensembl	CDS	832122	832201	.	+	0	Parent=GRMZM2G158729_T02;Name=CDS.169',
        '9	ensembl	CDS	832293	832368	.	+	2	Parent=GRMZM2G158729_T02;Name=CDS.170',
        '9	ensembl	CDS	832442	832652	.	+	0	Parent=GRMZM2G158729_T02;Name=CDS.171',
        '9	ensembl	CDS	832885	833078	.	+	1	Parent=GRMZM2G158729_T02;Name=CDS.172',
        '9	ensembl	CDS	833157	833312	.	+	0	Parent=GRMZM2G158729_T02;Name=CDS.173'
    ]

    # Create the top level locus
    locus = Locus.from_gff_line(lines[0])
    # add all subloci
    for l in lines[1:]:
        locus.add_sublocus(Locus.from_gff_line(l))
    # Return it
    return locus

@pytest.fixture(scope='module')
def simpleLoci():
    m80tools.delete('Loci','simpleLoci')
    # Create a Locus
    a = Locus(1,100,150, id='gene_a')
    # Create a couple more!
    b = Locus(1,160,175, id='gene_b')
    c = Locus(1,180,200, id='gene_c')
    d = Locus(1,210,300, id='gene_d')
    e = Locus(2,100,150, id='gene_e')

    x = Loci('simpleLoci')
    x.add_loci([a,b,c,d,e])
    return x

@pytest.fixture(scope="module")
def testRefGen():
    # We have to build it
    if m80tools.available('Loci','Zm5bFGS'):
        x =  Loci('Zm5bFGS')
    else:
        m80tools.delete('Loci','Zm5bFGS')
        gff = os.path.expanduser(
            os.path.join(
                'raw', 
                'ZmB73_5b_FGS.gff.gz'
            )
        )
        x = Loci('Zm5bFGS')
        if len(x) == 0:
            x.import_gff(  
                gff
            )
    x.set_primary_feature_type('gene')
    return x


@pytest.fixture(scope='module')
def m80_Fasta():
    '''
        Create a Fasta which doesn't get 
        returned. Access the Fasta through
        the m80 API
    '''
    # delete the onl
    m80tools.delete('Fasta','ACGT')
    f = Fasta.from_file('ACGT','raw/ACGT.fasta')
    return True


@pytest.fixture(scope='module')                                                               
def smpl_fasta():                                                                  
    ''' A simple fasta that agrees with smpl_annot'''                           
    m80tools.delete('Fasta','smpl_fasta')
    fasta = Fasta('smpl_fasta')                                                             
    chr1 = Chromosome('chr1','A'*500000)                                                  
    chr2 = Chromosome('chr2','C'*500000)                                               
    chr3 = Chromosome('chr3','G'*500000)                                               
    chr4 = Chromosome('chr4','T'*500000)                                                  
    fasta.add_chrom(chr1)                                                
    fasta.add_chrom(chr2)                                                
    fasta.add_chrom(chr3)                                                
    fasta.add_chrom(chr4)                                                
    fasta._add_nickname('chr1','CHR1')
    fasta._add_attribute('chr1','foobar')
    return fasta  
