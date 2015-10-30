"""

# gemini implementation uses 0, 1, 2, 3, 4 for HOM_REF ... HOM_ALT, but we use
# the strings here for illustrative testing.
>>> HOM_REF, HET, UNKNOWN, HOM_ALT = ('HOM_REF', 'HET', 'UNKNOWN', 'HOM_ALT')
>>> Sample, Family = make_classes(valid_gts, Filter, HOM_REF, HET, UNKNOWN, HOM_ALT)

>>> a, b = Sample(1, True), Sample(2, True)
>>> a.gt_types == HOM_REF
'gt_types[0] == HOM_REF'

>>> a.gt_types == b.gt_types
'gt_types[0] == gt_types[1]'

>>> a.gt_types != b.gt_types
'gt_types[0] != gt_types[1]'

>>> (a.gt_phred_ll_homref < 2) & (a.gt_phred_ll_homalt > 2)
'(gt_phred_ll_homref[0] < 2) and (gt_phred_ll_homalt[0] > 2)'

>>> (a.gt_phred_ll_homref < 2) | (b.gt_phred_ll_homalt > 2)
'(gt_phred_ll_homref[0] < 2) or (gt_phred_ll_homalt[1] > 2)'

>>> mom = Sample(1, affected=False)
>>> dad = Sample(2, affected=False)
>>> kid = Sample(3, affected=True)
>>> Family([mom, dad, kid], 'a').auto_rec()
'(gt_types[2] == HOM_ALT) and ((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT))'

# unknowns don't count
>>> kid2 = Sample(4, affected=None)

>>> Family([mom, dad, kid, kid2], 'a').auto_rec()
'(gt_types[2] == HOM_ALT) and ((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT))'

>>> Family([mom, dad, kid, kid2], 'a').auto_rec(gt_ll=1, min_depth=10)
'(((gt_types[2] == HOM_ALT) and (gt_phred_ll_homalt[2] <= 1)) and (((gt_types[0] != HOM_ALT) and (gt_types[1] != HOM_ALT)) and ((gt_phred_ll_homalt[0] > 1) and (gt_phred_ll_homalt[1] > 1)))) and ((gt_depths[2] >= 10) and ((gt_depths[0] >= 10) and (gt_depths[1] >= 10)))'

>>> f = Family([Sample("mom", False), Sample("dad", False),
...             Sample("kid", True)], "fam")
>>> f.subjects[2].mom = f.subjects[0]
>>> f.subjects[2].dad = f.subjects[1]
>>> f.famphase([HOM_REF, HOM_ALT, HET],
...            [False, False, False],
...            ["A/A", "A/C", "C/A"])
([False, False, True], ['A/A', 'A/C', 'A|C'])

>>> f.famphase([HOM_REF, HOM_REF, HET],
...            [False, False, False],
...            ["A/A", "A/C", "C/A"])
([False, False, False], ['A/A', 'A/C', 'C/A'])

>>> f.famphase([HOM_REF, HET, HET],
...            [False, False, False],
...            ["A/A", "A/C", "C/A"])
([False, False, True], ['A/A', 'A/C', 'A|C'])

>>> f.famphase([HET, HET, HET],
...            [False, False, False],
...            ["A/C", "A/C", "A/C"])
([False, False, False], ['A/C', 'A/C', 'A/C'])

>>> f.famphase([HOM_REF, HET, HET],
...            [False, False, False],
...            ["AA/A", "AA/C", "C/AA"])
([False, False, True], ['AA/A', 'AA/C', 'AA|C'])

>>> f.famphase([HOM_REF, HET, HET],
...            [False, False, False],
...            ['G/G', 'AA/C', 'A/C'])
([False, False, False], ['G/G', 'AA/C', 'A/C'])

>>> f.famphase([HOM_REF, HET, HET],
...            [False, False, False],
...            ['G/G', 'A/C', 'A/C'])
([False, False, False], ['G/G', 'A/C', 'A/C'])

>>> f = Family([Sample("mom", False), Sample("dad", False),
...             Sample("kid", True)], "fam")
>>> f.subjects[2].mom = f.subjects[0]
>>> f.subjects[2].dad = f.subjects[1]
>>> r = f.mendel_violations()
>>> for k in r:
...     print(k)
...     print(r[k])
uniparental disomy
(((((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_ALT)) or (((gt_types[kid] == HOM_REF) and (gt_types[dad] == HOM_REF)) and (gt_types[mom] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[dad] == HOM_REF)) and (gt_types[mom] == HOM_ALT))
plausible de novo
(((gt_types[kid] == HET) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_REF)) or (((gt_types[kid] == HET) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HOM_ALT))
loss of heterozygosity
(((((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HET)) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HET)) and (gt_types[dad] == HOM_ALT))) or (((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HET)) and (gt_types[dad] == HOM_REF))) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HET))
implausible de novo
(((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HOM_REF)) and (gt_types[dad] == HOM_REF)) or (((gt_types[kid] == HOM_REF) and (gt_types[mom] == HOM_ALT)) and (gt_types[dad] == HOM_ALT))

>>> f = Family([Sample("mom", False), Sample("dad", False),
...             Sample("kid", True)], "fam")
>>> f.gt_types = [HOM_ALT, HOM_ALT, HET]
>>> f.subjects[2].mom = f.subjects[0]
>>> f.subjects[2].dad = f.subjects[1]

>>> f.de_novo() #doctest: +ELLIPSIS
'(gt_types[kid] == HET...

>>> f.auto_rec()
'(((gt_types[kid] == HOM_ALT) and (gt_types[mom] == HET)) and (gt_types[dad] == HET)) and ((gt_types[mom] != HOM_ALT) and (gt_types[dad] != HOM_ALT))'

>>> f.auto_dom()
'False'
>>> f.auto_dom(strict=False)
'False'

# need an affected parent for autosomal dominant
>>> f.subjects[1].affected = True
>>> f.auto_dom(strict=False)
'((gt_types[dad] == HET) and (gt_types[kid] == HET)) and ((gt_types[mom] != HET) and (gt_types[mom] != HOM_ALT))'

>>> f.subjects[1].affected = False
>>> f.comp_het() #doctest: +ELLIPSIS
'(gt_types[kid] == HET)...'

>>> f2 = Family([Sample("mom", False), Sample("dad", False),
...              Sample("kid", True), Sample("kid2", False)], "fam")
>>> f2.subjects[2].mom = f.subjects[0]
>>> f2.subjects[2].dad = f.subjects[1]
>>> f2.subjects[3].mom = f.subjects[0]
>>> f2.subjects[3].dad = f.subjects[1]

>>> f2.gt_types = [HOM_REF, HOM_REF, HET, HET]

#>>> f2.de_novo()
#'False'

#>>> f2.de_novo(only_affected=False, strict=False)
#'True'
"""
# a specific implementation of inheritance models used by gemini that uses
# python eval statements
try:
    reduce
except NameError:
    from functools import reduce

from .inheritance import make_classes


valid_gts = (
    'gts',
    'gt_types',
    'gt_phases',
    'gt_depths',
    'gt_ref_depths',
    'gt_alt_depths',
    'gt_quals',
    'gt_copy_numbers',
    'gt_phred_ll_homref',
    'gt_phred_ll_het',
    'gt_phred_ll_homalt',
)

try:
    ints = (long, int)
except NameError:
    ints = (int,)

class Filter(object):
    def __init__(self, sample_id, gt_field):
        self.gt_field = gt_field
        if isinstance(sample_id, ints):
            self.sample0 = sample_id - 1
        else:
            self.sample0 = sample_id

    def in_(self, li):
        return reduce(op.or_, [self == i for i in li])

    def __lt__(self, o):
        return ostr("%s[%s] < %s" % (self.gt_field, self.sample0, o))

    def __le__(self, o):
        return ostr("%s[%s] <= %s" % (self.gt_field, self.sample0, o))

    def __gt__(self, o):
        return ostr("%s[%s] > %s" % (self.gt_field, self.sample0, o))

    def __ge__(self, o):
        return ostr("%s[%s] >= %s" % (self.gt_field, self.sample0, o))

    def __eq__(self, o):
        return ostr("%s[%s] == %s" % (self.gt_field, self.sample0, o))

    def __ne__(self, o):
        return ostr("%s[%s] != %s" % (self.gt_field, self.sample0, o))

    def __str__(self):
        return ostr("%s[%s]" % (self.gt_field, self.sample0))

    def __repr__(self):
        return ostr("%s[%s]" % (self.gt_field, self.sample0))

    def __and__(self, other):
        raise Exception("shouldn't be here. wrap &/| statements in parens")

    __or__ = __and__

    def __nonzero__(self):
        raise Exception("shouldn't be here. use & instead of 'and'")




def _bracket(other):
    o = str(other)
    if o in ("True", "False"): return o
    return "(%s)" % o

class ostr(str):
    def __and__(self, other):
        if other is None: return self
        if other is True: return self
        if other is False: return self
        return ostr("%s and %s" % (_bracket(self), _bracket(other)))

    def __or__(self, other):
        if other is None: return self
        if other is True: return True
        if other is False: return False
        return ostr("%s or %s" % (_bracket(self), _bracket(other)))

    def __nonzero__(self):
        raise Exception("shouldn't be here. use & instead of 'and'. and wrap in parens")

Sample, Family = make_classes(valid_gts, Filter, *range(4))
