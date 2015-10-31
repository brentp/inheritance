
from inheritance import Sample, Family, EvalFamily


mom = Sample('mom', affected=False)
dad = Sample('dad', affected=False)
kid = Sample('kid', affected=True)

kid.mom, kid.dad = mom, dad

fam = Family([mom, dad, kid], 'a')

def make_fam1():
    # only 1 affected kid.
    fam = Family.from_ped("""\
#family_id  sample_id   paternal_id maternal_id sex phenotype
1   dad   0   0   1  1
1   mom   grandpa   grandma   2  1
1   kid   dad   mom   1  2
1   kid2   dad   mom   1  1
1   grandma 0   0     2  1
1   grandpa 0   0     1  1""")
    return fam

def make_fam2():
    # 1 affected kid, parent, grandparent
    fam = Family.from_ped("""\
#family_id  sample_id   paternal_id maternal_id sex phenotype
1   dad   0   0   1  1
1   mom   grandpa   grandma   2  2
1   kid   dad   mom   1  2
1   kid2   dad   mom   1  1
1   grandma 0   0     2  2
1   grandpa 0   0     1  1""")
    return fam



def test_fam():
    assert fam.subjects == [mom, dad, kid]

def test_samples():
    assert repr(mom) != "", repr(mom)

def test_attrs():
    assert mom.affected is False
    assert kid.affected is True

    assert fam.family_id == 'a', fam.family_id

    assert dad.sample_id == 'dad'

def test_auto_rec():

    assert "gt_types[kid] == 3" in fam.auto_rec()

    efam = EvalFamily(fam)
    efam.gt_types = [Family.HET, Family.HET, Family.HOM_ALT]

    assert efam.auto_rec()

    efam.gt_types = [Family.HET, Family.HET, Family.HET]
    assert not efam.auto_rec()


def test_auto_rec_kid_unaffected():
    kid.affected = False
    efam = EvalFamily(fam)
    efam.gt_types = [Family.HET, Family.HET, Family.HOM_ALT]
    assert not efam.auto_rec()
    kid.affected = True


def test_auto_rec_extended():
    fam = make_fam1()

    efam = EvalFamily(fam)
    efam.gt_types = [Family.HET, Family.HET, Family.HOM_ALT, Family.HET, Family.HET, Family.HET]
    assert efam.auto_rec()

    #if grandpa is affected it is no longer autosomal recessive
    efam.subjects[5].affected = True
    assert not efam.auto_rec()
    efam.subjects[5].affected = False

    # set both kids to hom_alt:
    efam.gt_types[3] = Family.HOM_ALT
    assert not efam.auto_rec()
    assert efam.auto_rec(only_affected=False)

    # set unaffected kid back:
    efam.gt_types[3] = Family.HET
    # expected that is is auto_rec
    assert efam.auto_rec()
    # but it's not if we set all depths to 9 and have min_depth of 10
    efam.gt_depths = [9] * 6
    assert not efam.auto_rec(min_depth=10)

    efam.gt_depths[2] = 1000
    # even if affected has high depth
    assert not efam.auto_rec(min_depth=10)

    # passes if we meet the minimum depth
    efam.gt_depths = [100] * 6
    assert efam.auto_rec(min_depth=10)

    # check if we have no affecteds:
    efam.subjects[2].affected = False
    assert not efam.auto_rec(min_depth=10)

def test_auto_dom_extended():

    fam = make_fam2()
    efam = EvalFamily(fam)
    assert [f.affected for f in efam.subjects] == [False, True, True, False, True, False]

    efam.gt_types = [Family.HOM_REF, Family.HET, Family.HET, Family.HOM_REF, Family.HET, Family.HOM_REF]

    assert efam.auto_dom()

    # unaffected is het:
    efam.gt_types[0] = Family.HET
    assert not efam.auto_dom()
    assert efam.auto_dom(only_affected=False)
    efam.gt_types[0] = Family.HOM_REF

    assert efam.auto_dom()
    efam.gt_depths = [9] * 6
    assert not efam.auto_dom(min_depth=10)

    # check with no unaffecteds
    for f in efam.subjects:
        f.affected = True
    assert not efam.auto_dom()

    # check with all unaffecteds
    for f in efam.subjects:
        f.affected = False
    assert not efam.auto_dom()

def test_denovo():
    fam = make_fam2()
    efam = EvalFamily(fam)
    assert [f.affected for f in efam.subjects] == [False, True, True, False, True, False]

    efam.gt_types = [Family.HOM_REF, Family.HET, Family.HET, Family.HOM_REF, Family.HET, Family.HOM_REF]
    assert not efam.de_novo()

def test_comphet_pair():


    fam = make_fam2()
    efam = EvalFamily(fam)
    assert [f.affected for f in efam.subjects] == [False, True, True, False, True, False]
    efam.gt_types = [Family.HOM_REF, Family.HET, Family.HET, Family.HOM_REF, Family.HET, Family.HOM_REF]

    res = efam.comp_het()
    assert res

    # still ok with het
    efam.gt_types[0] = Family.HET
    res = efam.comp_het()
    assert res

    # not ok with hom_alt of unaffected
    efam.gt_types[0] = Family.HOM_ALT
    res = efam.comp_het()
    assert not res


def test_comphet_pattern():
    fam = make_fam2()
    efam = EvalFamily(fam)
    assert [f.affected for f in efam.subjects] == [False, True, True, False, True, False]
    gt_types1 = [Family.HOM_REF, Family.HET, Family.HET, Family.HOM_REF, Family.HET, Family.HOM_REF]
    gt_types2 = [Family.HOM_REF, Family.HET, Family.HET, Family.HOM_REF, Family.HET, Family.HOM_REF]
    efam.gt_types = gt_types1
    gt_bases1 = ["A/C", "A/A", "A/C", "A/A", "A/A", "A/C", "A/A"]
    gt_bases2 = ["A/C", "A/A", "A/C", "A/A", "A/A", "A/C", "A/A"]
    res = efam.comp_het_pair(gt_types1, gt_bases1, gt_types2, gt_bases2)
