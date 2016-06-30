"""
Create filters for given inheritance models.
See: https://github.com/arq5x/gemini/issues/388
"""
from __future__ import print_function
import sys
import os

from collections import defaultdict, Counter
import itertools as it
import operator as op
import re


try:
    reduce
except NameError:
    from functools import reduce

try:
    basestring
    zip = it.izip
except NameError:
    basestring = str
    unicode = str

def fix_sample_name(s):
    return s.replace("-", "_").replace(" ", "_")

def combine_and(*args):
    args = [a for a in args if not (a is True or a is None)]
    if any(a is False for a in args): return False

    cmd = None
    for i, a in enumerate(args):
        cmd = a if i == 0 else cmd & a
    return cmd

sex_lookup = {'1': 'male', '2': 'female',
               1: 'male', 2: 'female',
               'male': 'male', 'female': 'female'}

def make_classes(valid_gts, cfilter, HOM_REF, HET, UNKNOWN, HOM_ALT):
    class Sample(object):

        __slots__ = ('sample_id', 'name', 'affected', 'sex', 'mom', 'dad',
                     'family_id', '_i')

        valid_gts = None

        def __init__(self, sample_id, affected, sex=None, name=None,
                     family_id=None):
            #assert isinstance(sample_id, (long, int)), sample_id
            assert affected in (True, False, None)
            self.sample_id = sample_id
            self.name = name or sample_id
            self.affected = affected
            self.mom = None
            self.dad = None
            self.sex = sex_lookup.get(sex)
            self.family_id = family_id
            # _i is used to maintain the order in which they came in.
            self._i = None

        @property
        def genotype_lls(self):
            return [self.gt_phred_ll_homref,
                    self.gt_phred_ll_het,
                    self.gt_phred_ll_homalt]

        def __getattr__(self, gt_field):
            assert gt_field in valid_gts, gt_field
            return cfilter(self.sample_id, gt_field)

        def __repr__(self):
            c = self.__class__.__name__
            s = "%s(%s" % (c, self.name or self.sample_id)
            s += (";affected" if self.affected else (";unaffected"
                      if self.affected is False else ";unknown"))
            if self.sex is not None:
                s += ";%s" % self.sex
            return s + ")"

        def __str__(self):
            r = repr(self).split("(", 1)[1]
            return "%s(%s" % (self.name, r)
    Sample.valid_gts = valid_gts

    class Family(object):

        def __init__(self, subjects, fam_id):
            assert len(set(s.sample_id for s in subjects)) == len(subjects), subjects
            self.subjects = subjects
            self.family_id = fam_id
            # mostly for testing. this gets set when reading from db.
            if all(s._i is None for s in self.subjects):
                for i, s in enumerate(self.subjects):
                    s._i = i

        @classmethod
        def from_ped(klass, ped, order=None):
            """
            return a dict keyed by family_id with parent/kid relations defined from a ped file.
            """
            if isinstance(ped, basestring) and os.path.exists(ped):
                fh = open(ped)
            elif isinstance(ped, basestring):
                import io
                fh = io.StringIO(unicode(ped))
            else:
                fh = ped
            if order is not None:
                order = {fix_sample_name(k): v for k, v in order.items()}
            def agen():
                for toks in (l.rstrip().split() for l in fh if l[0] != "#"):
                    toks[1] = fix_sample_name(toks[1])
                    toks.append(toks[1]) # name
                    yield toks
            return klass._from_gen(agen(), order=order)

        def __len__(self):
            return len(self.subjects)

        def famphase(self, gt_types, gt_phases, gt_bases,
                     _splitter=re.compile("\||/"),
                     length_check=True):


            if length_check:
                assert len(self.subjects) == len(gt_types) == len(gt_bases)
            assert isinstance(gt_bases[0], basestring)
            #assert isinstance(gt_types[0], int) or int(gt_types[0]) == gt_types[0]

            for s in (subj for subj in self.subjects if subj.mom is not None and subj.dad is not None and gt_types[subj._i] == HET):
                # here we have a phaseable (HET) person with a mom and dad ..
                ## cant be same.
                if gt_types[s.mom._i] == gt_types[s.dad._i]: continue
                ## cant have unknown
                if UNKNOWN in (gt_types[s.mom._i], gt_types[s.dad._i]): continue

                kid_bases = set(_splitter.split(gt_bases[s._i]))
                mom_bases = _splitter.split(gt_bases[s.mom._i])
                dad_bases = _splitter.split(gt_bases[s.dad._i])

                parent_bases = set(mom_bases + dad_bases)
                # can't phase kid with de-novo

                if kid_bases - parent_bases:
                    sys.stderr.write("not phasing variant due to apparent de_novo in kid (%s) in family %s\n"
                                     % (s.name, self.family_id))
                    continue

                # no alleles from dad
                if len(kid_bases - set(dad_bases)) == len(kid_bases):
                    sys.stderr.write("not phasing variant due to no alleles in %s from dad %s (apparent mendelian error) in family %s\n"
                                      % (s.name, s.dad.name, self.family_id))
                    continue

                if len(kid_bases - set(mom_bases)) == len(kid_bases):
                    sys.stderr.write("not phasing variant due to no alleles in %s from mom %s (apparent mendelian error) in family %s\n"
                                     % (s.name, s.mom.name, self.family_id))
                    continue

                # should be able to phase here
                if gt_types[s.mom._i] in (HOM_REF, HOM_ALT):
                    assert gt_types[s.mom._i] in (HOM_REF, HOM_ALT)
                    mom_allele = mom_bases[0]
                    dad_alleles = dad_bases
                    dad_allele = next(d for d in dad_alleles if d != mom_allele)


                    gt_bases[s._i] = "%s|%s" % (mom_allele, dad_allele)
                else:
                    assert gt_types[s.dad._i] in (HOM_REF, HOM_ALT)

                    dad_allele = dad_bases[0]
                    mom_alleles = mom_bases
                    try:
                        mom_allele = next(m for m in mom_alleles if m != dad_allele)
                    except StopIteration:
                        gt_phases[s._i] = False
                        continue

                    gt_bases[s._i] = "%s|%s" % (mom_allele, dad_allele)

                gt_phases[s._i] = True

            return gt_phases, gt_bases

        def to_ped(self, fh=sys.stdout, header=True):
            if header:
                fh.write("#family_id sample_id paternal_id maternal_id sex phenotype\n")
            for s in self.subjects:
                paternal_id = (s.dad and s.dad.name) or "-9"
                maternal_id = (s.mom and s.mom.name) or "-9"
                phenotype = {True: '2', False: '1'}.get(s.affected, "-9")
                fh.write(" ".join((str(self.family_id), s.name, paternal_id, maternal_id, s.sex or '-9', phenotype)))
                fh.write("\n")

        @classmethod
        def from_cursor(klass, cursor):
            keys = "sample_id|family_id|name|paternal_id|maternal_id|sex|phenotype".split("|")
            def agen():
                for row in cursor.execute("select %s from samples" % ",".join(keys)):
                    if not isinstance(row, dict):
                        row = dict(zip(keys, row))
                    yield (row['family_id'], row['sample_id'], row['paternal_id'],
                           row['maternal_id'], row['sex'], str(row['phenotype']),
                           row['name'])
            return klass._from_gen(agen())

        @classmethod
        def _from_gen(klass, gen, order=None):
            fams = defaultdict(dict)
            pheno_lookup = {'1': False, '2': True}
            for i, line in enumerate(gen):
                (fam_id, indv, pat_id, mat_id, sex, pheno, name) = line[:7]

                assert indv not in fams[fam_id]
                s = fams[fam_id][name] = Sample(indv, pheno_lookup.get(pheno),
                                                sex=sex)
                s.mom = mat_id
                s.dad = pat_id
                # name in gemini is actually the id from the ped.
                # sample_id in gemini si the primary key
                s.name = name
                s.family_id = fam_id
                s._i = i if order is None else order[name]

            ofams = {}
            for fam_id, fam_dict in fams.items():
                ofams[fam_id] = []
                for name in fam_dict:
                    sample = fam_dict[name]
                    # convert sample_id to dad or None
                    sample.dad = fam_dict.get(sample.dad)
                    sample.mom = fam_dict.get(sample.mom)
                    ofams[fam_id].append(sample)

                ofams[fam_id] = Family(ofams[fam_id], fam_id)
                # maintain the order in which they came in.
                ofams[fam_id].subjects.sort(key=op.attrgetter('_i'))
            return ofams

        def __repr__(self):
            return "%s([%s])" % (self.__class__.__name__,
                               ", ".join(repr(s) for s in self.subjects))

        @property
        def males(self):
            return [s for s in self.subjects if s.sex in ('1', 1, 'male')]

        @property
        def females(self):
            return [s for s in self.subjects if s.sex in ('2', 2, 'female')]

        @property
        def affecteds(self):
            return [s for s in self.subjects if s.affected]

        @property
        def affecteds_with_parent(self):
            return [s for s in self.affecteds if not None in (s.mom, s.dad)]

        @property
        def samples_with_parent(self):
            return [s for s in self.subjects if not None in (s.mom, s.dad)]

        @property
        def unaffecteds(self):
            return [s for s in self.subjects if s.affected is False]

        @property
        def unknown(self):
            return [s for s in self.subjects if s.affected is None]

        # so we can do, e.g. fam.gts and get the list of required strings.
        def __getattr__(self, gt_field):
            if len(self.subjects) > 0:
                assert gt_field in self.subjects[0].valid_gts, gt_field
            return [getattr(s, gt_field) for s in self.subjects]

        def _restrict_to_min_depth(self, min_depth, unknowns=False):
            if min_depth is not None and min_depth > 0:
                if len(self.affecteds):
                    af = reduce(op.and_, [s.gt_depths >= min_depth for s in self.affecteds])
                else:
                    af = None
                if len(self.unaffecteds):
                    un = reduce(op.and_, [s.gt_depths >= min_depth for s in self.unaffecteds])
                else:
                    un = None
                if unknowns and len(self.unknown):
                    return combine_and(reduce(op.and_, [s.gt_depths >= min_depth for s in self.unknown]), af, un)
                return combine_and(af, un)
            else:
                return None

        def _restrict_to_min_gq(self, min_gq, unknowns=False):
            if min_gq is None or min_gq <= 0:
                return None

            if len(self.affecteds):
                af = reduce(op.and_, [s.gt_quals >= min_gq for s in self.affecteds])
            else:
                af = None
            if len(self.unaffecteds):
                un = reduce(op.and_, [s.gt_quals >= min_gq for s in self.unaffecteds])
            else:
                un = None
            if unknowns and len(self.unknown):
                return combine_and(reduce(op.and_, [s.gt_quals >= min_gq for s in self.unknown]), af, un)
            return combine_and(af, un)

        def auto_dom(self, min_depth=0, gt_ll=False, strict=True,
                only_affected=True, min_gq=0):
            """
            If strict then all affected kids must have at least 1 affected parent.
            parent.
            Parents of affected can't have unknown phenotype (for at least 1 kid)
            """
            if len(self.affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
                if strict:
                    return 'False'
            af = reduce(op.and_, [s.gt_types == HET for s in self.affecteds]) if len(self.affecteds) else True
            if len(self.unaffecteds) and only_affected:
                un = reduce(op.and_, [(s.gt_types != HET) & (s.gt_types != HOM_ALT) for s in self.unaffecteds])
            else:
                un = None
            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            if gt_ll:
                if len(self.affecteds):
                    af = reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in self.affecteds]) & af
                if len(self.unaffecteds) and only_affected:
                    un = reduce(op.and_, [s.gt_phred_ll_het > gt_ll for s in
                                           self.unaffecteds]) & un
            # need at least 1 kid with parent who has the mutation
            # parents can't have unkown phenotype.
            kid_with_known_parents = False
            # all affected kids must have at least 1 affected parent (or no parents)
            kid_with_parents = False
            for kid in self.affecteds:
                # if they have no parents, don't require it
                if kid.mom is None and kid.dad is None:
                    continue
                # mom or dad must be affected.
                kid_with_parents = True
                if strict and not any(p is not None and p.affected for p in (kid.mom, kid.dad)):
                    return 'False'
                # parents can't have unknown phenotype.
                if (kid.mom and kid.dad):
                    if (kid.mom.affected is not None) and (kid.dad.affected is not None):
                        kid_with_known_parents = True
                    # if he has a mom and dad that arent unknown, at least one of them must be affected
                    if not None in (kid.mom.affected, kid.dad.affected):
                        if not (kid.mom.affected or kid.dad.affected): return 'False'

            if strict and not kid_with_known_parents:
                return 'False'

            if not kid_with_parents:
                if len(self.affecteds) > 0:
                    sys.stderr.write("WARNING: using affected without parents for family %s for autosomal dominant test. Use strict to prevent this.\n" % self.family_id)
            return combine_and(af, un, depth, quals)

        def x_denovo(self, min_depth=0, min_gq=0):
            """
            #. affected female child must be het
            #. affected male child must be hom_alt (or het)
            #. parents should be unaffected and hom_ref
            """
            affecteds = self.affecteds
            if len(affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s. skipping\n" % self.family_id)
                return False

            female_af = [s.gt_types == HET for s in affecteds if s.sex == 'female']
            male_af = [(s.gt_types == HOM_ALT) | (s.gt_types == HET) for s in affecteds if s.sex == 'male']

            try:
                af = reduce(op.and_, female_af + male_af)
            except TypeError:
                assert len([x for x in affecteds if x.sex in ('male', 'female')]) == 0
                af = None

            un = None
            for kid in affecteds:
                if (kid.mom is not None and kid.mom.affected) or (kid.dad is not None and kid.dad.affected):
                    sys.stderr.write("WARNING: affected kid with affected parent in family: %s. not x-linked de novo. skipping\n" % self.family_id)
                    return 'False'
                for parent in (kid.mom, kid.dad):
                    if parent is None: continue
                    if un is None:
                        un = (parent.gt_types == HOM_REF)
                    else:
                        un &= (parent.gt_types == HOM_REF)

                if kid.mom is None or kid.dad is None:
                    sys.stderr.write("WARNING: running x-linked de novo on kid with no parents for family: %s\n" % self.family_id)

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            return combine_and(af, un, depth, quals)

        def x_dom(self, min_depth=0, min_gq=0):
            """
            X-linked dominant
            #. unaffecteds are hom_ref
            #. affected males must be het (PAR) or hom_alt
            #. affected females must be het
            #. boys of affected dad must be unaffected
            #. girls of affected dad must be affected
            """
            affecteds = self.affecteds
            if len(affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s. skipping\n" % self.family_id)
                return False

            try:
                # affected males het or hom_alt
                af = reduce(op.and_, [(s.gt_types == HET) | (s.gt_types == HOM_ALT) for s in affecteds if s.sex == 'male'])
                # affected females must be het
                af = reduce(op.and_, [s.gt_types == HET for s in affecteds if s.sex == 'female'], af)
            except TypeError:
                # affected females must be het
                try:
                    af = reduce(op.and_, [s.gt_types == HET for s in affecteds if s.sex == 'female'])
                except TypeError:
                    af = None

            # unaffecteds must be hom_ref
            if len(self.unaffecteds) > 0:
                un = reduce(op.and_, [s.gt_types == HOM_REF for s in self.unaffecteds])
            else:
                un = None

            def is_parent(a, b):
                if a is None or b is None: return False
                return a.sample_id == b.sample_id and a.family_id == b.family_id

            for parent in self.affecteds:
                kids = [k for k in self.subjects if is_parent(parent, k.dad)]
                for kid in kids:
                    # girls of affected dad must be affected
                    if kid.sex == 'female' and parent.sex == 'male' and not kid.affected:
                        sys.stderr.write("WARNING: unaffected female kid of affected dad "  +
                                "in family: %s. not x-linked dominant. skipping %s\n" % (kid.family_id, parent.name))
                        return False

                    #. boys of affected dad must be unaffected
                    if kid.sex == 'male' and parent.sex == 'male' and kid.affected:
                        sys.stderr.write("WARNING: affected male kid of affected dad in " +
                                "family: %s. not x-linked dominant. skipping %s\n" % (kid.family_id, parent.name))
                        return False

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            return combine_and(af, un, depth, quals)

        def x_rec(self, min_depth=0, min_gq=0):
            """
            X-linked recessive.
            #. affected females are hom_alt
            #. unaffected females are het or hom_ref
            #. affected males are het or hom_alt (issue warning for het males for PAR)

            NOTE: could add something for obligate carriers -- where unaffected females are het
            """
            affecteds = self.affecteds
            if len(affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
                return False

            try:
                female_af = [s.gt_types == HOM_ALT for s in affecteds if s.sex == 'female']
                female_af = reduce(op.and_, female_af)
            except TypeError:
                female_af = None
            try:
                female_un = reduce(op.and_, [s.gt_types == HET for s in self.unaffecteds if s.sex == 'female'])
            except TypeError:
                female_un = None

            try:
                male_af = reduce(op.and_, [(s.gt_types != UNKNOWN) & (s.gt_types != HOM_REF) for s in self.males if s.affected])
            except TypeError:
                male_af = None
            try:
                male_un = reduce(op.and_, [s.gt_types == HOM_REF for s in self.males if not s.affected])
            except TypeError:
                male_un = None

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            v = combine_and(female_af, female_un, male_af, male_un, depth, quals)
            return v

        def auto_rec(self, min_depth=0, gt_ll=False, strict=True,
                     only_affected=True, min_gq=0):
            """
            If strict, then if parents exist, they must be het for all affecteds
            """
            if len(self.affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
                return False
            af = reduce(op.and_, [s.gt_types == HOM_ALT for s in self.affecteds])
            if only_affected and len(self.unaffecteds) != 0:
                un = reduce(op.and_, [s.gt_types != HOM_ALT for s in self.unaffecteds])
            else:
                un = None
            if strict:
                # if parents exist, they must be het or affected for all affecteds
                # if both parents are not het then it's a de novo.
                usable_kid = False
                for kid in self.affecteds:
                    usable_kid = usable_kid or (kid.mom and kid.dad)
                    for parent in (kid.mom, kid.dad):
                        if parent is not None:
                            af &= parent.gt_types == HET
                            if parent.affected:
                                sys.stderr.write("WARNING: auto-recessive called on family "
                                        "%s where affected has affected parents\n" % self.family_id)
                                return "False"
                    if not usable_kid:
                        sys.stderr.write("WARNING: auto-recessive called on family "
                                "%s where no affected has parents\n" % self.family_id)

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            if gt_ll:
                af &= reduce(op.and_, [s.gt_phred_ll_homalt <= gt_ll for s in
                                       self.affecteds])
                if only_affected and len(self.unaffecteds) > 0:
                    un &= reduce(op.and_, [s.gt_phred_ll_homalt > gt_ll for s in
                                           self.unaffecteds])

            return combine_and(af, un, depth, quals)

        def de_novo(self, min_depth=0, gt_ll=False, strict=True,
                only_affected=True, min_gq=0):
            """
            all affected must be het.
            all unaffected must be homref.
            if strict, all affected kids must have unaffected parents.

            """
            if len(self.affecteds) == 0:
                sys.stderr.write("WARNING: no affecteds in family %s\n" % self.family_id)
                return 'False'
            af = reduce(op.and_, [s.gt_types == HET for s in self.affecteds])
            un = True
            have_un = len(self.unaffecteds) != 0
            if only_affected and have_un:
                un = reduce(op.and_, [s.gt_types == HOM_REF for s in self.unaffecteds])
            if gt_ll:
                af &= reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in self.affecteds])
                if only_affected and have_un:
                    un &= reduce(op.and_, [s.gt_phred_ll_homref <= gt_ll for s in self.unaffecteds])

            if only_affected and have_un:
                un2 = reduce(op.and_, [s.gt_types == HOM_ALT for s in
                                       self.unaffecteds])
                if gt_ll and have_un:
                    un2 &= reduce(op.and_, [s.gt_phred_ll_homalt <= gt_ll for s
                                            in self.unaffecteds])
                un |= un2

            # at least 1 affected kid must have unaffected parents
            un_parents = False
            for kid in self.affecteds:
                if kid.mom and kid.mom.affected is False and kid.dad and kid.dad.affected is False:
                    un_parents = True
                # can't have a parent with the variant
                for parent in (kid.mom, kid.dad):
                    if parent is None: continue
                    un = ((parent.gt_types == HOM_REF) | (parent.gt_types == HOM_ALT)) & un
                if kid.mom and kid.dad:
                    # otherwise, they could be HOM_REF, HOM_ALT and a het kid is not
                    # de_novo
                    un = (kid.mom.gt_types == kid.dad.gt_types) & un

            if not un_parents:
                return "False"

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            if strict:
                # if a parent is affected it's not de novo.
                for kid in self.affecteds:
                    for parent in (kid.mom, kid.dad):
                        if parent is not None and parent.affected:
                            return 'False'
            return combine_and(af, un, depth, quals)

        def mendel_plausible_denovo(self, min_depth=0, gt_ll=False,
                only_affected=False, min_gq=0):
            """
            kid == HET and dad, mom == HOM_REF or dad, mom == HOM_ALT
            only use kids with both parents present.
            """
            if only_affected:
                subset = self.affecteds_with_parent
            else:
                subset = self.samples_with_parent
            if len(subset) == 0: return 'False'

            depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
            quals = self._restrict_to_min_gq(min_gq, unknowns=not only_affected)
            # join exprs by or. so just look for any kid that meets these within a family.
            exprs = []
            for s in subset:
                # kid het, parents hom_ref
                expr = (s.gt_types == HET) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_REF)
                if gt_ll:
                    expr &= (s.gt_phred_ll_het <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
                exprs.append(expr)
                # kid het, parents hom_alt
                expr = (s.gt_types == HET) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HOM_ALT)
                if gt_ll:
                    expr &= (s.gt_phred_ll_het <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
                exprs.append(expr)
            return combine_and(reduce(op.or_, exprs), depth, quals)

        def mendel_implausible_denovo(self, min_depth=0, gt_ll=False,
                only_affected=False, min_gq=0):
            # everyone is homozygote. kid is opposit of parents.
            # only use kids with both parents present.
            if only_affected:
                subset = self.affecteds_with_parent
            else:
                subset = self.samples_with_parent
            if len(subset) == 0: return 'False'

            depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
            quals = self._restrict_to_min_gq(min_gq, unknowns=not only_affected)
            exprs = []
            for s in subset:
                # kid hom_alt, parents hom_ref
                expr = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_REF)
                if gt_ll:
                    expr &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
                exprs.append(expr)
                # parents hom_alt kid homref
                expr = (s.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HOM_ALT)
                if gt_ll:
                    expr &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
                exprs.append(expr)
            return combine_and(reduce(op.or_, exprs), depth, quals)

        def mendel_uniparental_disomy(self, min_depth=0, gt_ll=False,
                                      only_affected=False, min_gq=0):
            # parents are opposite homs, kid matches one of them (but should be
            # het).
            if only_affected:
                subset = self.affecteds_with_parent
            else:
                subset = self.samples_with_parent
            if len(subset) == 0: return 'False'
            depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
            quals = self._restrict_to_min_gq(min_gq, unknowns=not only_affected)

            # join exprs with or
            exprs = []
            for s in subset:
                for gtkid in (HOM_REF, HOM_ALT):
                    # mom homref, dad hom_alt.
                    expr = (s.gt_types == gtkid) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HOM_ALT)
                    if gt_ll:
                        if gtkid == HOM_REF:
                            expr &= (s.gt_phred_ll_homref <= gt_ll)
                        else:
                            expr &= (s.gt_phred_ll_homalt <= gt_ll)
                        expr &= (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
                    exprs.append(expr)

                    # mom homalt, dad hom_ref.
                    expr = (s.gt_types == gtkid) & (s.dad.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT)
                    if gt_ll:
                        if gtkid == HOM_REF:
                            expr &= (s.gt_phred_ll_homref <= gt_ll)
                        else:
                            expr &= (s.gt_phred_ll_homalt <= gt_ll)
                        expr &= (s.dad.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll)
                    exprs.append(expr)

            return combine_and(reduce(op.or_, exprs), depth, quals)

        def mendel_LOH(self, min_depth=0, gt_ll=False, only_affected=False,
                min_gq=0):
            # kid and one parent are opposite homozygotes other parent is het.
            if only_affected:
                subset = self.affecteds_with_parent
            else:
                subset = self.samples_with_parent
            if len(subset) == 0: return 'False'

            depth = self._restrict_to_min_depth(min_depth, unknowns=not only_affected)
            quals = self._restrict_to_min_gq(min_gq, unknowns=not only_affected)
            exprs = []  # joined with or
            for s in subset:
                # kid hom_alt, mom hom_ref, dad het.
                e = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HOM_REF) & (s.dad.gt_types == HET)
                if gt_ll:
                    e &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_homref <= gt_ll) & (s.dad.gt_phred_ll_het <= gt_ll)
                exprs.append(e)

                # kid hom_ref, mom het, dad hom_alt
                e = (s.gt_types == HOM_REF) & (s.mom.gt_types == HET) & (s.dad.gt_types == HOM_ALT)
                if gt_ll:
                    e &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_het <= gt_ll) & (s.dad.gt_phred_ll_homalt <= gt_ll)
                exprs.append(e)

                # kid hom_alt, mom het, dad hom_ref
                e = (s.gt_types == HOM_ALT) & (s.mom.gt_types == HET) & (s.dad.gt_types == HOM_REF)
                if gt_ll:
                    e &= (s.gt_phred_ll_homalt <= gt_ll) & (s.mom.gt_phred_ll_het <= gt_ll) & (s.dad.gt_phred_ll_homref <= gt_ll)
                exprs.append(e)

                # kid hom_ref, mom hom_alt, dad het
                e = (s.gt_types == HOM_REF) & (s.mom.gt_types == HOM_ALT) & (s.dad.gt_types == HET)
                if gt_ll:
                    e &= (s.gt_phred_ll_homref <= gt_ll) & (s.mom.gt_phred_ll_homalt <= gt_ll) & (s.dad.gt_phred_ll_het <= gt_ll)
                exprs.append(e)

            return combine_and(reduce(op.or_, exprs), depth, quals)

        def mendel_violations(self, min_depth=0, gt_ll=False,
                              only_affected=False, min_gq=0):
            return {'plausible de novo': self.mendel_plausible_denovo(min_depth,
                                                                      gt_ll,
                                                                      only_affected,
                                                                      min_gq=min_gq),
                    'implausible de novo': self.mendel_implausible_denovo(min_depth,
                                                                          gt_ll,
                                                                          only_affected,
                                                                          min_gq=min_gq),
                    'uniparental disomy': self.mendel_uniparental_disomy(min_depth,
                                                                         gt_ll,
                                                                         only_affected,
                                                                         min_gq=min_gq),
                    'loss of heterozygosity': self.mendel_LOH(min_depth, gt_ll,
                                                              only_affected,
                                                              min_gq=min_gq)
                    }

        def _get_ref_alt(self, gt_types, gt_bases,
                         _splitter=re.compile("\||/")):
            """
            Guess the ref and alt. Mostly for convenience for comp_het functions,
            as we should know these anyway.
            """
            s = self.subjects[0]
            HOM_REF, HET, HOM_ALT = self.HOM_REF, self.HET, self.HOM_ALT
            ref, alt = None, None
            for i, gt in enumerate(gt_types):
                if gt == HOM_REF:
                    ref = _splitter.split(gt_bases[i])[0]
                elif gt == HOM_ALT:
                    alt = _splitter.split(gt_bases[i])[0]
                elif "/" in gt_bases[i]:
                    _ref, _alt = gt_bases[i].split("/")
                    if ref is None:
                        ref = _ref
                    if alt is None and _ref != _alt:
                        alt = _alt
            # fall back to allele frequency
            if ref is None or alt is None or ref == alt:
                c = Counter()
                for b in gt_bases:
                    c.update(_splitter.split(b))
                if ref is None:
                    ref = c.most_common(1)[0][0]
                    if ref == alt:
                        ref = c.most_common(2)[1][0]
                if alt is None:
                    alt = c.most_common(2)[1][0]
                    if ref == alt:
                        alt = c.most_common(1)[0][0]
            return ref, alt

        def _comp_het_pair_pattern(self,
                                   gt_types1, gt_nums1,
                                   gt_types2, gt_nums2,
                                   gt_phases1, gt_phases2):
            """
            + kid has to be phased het at both sites.
            + kid has to have alts on different chroms.
            + neither parent can be hom_alt at either site.
            + if either parent is phased at both sites and matches the kid, exclude.
            + if either parent is het at both sites, priority is reduced
            """

            # already phased before sending here.
            ret = {'candidates': [], 'priority': 4}
            for kid in self.samples_with_parent:
                if gt_nums1[kid._i] == gt_nums2[kid._i]: continue
                if not (gt_types1[kid._i] == HET and gt_types2[kid._i] == HET): continue
                #if not (gt_phases1[kid._i] and gt_phases2[kid._i]): continue
                if gt_types1[kid.mom._i] == HOM_ALT or gt_types2[kid.dad._i] == HOM_ALT: continue
                mom, dad = kid.mom, kid.dad

                kid_phased = gt_phases1[kid._i] and gt_phases2[kid._i]
                dad_phased = gt_phases1[dad._i] and gt_phases2[dad._i]
                mom_phased = gt_phases1[mom._i] and gt_phases2[mom._i]

                if kid_phased and dad_phased and (gt_nums1[dad._i] == gt_nums1[kid._i]) and (gt_nums2[dad._i] == gt_nums2[kid._i]):
                    continue
                if kid_phased and mom_phased and (gt_nums1[mom._i] == gt_nums1[kid._i]) and (gt_nums2[mom._i] == gt_nums2[kid._i]):
                    continue

                if kid_phased and dad_phased and mom_phased and gt_types1[dad._i] != gt_types2[dad._i] and gt_types1[mom._i] != gt_types2[mom._i]:
                    priority = 1

                elif kid_phased and gt_types1[dad._i] != gt_types1[mom._i] and gt_types2[dad._i] != gt_types2[mom._i]:
                    # parents are unphased hets at different sites.
                    priority = 1
                else:
                    priority = 2
                    for parent in (kid.mom, kid.dad):
                        # unphased het
                        if gt_types2[parent._i] == gt_types1[parent._i] == HET:
                            priority += 1

                ret['candidates'].append(kid)
                ret['priority'] = min(ret['priority'], priority)
            ret['candidate'] = len(ret['candidates']) > 0
            return ret


        def comp_het_pair(self, gt_types1, gt_bases1,
                          gt_types2, gt_bases2,
                          gt_phases1=None,
                          gt_phases2=None,
                          ref1=None, alt1=None,
                          ref2=None, alt2=None,
                          allow_unaffected=False,
                          fast_mode=False,
                          pattern_only=False):
            """
            Each of the sites here must have passed the comp_het() filter.
            This further checks that a give pair is comp_het.

            if pattern_only is False, affected/unaffected status is ignored.
            A priority is returned such that:
            1. It's a true, canonical compound het where everyone is phased (or
               it's unambiguous
            2. It's a singleton (phased or unphased) sample with no parents
            3. It's an unphased sample and unphased parents where all are hets
               at the pair.
            """
            if gt_phases1 is None:
                gt_phases1 = ["|" in b for b in gt_bases1]
            if gt_phases2 is None:
                gt_phases2 = ["|" in b for b in gt_bases2]

            if ref1 is None and alt1 is None:
                ref1, alt1 = self._get_ref_alt(gt_types1, gt_bases1)

            if ref2 is None and alt2 is None:
                ref2, alt2 = self._get_ref_alt(gt_types2, gt_bases2)

            idxs = set(s._i for s in self.subjects)

            self.famphase(gt_types1, gt_phases1, gt_bases1,
                          length_check=False)
            self.famphase(gt_types2, gt_phases2, gt_bases2,
                          length_check=False)

            # we index by sample._i, but we only need to split for the current
            # samples so we check if i in idxs.
            gt_bases1 = [b.split("|" if p else "/") if i in idxs else b for i,
                    (b, p) in enumerate(zip(gt_bases1, gt_phases1))]
            gt_bases2 = [b.split("|" if p else "/") if i in idxs else b for i,
                    (b, p) in enumerate(zip(gt_bases2, gt_phases2))]

            # get in (0, 1) format instead of (A, T)
            try:
                ra = {ref1: 0, alt1: 1, ".": 2}
                gt_nums1 = [(ra[b[0]], ra[b[1]]) if i in idxs else None for i, b in enumerate(gt_bases1)]
                ra = {ref2: 0, alt2: 1, ".": 2}
                gt_nums2 = [(ra[b[0]], ra[b[1]]) if i in idxs else None for i, b in enumerate(gt_bases2)]
            except KeyError:
                sys.stderr.write("can't phase sites with multiple alternate alleles\n")
                return {'candidate': False}
            except IndexError:
                # alternate is unknown, e.g. "A/."
                return {'candidate': False}

            if pattern_only:
                return self._comp_het_pair_pattern(gt_types1, gt_nums1,
                                                   gt_types2, gt_nums2,
                                                   gt_phases1, gt_phases2)

            unaffecteds = self.unaffecteds
            for un in unaffecteds:
                if gt_types2[un._i] == HOM_ALT or gt_types1[un._i] == HOM_ALT:
                    return {'candidate': False}

            ret = {'affected_phased': [], 'unaffected_phased': [],
                   'unaffected_unphased': [], 'affected_unphased': [],
                   'affected_skipped': [], 'candidates': []}

            aff = None
            for aff in self.affecteds:
                if gt_types1[aff._i] != HET or gt_types2[aff._i] != HET:
                    ret['affected_skipped'].append(aff)
                    ret['candidate'] = False
                    if fast_mode: break
                    continue

                aff_phased = gt_phases1[aff._i] and gt_phases2[aff._i]
                # on same chrom.
                if aff_phased and gt_nums1[aff._i] == gt_nums2[aff._i]:
                    ret['affected_skipped'].append(aff)
                    # Remove candidates where an affected from the same family does
                    # NOT share the same het pair.
                    ret['candidate'] = False
                    if fast_mode: break
                    continue

                if not 'candidate' in ret: ret['candidate'] = True
                if aff_phased:
                    ret['affected_phased'].append(aff)
                else:
                    ret['affected_unphased'].append(aff)
                ret['candidates'].append(aff)

            del aff
            if ret['candidates'] != []:

                for un in unaffecteds:
                    if gt_types1[un._i] != HET or gt_types2[un._i] != HET:
                        continue

                    is_phased = gt_phases1[un._i] and gt_phases2[un._i]
                    # unaffected has the candidate pair on the same chromosome
                    if is_phased and gt_nums1[un._i] == gt_nums2[un._i]:
                        continue

                    if is_phased:
                        # found an unaffected with the same het-pair.
                        ret['unaffected_phased'].append(un)
                        if not allow_unaffected:
                            ret['candidate'] = False
                            if fast_mode: break
                    else:
                        ret['unaffected_unphased'].append(un)
            if not 'candidate' in ret:
                ret['candidate'] = False
                ret['priority'] = None
            elif ret['candidate']:

                ret['priority'] = 3
                if len(ret['affected_phased']) and len(ret['unaffected_unphased']) == 0:
                    ret['priority'] = 1
                    # priority 2 for a single unphased affected.
                elif len(ret['affected_unphased']) and len(ret['unaffected_unphased']) == 0:
                    ret['priority'] = 2

            return ret

        def comp_het(self, min_depth=0, gt_ll=False,
                     only_affected=True,
                     pattern_only=False,
                     min_gq=0):

            if pattern_only:
                af, un = None, None
                for i, kid in enumerate(self.samples_with_parent):

                    icmp = (kid.gt_types == HET) & (kid.mom.gt_types != HOM_ALT) \
                                                 & (kid.dad.gt_types != HOM_ALT) \
                                                 & (kid.mom.gt_types != UNKNOWN) \
                                                 & (kid.dad.gt_types != UNKNOWN)
                    if i == 0:
                        af = icmp
                    else:
                        af |= icmp
            else:
                # all affecteds must be het at both sites
                af = None
                if len(self.affecteds):
                    af = reduce(op.or_, [s.gt_types == HET for s in self.affecteds])

                # no unaffected can be homozygous alt at either site.
                un = None
                if len(self.unaffecteds):
                    un = reduce(op.and_, [s.gt_types != HOM_ALT for s in self.unaffecteds])

                for kid in self.samples_with_parent:
                    if not kid.affected: continue
                    un = (kid.mom.gt_types != UNKNOWN) & un
                    un = (kid.dad.gt_types != UNKNOWN) & un

                if gt_ll:
                    if len(self.affecteds):
                        af &= reduce(op.and_, [s.gt_phred_ll_het <= gt_ll for s in
                            self.affecteds])
                    if len(self.unaffecteds):
                        un = reduce(op.and_, [s.gt_phred_ll_homalt > gt_ll for s in
                            self.unaffecteds], un)

            depth = self._restrict_to_min_depth(min_depth)
            quals = self._restrict_to_min_gq(min_gq)
            return combine_and(af, un, depth, quals)

    Family.HOM_REF, Family.HET, Family.UNKNOWN, Family.HOM_ALT = HOM_REF, HET, UNKNOWN, HOM_ALT
    return Sample, Family
