from __future__ import print_function

import os
import sys
import tempfile
import atexit

import itertools as it

try:
    basestring
    ints = (int, long)
except NameError:
    basestring = str
    ints = int

def tmp(pedstr, suf=".ped"):
    t = tempfile.mktemp(suffix=suf)
    atexit.register(os.unlink, t)
    with open(t, "w") as fh:
        for line in pedstr.split("\n"):
            if not line.strip(): continue
            print(line.strip(), file=fh)
    return t


class EvalFamily(object):

    __slots__ = ('ped', 'family', '_gt_types',
                 '_gt_phred_ll_homref',
                 '_gt_phred_ll_het',
                 '_gt_phred_ll_homalt',
                 '_gt_quals',
                 '_gt_depths', 'strict', 'subjects')

    def __init__(self, family, fam_id=None, gt_types=None, gt_depths=None):
        # can send in a family.
        if isinstance(family, dict):
            if len(family) != 1:
                raise Exception("only single families supported in EvalFamily")
            family = next(iter(family.values()))
        self.family = family
        for s in self.family.subjects:
            if s.sample_id[0].isdigit(): s.sample_id = "s" + s.sample_id
        self.subjects = self.family.subjects

        self._gt_types = None
        self.gt_types = gt_types
        self._gt_depths = None
        self.gt_depths = gt_depths
        self._gt_phred_ll_homalt = self._gt_phred_ll_homref = self._gt_phred_ll_het = None
        self._gt_quals = None

    def draw(self, tests=('auto_rec', 'auto_dom')):
        from IPython.display import Image, display
        if isinstance(tests, basestring):
            tests = (tests,)
        img = self.dot(tests=tests)
        return display(Image(filename=img))

    def dot(self, comment=None, path="test.gv", view=False, tests=('auto_rec', 'auto_dom')):
        from graphviz import Digraph
        viz = Digraph(comment=comment)
        subjects = self.family.subjects
        lookup = ["HOM_REF", "HET", "UNKOWN", "HOM_ALT"]
        for i, s in enumerate(subjects):

            attrs = dict(style="filled", fontcolor="white")
            attrs["fillcolor"] = {True: 'black', False: 'white', None: 'gray'}[s.affected]
            attrs["shape"] = {'male': 'square', 'female': 'circle', None: 'octagon'}[s.gender]

            if attrs["fillcolor"] == "black":
                attrs["fontcolor"] = "white"
            elif attrs["fillcolor"] == "white":
                attrs["fontcolor"] = "black"

            gt = lookup[self.gt_types[i]]
            label = s.name
            viz.node(s.name, label + "\n" + gt, **attrs)
        for s in subjects:
            if s.dad is not None:
                viz.edge(s.dad.name, s.name)
            if s.mom is not None:
                viz.edge(s.mom.name, s.name)
        for test in tests:
            res = {}
            res['default'] = getattr(self, test)()
            res['strict=False'] = getattr(self, test)(strict=False)
            res['only_affected=False'] = getattr(self, test)(only_affected=False)
            res['both False'] = getattr(self, test)(only_affected=False, strict=False)
            print("\n" + test)
            print("-" * len(test))
            for k in ("default", "strict=False", "only_affected=False", "both False"):
                print("%-20s\t%s" % (k, res[k]))

        viz._format = "png"
        return viz.render(path, view=view)

    @property
    def gt_types(self):
        return self._gt_types

    @gt_types.setter
    def gt_types(self, gt_types):
        if gt_types is not None:
            assert len(gt_types) == len(self.family)
            self._gt_types = gt_types

    @property
    def gt_quals(self):
        return self._gt_quals

    @gt_quals.setter
    def gt_quals(self, gt_quals):
        if gt_quals is not None:
            assert len(gt_quals) == len(self.family)
            self._gt_quals = gt_quals

    @property
    def gt_depths(self):
        return self._gt_depths

    @gt_depths.setter
    def gt_depths(self, gt_depths):
        if gt_depths is not None:
            assert len(gt_depths) == len(self.family)
            self._gt_depths = gt_depths

    @property
    def gt_phred_ll_homref(self):
        return self._gt_phred_ll_homref

    @gt_phred_ll_homref.setter
    def gt_phred_ll_homref(self, gt_phred_ll_homref):
        if gt_phred_ll_homref is not None:
            assert len(gt_phred_ll_homref) == len(self.family)
            self._gt_phred_ll_homref = gt_phred_ll_homref

    @property
    def gt_phred_ll_homalt(self):
        return self._gt_phred_ll_homalt

    @gt_phred_ll_homalt.setter
    def gt_phred_ll_homalt(self, gt_phred_ll_homalt):
        if gt_phred_ll_homalt is not None:
            assert len(gt_phred_ll_homalt) == len(self.family)
            self._gt_phred_ll_homalt = gt_phred_ll_homalt

    @property
    def gt_phred_ll_het(self):
        return self._gt_phred_ll_het

    @gt_phred_ll_het.setter
    def gt_phred_ll_het(self, gt_phred_ll_het):
        if gt_phred_ll_het is not None:
            assert len(gt_phred_ll_het) == len(self.family)
            self._gt_phred_ll_het = gt_phred_ll_het

    def __getattr__(self, gt):
        assert self._gt_types is not None
        def func(*args, **kwargs):
            if 'min_depth' in kwargs:
                assert self._gt_depths is not None
            debug = kwargs.pop('debug', False)
            flt = getattr(self.family, gt)(*args, **kwargs)
            if flt is False or flt is None:
                return False
            if gt == "comp_het_pair":
                return flt
            env = {s.sample_id: i for i, s in enumerate(self.family.subjects)}
            if debug:
                print(flt, file=sys.stderr)
            env['gt_types'] = self.gt_types
            env['gt_quals'] = self.gt_quals
            env['gt_depths'] = self.gt_depths
            env['gt_phred_ll_homref'] = self.gt_phred_ll_homref
            env['gt_phred_ll_het'] = self.gt_phred_ll_het
            env['gt_phred_ll_homalt'] = self.gt_phred_ll_homalt
            try:
                return eval(flt, env)
            except:
                print("attempted eval:", flt)
                raise
        return func

    def to_vcf(self, fh, var_dict=None, header=True, _POS=[100001]):
        if header:
            fh.write("##fileformat=VCFv4.1\n")
            fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
            fh.write("\t".join(s.name for s in self.subjects) + "\n")
        if var_dict is None:
            var_dict = {}

        for k in ("ID", "QUAL", "INFO"):
            if k not in var_dict:
                var_dict[k] = "."
        var_dict["FILTER"] = "PASS"
        var_dict["FORMAT"] = "GT"
        if not "CHROM" in var_dict:
            var_dict["CHROM"] = "1"
        if not "POS" in var_dict:
            var_dict["POS"] = _POS[0]
            _POS[0] += 1

        if not "REF" in var_dict:
            var_dict["REF"] = "A"
        if not "ALT" in var_dict:
            var_dict["ALT"] = "G"

        # convert from number back to repr
        x = ["0/0", "0/1", "./.", "1/1"]
        formats = [x[t] for t in self.gt_types]
        if self.gt_depths:
            var_dict["FORMAT"] += ":DP"
            for i, d in enumerate(self.gt_depths):
                formats[i] += (":%d" % d)

        fh.write("{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t".format(**var_dict))
        fh.write("\t".join(formats) + "\n")
