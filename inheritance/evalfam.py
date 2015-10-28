from __future__ import print_function

import os
import sys
import tempfile
import atexit

import itertools as it

def tmp(pedstr, suf=".ped"):
    t = tempfile.mktemp(suffix=suf)
    atexit.register(os.unlink, t)
    with open(t, "w") as fh:
        for line in pedstr.split("\n"):
            if not line.strip(): continue
            print(line.strip(), file=fh)
    return t


class EvalFamily(object):

    __slots__ = ('ped', 'family', 'gt_types', '_gt_types', 'gt_depths',
                 '_gt_depths', 'strict', 'subjects')

    def __init__(self, family, fam_id=None, gt_types=None, gt_depths=None):
        # can send in a family.
        if isinstance(family, dict):
            if len(family) != 1:
                raise Exception("only single families supported in EvalFamily")
            family = family.values()[0]
        self.family = family
        for s in self.family.subjects:
            if s.sample_id[0].isdigit(): s.sample_id = "s" + s.sample_id
        self.subjects = self.family.subjects

        self._gt_types = None
        self.gt_types = gt_types
        self._gt_depths = None
        self.gt_depths = gt_depths

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
    def gt_depths(self):
        return self._gt_depths

    @gt_depths.setter
    def gt_depths(self, gt_depths):
        if gt_depths is not None:
            assert len(gt_depths) == len(self.family)
            self._gt_depths = gt_depths


    def __getattr__(self, gt):
        assert self._gt_types
        def func(*args, **kwargs):
            if 'min_depth' in kwargs:
                assert self._gt_depths is not None
            debug = kwargs.pop('debug', False)
            flt = getattr(self.family, gt)(*args, **kwargs)
            if flt is False:
                return False
            if gt == "comp_het_pair":
                return flt
            env = {s.sample_id: i for i, s in enumerate(self.family.subjects)}
            if debug:
                print(flt, file=sys.stderr)
            env['gt_types'] = self.gt_types
            env['gt_depths'] = self.gt_depths
            try:
                return eval(flt, env)
            except:
                print(flt)
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
