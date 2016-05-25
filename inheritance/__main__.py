import argparse
import re
import itertools as it
from collections import defaultdict

from .pyeval import Family
from .evalfam import EvalFamily
from geneimpacts import Effect

def main(args):

    p = argparse.ArgumentParser()
    p.add_argument("--inheritance_model", choices=("comp_het", "auto_dom", "auto_rec", "de_novo"))
    p.add_argument("--min-depth", type=int, default=5)
    p.add_argument("--min-gq", type=int, default=5)
    p.add_argument("--min-kindreds", type=int, default=0)
    p.add_argument("--min-severity", default=None, choices=("MED", "HIGH"))
    p.add_argument("ped")
    p.add_argument("vcf")

    a = p.parse_args(args)

    severity = a.min_severity
    if severity == "MED":
        # if they said MED, it means MED or HIGH(er).
        severity = ("MED", "HIGH")
    elif severity == "HIGH":
        severity = ("HIGH", )

    run(a.inheritance_model, a.ped, a.vcf, a.min_depth, a.min_gq, a.min_kindreds, severity)


def run(inheritance_model, ped, vcf, min_depth, min_gq, min_kindreds, severity):
    from cyvcf2 import VCF, Writer
    vcf = VCF(vcf, samples="-")

    annos = {}
    if "ANN" in vcf:
        desc = vcf["ANN"]["Description"]
        parts = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]
        annos["ANN"] = desc
    if "EFF" in vcf:
        desc = vcf["EFF"]["Description"]
        parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        annos["EFF"] = parts
    if "CSQ" in vcf:
        desc = vcf["CSQ"]["Description"]
        parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        annos["CSQ"] = parts

    vcf.update(id="inheritance", type="String", number="1", description="inheritance stuffs")
    out = Writer("-", vcf)

    vcf_order = dict((n, i) for i, n in (enumerate(vcf.samples)))
    fams = Family.from_ped(ped, order=vcf_order)
    for fam_id in fams:
        fams[fam_id] = (EvalFamily(fams[fam_id]), [s._i for s in fams[fam_id].subjects])

    def get_gene(variant):
        for anno in annos:
            consequences = variant.INFO[anno].split(",")
            effs = (Effect.new(anno, c, annos[anno]) for c in consequences)
            # limit to requested severity
            if severity is not None:
                effs = [e for e in effs if e.impact_severity in severity]
            effs = sorted(effs, reverse=True)
            for eff in effs:
                if eff.gene:
                    return eff.gene

    # TODO: more flexible groupby
    for gene, variants in it.groupby(vcf, get_gene):

        matching_fams = defaultdict(list)
        saved_vars = []
        uniq_fams = []

        for i, variant in enumerate(variants):
            saved_vars.append(variant)

            for family_id, (fam, idxs) in fams.items():
                fam.gt_types = variant.gt_types[idxs]
                fam.gt_depths = variant.gt_depths[idxs]
                fam.gt_quals = variant.gt_quals[idxs]
                # this dispatches to fam.auto_rec/auto_dom/de_novo/, etc. by the string
                # in inheritance model
                res = getattr(fam, inheritance_model)(min_depth=min_depth, min_gq=min_gq)

                # matched the inheritance model.
                if res: # can add custom logic here, e.g. and v.call_rate > 0.9:
                    matching_fams[i].append(family_id)
                    uniq_fams.append(family_id)

        if 0 < len(set(uniq_fams)) >= min_kindreds:

            if inheritance_model == 'comp_het':
                # TODO: idxs = matching_fams.keys()
                # run idxs[1:] vs idxs[:-1] for variants
                pass
            for i, family_ids in sorted(matching_fams.items()):
                variant = saved_vars[i]
                variant.INFO["inheritance"] = "%s:%s" % (gene, ",".join(set(family_ids)))

                out.write_record(variant)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
