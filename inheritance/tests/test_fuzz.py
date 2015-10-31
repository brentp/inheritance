from inheritance import Sample, Family, EvalFamily
import random

def make_fam(n_affecteds, n_unaffecteds, n_unknowns, id="xxx"):

    samples = []
    for i in range(n_affecteds):
        samples.append(Sample('affected_%d' % i, affected=True))
    for i in range(n_unaffecteds):
        samples.append(Sample('unaffected_%d' % i, affected=False))
    for i in range(n_unknowns):
        samples.append(Sample('unknown_%d' % i, affected=None))

    for i in range((n_affecteds + n_affecteds + n_unknowns)/ 4):

        sample = random.choice(samples)
        if random.random() < 0.85:
            sample.dad = random.choice([s for s in samples if not s == sample])

        if random.random() < 0.85:
            sample.mom = random.choice([s for s in samples if not s == sample])

    fam = EvalFamily(Family(samples, 'fam_%s' % id))
    fam.gt_types = [random.randrange(0, 4) for _ in range(len(samples))]
    fam.gt_depths = [random.randrange(0, 100) for _ in range(len(samples))]
    return fam

def test_fuzz():

    for i in range(2000):
        n_affecteds, n_unaffecteds, n_unknowns = random.randint(0, 10), random.randint(0, 10), random.randint(0, 5)
        if n_affecteds + n_unaffecteds + n_unknowns == 0: continue

        f = make_fam(n_affecteds, n_unaffecteds, n_unknowns, str(i))

        for min_depth in (0, 10, 100):
            for only_affected in (True, False):
                for strict in (True, False):

                    f.auto_rec(min_depth=min_depth, only_affected=only_affected, strict=strict)
                    f.auto_dom(min_depth=min_depth, only_affected=only_affected, strict=strict)
                    f.de_novo(min_depth=min_depth, only_affected=only_affected, strict=strict)
                f.comp_het(min_depth=min_depth, only_affected=only_affected)
                f.mendel_plausible_denovo(min_depth=min_depth, gt_ll=False, only_affected=only_affected)
                f.mendel_implausible_denovo(min_depth=min_depth, gt_ll=False, only_affected=only_affected)
                f.mendel_uniparental_disomy(min_depth=min_depth, gt_ll=False, only_affected=only_affected)
                f.mendel_LOH(min_depth=min_depth, gt_ll=False, only_affected=only_affected)
