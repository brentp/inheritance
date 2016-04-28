inheritance models for mendelian diseases
-----------------------------------------

[![Build Status](https://travis-ci.org/brentp/inheritance.svg?branch=master)](https://travis-ci.org/brentp/inheritance)

This module is a general-purpose framework for evaluating if a family exihibits, for example, and autosomal dominant pattern.
The logic for this was tuned in [gemini](https://github.com/arq5x/gemini) but we make it available here as a more general purpose library to encourage:
1. community driven improvements
2. use outside of gemini
3. more comprehensive testing

Finding variants that match autosomal dominance in a trio, for example is very simple to find,
however, after considering multiple generations, arbitrary family sizes, depth cutoffs, and unknown
genotypes and phenotypes to support to real-world datasets it becomes tedious and error-prone.

Supported inheritance tests
===========================

+ autosomal dominant
+ autosomal recessive
+ de novo
+ X-linked dominant, recessive, and de novo
+ compound heterozygote
+ mendelian violation

Usage
=====

For now, the use is via api only. Users can look at the tests to see how to use. Most functions have a signature like:

```Python
 auto_dom(self, min_depth=0, gt_ll=False, strict=True, only_affected=True)
```

where the arguments enforce a minimum depth, a maximum genotype likelihood, strictness (mostly related to parent-offspring requirements)
and wether to allow unaffecteds to have the variant (or be homozygous alt).

ToDo
====
1. add support for X-linked soon.
2. add a simple example of running on a VCF+PED
3. code coverage
4. code documentation

Testing
=======

Tests can be run as:

```
nosetests --with-coverage -x --with-doctest --cover-package inheritance
```

Overview
========

the generic code is in `inheritance/inheritance.py` and a specific implementation that we use in gemini is in `inheritance/pyeval.py`.
To make this available for a new resource, for example to `bcftools` we would look at the implementation of Filter in `inheritance/pyval.py`
and make the generated strings match those expected by `bcftools`.
