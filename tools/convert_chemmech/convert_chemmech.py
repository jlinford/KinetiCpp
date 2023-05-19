#!/usr/bin/env python3
#
# Copyright 2023, John Linford <john@redhpc.com>
# SPDX-License-Identifier: Apache-2.0
#

import sys

snip = """
<1> NO2 = NO + O3P # 1.0/<NO2_06>;
<2> O3P + O2 + M = O3 # 5.68e-34^-2.60;
<3> O3P + O3 =  # 8.00e-12@2060;
<4> O3P + NO = NO2 # 9.00e-32^-1.50&3.00e-11&0.60&1.0;
<5> O3P + NO2 = NO # 5.50e-12@-188;
"""


def build_parser():
    import pyparsing as pp

    # Don't typecast number strings
    sci_num = pp.common.sci_real.copy().set_parse_action(None)
    real_num = pp.common.real.copy().set_parse_action(None)
    int_num = pp.common.integer.copy().set_parse_action(None)

    label = pp.Word(pp.alphanums+"<>").set_name("label")
    species = pp.Word(pp.identchars, pp.identbodychars+":").set_name("species")
    coef = (sci_num | real_num | int_num).set_name("coefficient")
    term = ((coef + "*") + species) | (species + ("*" + coef)) | species 
    expr = term + pp.ZeroOrMore((pp.Word("+") | "-") + term)
    
    rki_begin = "#" | ("%" + pp.one_of("1234H") + "#")
    rki_expr = pp.Word("".join([c for c in pp.printables if c not in ";"]))
    rki = rki_begin + rki_expr + ";"
    
    reaction = pp.Opt(label) + expr + pp.Opt("=" + expr) + rki

    parsed = reaction.parse_string("<1> N:O2 = 3e6*NO + O3P # 1.0/<NO2_06>;")
    print(parsed)


def main(argv):
    parser = build_parser()


if __name__ == "__main__":
    main(sys.argv)
