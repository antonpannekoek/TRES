lib_inner_primary_mass_distr = {
    0: "Kroupa",  # default
    1: "Scalo",
    2: "Miller & Scalo",
    3: "Salpeter",
    4: "Logarithmically flat",
    5: "Eggleton",
    6: "Kroupa for massive stars M>0.5 powerlaw with exp=-2.3",
}

# draws from mass distribution instead of mass ratio distribution, appropriate
# for planets
lib_inner_mass_ratio_distr = {
    0: "Uniform distribution",  # default
    1: "Kroupa IMF",  # draws from mass distribution instead of mass ratio distribution
    2: "Galicher et al. 2016 powerlaw (M^-1.31)",
}

# draws from mass distribution instead of mass ratio distribution, appropriate
# for planets
lib_outer_mass_ratio_distr = {
    0: "Uniform distribution",  # default
    1: "Kroupa IMF",  # draws from mass distribution instead of mass ratio distribution
    2: "Galicher et al. 2016 powerlaw (M^-1.31)",
}

# appropriate for planets
lib_inner_semi_distr = {
    0: "Log Uniform distribution",  # default
    1: "Constant semi-major axis",
    2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
    3: "Lognormal mu = 10^3.5d, sigma = 2.3",
    4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
    5: "Sana et al. 2012",
    6: "flat distribution",
    7: "Galicher et al. 2016 powerlaw (a^-0.61)",
}

# appropriate for planets
lib_outer_semi_distr = {
    0: "Log Uniform distribution",  # default
    1: "Constant semi-major axis",
    2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
    3: "Lognormal mu = 10^3.5d, sigma = 2.3",
    4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
    5: "Sana et al. 2012",
    6: "flat distribution",
    7: "Galicher et al. 2016 powerlaw (a^-0.61)",
}

# appropriate for planets
lib_inner_ecc_distr = {
    0: "Thermal",  # default
    1: "Constant eccentricity",
    2: "Sana et al. 2012 e^-0.45",  # -> close binaries
    3: "Flat distribution",
    4: "Powerlaw e^0.5",
    5: "Bowler et al. 2020 Beta distribution",
}

# appropriate for planets
lib_outer_ecc_distr = {
    0: "Thermal",  # default
    1: "Constant eccentricity",
    2: "Sana et al. 2012 e^-0.45",  # -> close binaries
    3: "Flat distribution",
    4: "Powerlaw e^0.5",
    5: "Bowler et al. 2020 Beta distribution",
}

lib_incl_distr = {
    0: "Circular uniform distribution",  # default
    1: "Constant inclination",
}

lib_inner_aop_distr = {
    0: "Uniform distribution",  # default
    1: "Constant argument of pericenter",
}

lib_outer_aop_distr = {
    0: "Uniform distribution",  # default
    1: "Constant argument of pericenter",
}

# default
lib_inner_loan_distr = {
    0: "Circular niform distribution",
    1: "Constant longitude of ascending nodes",
}

# default
lib_SN_kick_distr = {
    0: "No kick",
    1: "Hobbs",  # Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974
    2: "Arzoumanian",  # Arzoumanian ea 2002, 568, 289
    3: "Hansen",  # Hansen & Phinney 1997, 291, 569
    4: "Paczynski",  # Paczynski 1990, 348, 485
    5: "Verbunt",  # Verbunt, Igoshev & Cator, 2017, 608, 57
}

lib_CE = {
    0: "alpha-ce + alpha-dce",
    1: "gamma-ce + alpha-dce",
    2: "seba style; combination of gamma-ce, alpha-ce & alpha-dce",
}
