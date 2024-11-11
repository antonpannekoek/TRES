# to do
# min teken in mean anomaly
"""
Triple:      Triple evolution
             computes the evolution of a given triple
             given any initial conditions (M, m, l, A, a, E, e, i, G, g, O, o, T, z).
"""
import sys
import argparse
import numpy as np

from amuse.units import units
from amuse.support.console import set_printing_strategy
from amuse.community.seba import Seba

from tres.seculartriple import Seculartriple
from tres.triple_class import Triple
from tres.plotting import PlotDataContainer, plot_function
from tres.setup import make_particle_sets, setup_stellar_code
from tres.options import (
    REPORT_DEBUG,
    REPORT_TRIPLE_EVOLUTION,
    MAKE_PLOTS,
    REPORT_USER_WARNINGS,
)
from tres.interactions import (
    corotating_spin_angular_frequency_binary,
    lang_spin_angular_frequency,
    break_up_angular_frequency,
    criticial_angular_frequency_CHE,
)

# from interactions import *
# from tidal_friction_constant import *


def initialize_triple_class(
    stars,
    bins,
    correct_params,
    stellar_code,
    secular_code,
    relative_inclination=80.0 | units.deg,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.hdf",
    file_type="hdf5",
    dir_plots="",
):

    triple = Triple(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )
    triple.stellar_code.parameters.metallicity = metallicity

    return triple


# -----
# for running TRES.py from other routines
def tres_main(
    inner_primary_mass=1.3 | units.MSun,
    inner_secondary_mass=0.5 | units.MSun,
    outer_mass=0.5 | units.MSun,
    inner_semimajor_axis=1.0 | units.au,
    outer_semimajor_axis=100.0 | units.au,
    inner_eccentricity=0.1,
    outer_eccentricity=0.5,
    relative_inclination=80.0 | units.deg,
    inner_argument_of_pericenter=0.1,
    outer_argument_of_pericenter=0.5,
    inner_longitude_of_ascending_node=0.0,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.hdf",
    file_type="hdf5",
    dir_plots="",
    stellar_code=None,
    secular_code=None,
):

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    relative_inclination = relative_inclination
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)

    stars, bins, correct_params = make_particle_sets(
        inner_primary_mass,
        inner_secondary_mass,
        outer_mass,
        inner_semimajor_axis,
        outer_semimajor_axis,
        inner_eccentricity,
        outer_eccentricity,
        relative_inclination,
        inner_argument_of_pericenter,
        outer_argument_of_pericenter,
        inner_longitude_of_ascending_node,
    )

    clean_up_stellar_code = False
    clean_up_secular_code = False
    if stellar_code is None:
        stellar_code = Seba()
        # stellar_code = Seba(redirection='none')
        # stellar_code = Seba(redirection='file', redirect_file='output_SeBa_TRES.txt')
        clean_up_stellar_code = True

    stellar_code.parameters.metallicity = metallicity
    if secular_code is None:
        secular_code = Seculartriple()
        # secular_code = Seculartriple(redirection='none')
        # secular_code = Seculartriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')
        clean_up_secular_code = True

    triple_class_object = Triple(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )

    if triple_class_object.correct_params is False:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The parameters of the given triple are incorrect."
            )
        return triple_class_object  # no codes initialized yet
    elif (
        stop_at_semisecular_regime is True
        and triple_class_object.semisecular_regime_at_initialisation is True
    ):
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is in the semisecular regime at initialization."
            )
    elif triple_class_object.dynamical_instability_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is dynamically unstable at initialization."
            )
    elif triple_class_object.mass_transfer_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. There is mass transfer in the given triple at initialization."
            )
    elif stop_at_no_CHE is True and triple_class_object.CHE_at_initialisation is False:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. No chemically homogeneous evolution at initialization"
            )
    else:
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()

    stellar_code.particles.remove_particles(stars)
    triple_set = triple_class_object.triple.as_set()
    secular_code.triples.remove_particles(triple_set)
    del stars, bins, triple_set

    if clean_up_stellar_code:
        triple_class_object.stellar_code.stop()
        print("cleaning se")
    if clean_up_secular_code:
        triple_class_object.secular_code.stop()
        print("cleaning sec")

    return triple_class_object


def tres_main_developer(
    stars,
    bins,
    correct_params,
    stellar_code,
    secular_code,
    relative_inclination=80.0 | units.deg,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.hdf",
    file_type="hdf5",
    dir_plots="",
):

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    bins.eccentricity[0] = float(bins.eccentricity[0])
    bins.eccentricity[1] = float(bins.eccentricity[1])
    bins.argument_of_pericenter[0] = float(bins.argument_of_pericenter[0])
    bins.argument_of_pericenter[1] = float(bins.argument_of_pericenter[1])
    bins.longitude_of_ascending_node[0] = float(bins.longitude_of_ascending_node[0])
    bins.longitude_of_ascending_node[1] = float(bins.longitude_of_ascending_node[1])
    relative_inclination = relative_inclination

    triple_class_object = Triple(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )

    if triple_class_object.correct_params is False:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The parameters of the given triple are "
                "incorrect."
            )
        return triple_class_object  # no codes initialized yet
    if (
        stop_at_semisecular_regime is True
        and triple_class_object.semisecular_regime_at_initialisation is True
    ):
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is in the semisecular "
                "regime at initialization."
            )
    elif triple_class_object.dynamical_instability_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is dynamically unstable "
                "at initialization."
            )
    elif triple_class_object.mass_transfer_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. There is mass transfer in the given triple "
                "at initialization."
            )
    elif stop_at_no_CHE is True and triple_class_object.CHE_at_initialisation is False:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. No chemically homogeneous evolution at "
                "initialization"
            )
    else:
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()

    return triple_class_object


# -----


# -----
# for running triple.py from the commandline
def parse_arguments():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-M",
        "--M1",
        type=units.MSun,
        dest="inner_primary_mass",
        default=1.3 | units.MSun,
        help="inner primary mass",
    )
    parser.add_argument(
        "-m",
        "--M2",
        type=units.MSun,
        dest="inner_secondary_mass",
        default=0.5 | units.MSun,
        help="inner secondary mass",
    )
    parser.add_argument(
        "-l",
        "--M3",
        type=units.MSun,
        dest="outer_mass",
        default=0.5 | units.MSun,
        help="outer mass",
    )

    parser.add_argument(
        "-A",
        "--Ain",
        type=units.RSun,
        dest="inner_semimajor_axis",
        default=200.0 | units.RSun,
        help="inner semi major axis",
    )
    parser.add_argument(
        "-a",
        "--Aout",
        type=units.RSun,
        dest="outer_semimajor_axis",
        default=20000.0 | units.RSun,
        help="outer semi major axis",
    )
    parser.add_argument(
        "-E",
        "--Ein",
        dest="inner_eccentricity",
        type=float,
        default=0.1,
        help="inner eccentricity",
    )
    parser.add_argument(
        "-e",
        "--Eout",
        dest="outer_eccentricity",
        type=float,
        default=0.5,
        help="outer eccentricity",
    )
    parser.add_argument(
        "-i",
        "-I",
        dest="relative_inclination",
        type=units.rad,
        default=80.0 | units.deg,
        help="relative inclination",
    )
    parser.add_argument(
        "-G",
        "--Gin",
        dest="inner_argument_of_pericenter",
        type=units.rad,
        default=0.1 | units.rad,
        help="inner argument of pericenter",
    )
    parser.add_argument(
        "-g",
        "--Gout",
        dest="outer_argument_of_pericenter",
        type=units.rad,
        default=0.5 | units.rad,
        help="outer argument of pericenter",
    )
    parser.add_argument(
        "-O",
        "--Oin",
        dest="inner_longitude_of_ascending_node",
        type=units.rad,
        default=0.0 | units.rad,
        help="inner longitude of ascending node [rad]",
    )
    # outer longitude of ascending nodes = inner - pi
    # parser.add_argument(
    #     "-o",
    #     dest="outer_longitude_of_ascending_node",
    #     type=units.rad,
    #     default=0.0 | units.rad,
    #     help="outer longitude of ascending node"
    # )

    parser.add_argument(
        "-z",
        "-Z",
        dest="metallicity",
        type=float,
        default=0.02,
        help="metallicity",
    )
    parser.add_argument(
        "-t",
        "-T",
        type=units.Myr,
        dest="tend",
        default=5.0 | units.Myr,
        help="end time",
    )
    parser.add_argument(
        "--initial_time",
        type=units.Myr,
        dest="tinit",
        default=0.0 | units.Myr,
        help="initial time",
    )
    parser.add_argument(
        "-N",
        dest="number",
        type=int,
        default=0,
        help="number ID of system",
    )
    parser.add_argument(
        "-r",
        dest="maximum_radius_change_factor",
        type=float,
        default=0.01,
        help="maximum_radius_change_factor",
    )

    # parser.add_argument(
    #     "--tidal",
    #     dest="tidal_terms",
    #     action="store_false",
    #     default=True,
    #     help="tidal terms included"
    # )

    parser.add_argument(
        "--no_stop_at_mass_transfer",
        dest="stop_at_mass_transfer",
        action="store_false",
        default=True,
        help="stop at mass transfer",
    )
    parser.add_argument(
        "--no_stop_at_init_mass_transfer",
        dest="stop_at_init_mass_transfer",
        action="store_false",
        default=True,
        help="stop if initially mass transfer",
    )
    parser.add_argument(
        "--no_stop_at_outer_mass_transfer",
        dest="stop_at_outer_mass_transfer",
        action="store_false",
        default=True,
        help="stop at triple mass transfer",
    )

    #   if stop_at_mass_transfer is False, the following 4 stopping conditions
    #   can be used to further specify.
    #   if stop_at_mass_transfer is True, the following 4 are ignored.
    parser.add_argument(
        "--stop_at_stable_mass_transfer",
        dest="stop_at_stable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at stable mass transfer",
    )
    parser.add_argument(
        "--stop_at_eccentric_stable_mass_transfer",
        dest="stop_at_eccentric_stable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at eccentric stable mass transfer",
    )
    # unstable mass transfer leads to common-envelope evolution
    parser.add_argument(
        "--stop_at_unstable_mass_transfer",
        dest="stop_at_unstable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at unstable mass transfer",
    )
    parser.add_argument(
        "--stop_at_eccentric_unstable_mass_transfer",
        dest="stop_at_eccentric_unstable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at eccentric unstable mass transfer",
    )
    # 0  alpha-ce + alpha-dce
    # 1  gamma-ce + alpha-dce
    # 2  seba style; combination of gamma-ce, alpha-ce & alpha-dce
    parser.add_argument(
        "--CE",
        dest="which_common_envelope",
        type=int,
        default=2,
        help="which common envelope modeling",
    )

    parser.add_argument(
        "--stop_at_no_CHE",
        dest="stop_at_no_CHE",
        action="store_true",
        default=False,
        help="stop if no chemically homogeneous evolution",
    )
    parser.add_argument(
        "--include_CHE",
        dest="include_CHE",
        action="store_true",
        default=False,
        help="include chemically homogeneous evolution in the stellar evolution",
    )

    parser.add_argument(
        "--no_stop_at_merger",
        dest="stop_at_merger",
        action="store_false",
        default=True,
        help="stop at merger",
    )
    parser.add_argument(
        "--no_stop_at_disintegrated",
        dest="stop_at_disintegrated",
        action="store_false",
        default=True,
        help="stop at disintegrated",
    )
    parser.add_argument(
        "--no_stop_at_inner_collision",
        dest="stop_at_inner_collision",
        action="store_false",
        default=True,
        help="stop at collision in inner binary",
    )
    parser.add_argument(
        "--no_stop_at_outer_collision",
        dest="stop_at_outer_collision",
        action="store_false",
        default=True,
        help="stop at collision in outer binary",
    )
    parser.add_argument(
        "--no_stop_at_dynamical_instability",
        dest="stop_at_dynamical_instability",
        action="store_false",
        default=True,
        help="stop at dynamical instability",
    )
    parser.add_argument(
        "--stop_at_semisecular_regime",
        dest="stop_at_semisecular_regime",
        action="store_true",
        default=False,
        help="stop at semisecular regime",
    )

    parser.add_argument(
        "--stop_at_SN",
        dest="stop_at_SN",
        action="store_true",
        default=False,
        help="stop at supernova",
    )
    # 0  No kick
    # 1  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974
    # 2  Arzoumanian ea 2002, 568, 289
    # 3  Hansen & Phinney 1997, 291, 569
    # 4  Paczynski 1990, 348, 485
    # 5  Verbunt, Igoshev & Cator, 2017, 608, 57
    parser.add_argument(
        "--SN_kick_distr",
        dest="SN_kick_distr",
        type=int,
        default=5,
        help="which supernova kick distribution",
    )
    parser.add_argument(
        "--no_impulse_kick_for_black_holes",
        dest="impulse_kick_for_black_holes",
        action="store_false",
        default=True,
        help="do not rescale the BH SN kick by mass -> impulse kick",
    )
    parser.add_argument(
        "--no_fallback_kick_for_black_holes",
        dest="fallback_kick_for_black_holes",
        action="store_false",
        default=True,
        help="do not rescale the BH SN kick with fallback ",
    )

    parser.add_argument(
        "--stop_at_CPU_time",
        dest="stop_at_CPU_time",
        action="store_true",
        default=False,
        help="stop at CPU time",
    )
    parser.add_argument(
        "--max_CPU_time",
        dest="max_CPU_time",
        type=float,
        default=3600.0,
        help="max CPU time",
    )

    parser.add_argument(
        "-f",
        dest="file_name",
        type=str,
        default="TRES.hdf",  # "TRES.txt"
        help="file name",
    )
    parser.add_argument(
        "-F",
        dest="file_type",
        type=str,
        default="hdf5",  # "txt"
        help="file type",
    )
    parser.add_argument(
        "--dir_plots",
        dest="dir_plots",
        type=str,
        default="",
        help="directory for plots for debugging mode",
    )

    args = parser.parse_args()
    return args.__dict__


def main():
    arg = parse_arguments()

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    stars, bins, correct_params = make_particle_sets(
        arg["inner_primary_mass"],
        arg["inner_secondary_mass"],
        arg["outer_mass"],
        arg["inner_semimajor_axis"],
        arg["outer_semimajor_axis"],
        arg["inner_eccentricity"],
        arg["outer_eccentricity"],
        arg["relative_inclination"],
        arg["inner_argument_of_pericenter"],
        arg["outer_argument_of_pericenter"],
        arg["inner_longitude_of_ascending_node"],
    )

    stellar_code = Seba()
    # stellar_code = Seba(redirection='none')
    # stellar_code = Seba(redirection='file', redirect_file='output_SeBa_TRES.txt')
    stellar_code.parameters.metallicity = arg["metallicity"]
    secular_code = Seculartriple()
    # secular_code = Seculartriple(redirection='none')
    # secular_code = Seculartriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')

    triple_class_object = Triple(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        arg["relative_inclination"],
        arg["tend"],
        arg["tinit"],
        arg["number"],
        arg["maximum_radius_change_factor"],
        arg["stop_at_mass_transfer"],
        arg["stop_at_init_mass_transfer"],
        arg["stop_at_outer_mass_transfer"],
        arg["stop_at_stable_mass_transfer"],
        arg["stop_at_eccentric_stable_mass_transfer"],
        arg["stop_at_unstable_mass_transfer"],
        arg["stop_at_eccentric_unstable_mass_transfer"],
        arg["which_common_envelope"],
        arg["stop_at_no_CHE"],
        arg["include_CHE"],
        arg["stop_at_merger"],
        arg["stop_at_disintegrated"],
        arg["stop_at_inner_collision"],
        arg["stop_at_outer_collision"],
        arg["stop_at_dynamical_instability"],
        arg["stop_at_semisecular_regime"],
        arg["stop_at_SN"],
        arg["SN_kick_distr"],
        arg["impulse_kick_for_black_holes"],
        arg["fallback_kick_for_black_holes"],
        arg["stop_at_CPU_time"],
        arg["max_CPU_time"],
        arg["file_name"],
        arg["file_type"],
        arg["dir_plots"],
    )

    if triple_class_object.correct_params is False:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The parameters of the given triple are incorrect."
            )
        # no codes initialized yet
        sys.exit(
            "Choose a different system. The parameters of the given triple are incorrect."
        )
    elif (
        arg["stop_at_semisecular_regime"] is True
        and triple_class_object.semisecular_regime_at_initialisation is True
    ):
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is in the semisecular regime at initialization."
            )
    elif triple_class_object.dynamical_instability_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. The given triple is dynamically unstable at initialization."
            )
    elif triple_class_object.mass_transfer_at_initialisation is True:
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. There is mass transfer in the given triple at initialization."
            )
    elif (
        arg["stop_at_no_CHE"] is True
        and triple_class_object.CHE_at_initialisation is False
    ):
        if REPORT_USER_WARNINGS:
            print(
                "Choose a different system. No chemically homogeneous evolution at initialization"
            )
    else:
        triple_class_object.evolve_model(arg["tend"])
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, arg["dir_plots"])
            triple_class_object.print_stellar_system()

        if REPORT_TRIPLE_EVOLUTION:
            print("Simulation has finished succesfully")

    print("\nYou have used the TRES triple evolution code. Literature reference:")
    print("** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:")
    print('... "The evolution of hierarchical triple star-systems" ')

    triple_class_object.stellar_code.stop()
    triple_class_object.secular_code.stop()


if __name__ == "__main__":
    # main()
    tres_main()
