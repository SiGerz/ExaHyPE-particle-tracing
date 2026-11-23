"""
This file is can be used to play around and find out how the particle framework works.
"""

import argparse
import math
from pathlib import Path

import exahype2


def generate_grid_particles(
    file_path,
    num_rows=5,
    num_cols=4,
    x_width=2.0,
    y_width=1.0,
    centre_x=0.0,
    centre_y=0.0,
):
    """
    Generates a .dat file with particle positions on rays around (centre_x, centre_y).
    No ray lies on an axis. A tracer on an axis of a Hilbert grid reads only the
    cell below its face.
    """
    with open(file_path, "w") as file:
        for i in range(num_rows):
            for j in range(num_cols):
                x = centre_x - x_width/2 + (i * x_width/(num_rows - 1) if num_rows > 1 else x_width/2)
                y = centre_y - y_width/2 + (j * y_width/(num_cols - 1) if num_cols > 1 else y_width/2)
                file.write(f"{x:.6f} {y:.6f}\n")



interpolation_methods = {
    "None": 0,
    "Constant": 1,
    "Linear": 2,
}
integration_schemes = {
    "Static": 0,
    "Euler": 1,
    "SemiEuler": 2,
    "Verlet": 3,
    "Midpoint": 4,
    "RK3": 5,
    "RK4": 6,
}


def parse_arguments():
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "-h",
        "--help",
        action="store_true",
        help="Show this help message and continue with project help.",
    )
    parser.add_argument(
        "-interp",
        "--interpolation",
        choices=interpolation_methods.keys(),
        default="Linear",
    )
    parser.add_argument(
        "-integ", "--integration", choices=integration_schemes.keys(), default="RK4"
    )
    parser.add_argument("-et", "--end-time", type=float, default=0.8)
    parser.add_argument("-dofs", "--degrees-of-freedom", type=int, default=16)
    parser.add_argument("-l", "--grid-levels", type=int, default=3)
    args, _ = parser.parse_known_args()
    if args.help:
        parser.print_help()
        print()
    return args


args = parse_arguments()

project = exahype2.Project(
    namespace=["tests", "exahype2", "fv"],
    project_name="solution",
    dimensions=2,
    size=[4.0, 2.0, 2.0],
    offset=[-2.0, -1.0, -1.0],
    periodic_bc=[True, False, False],
    min_end_time=args.end_time,
    max_end_time=args.end_time,
    min_grid_levels=args.grid_levels,
    max_grid_levels=args.grid_levels,
)

initial_conditions = """
  constexpr double Ux = 1.0;
  constexpr double CentreY = 0.0;

  const double dy = x(1) - CentreY;
  const double ux = Ux * (1.0 - std::abs(dy));

  Q[Shortcuts::rho] = 1.0;
  Q[Shortcuts::u]   = ux;
  Q[Shortcuts::v]   = 0.0;
  Q[Shortcuts::e]   = 1.0;
"""



boundary_conditions = """
  for (int i = 0; i < NumberOfUnknowns; i++) {
    Qoutside[i] = Qinside[i];
  }
"""


flux = """
  constexpr double Ux = 1.0;
  constexpr double CentreY = 0.0;

  const double dy = x(1) - CentreY;
  const double speed = Ux * (1.0 - std::abs(dy));

  for (int i = 0; i < NumberOfUnknowns; i++) {
    if (normal == 0) {
      F[i] = speed * Q[i];  // x-Flux
    } else {
      F[i] = 0.0;           // kein y-Flux
    }
  }
"""


max_eigenvalue = """
  return 1.0;
"""

rkdg_solver = exahype2.dg.ader.Solver(
    name="ADERDGSolver",
    order = 4,
    unknowns={"rho": 1, "u": 1, "v": 1, "e": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.5),
    #scheme=exahype2.fv.Godunov(),
)

rkdg_solver.set_pde_terms(
    initial_condition=initial_conditions,
    boundary_conditions=boundary_conditions,
    max_eigenvalue=max_eigenvalue,
    flux=flux,
)

project.add_solver(rkdg_solver)

tracer_attributes = {
    "unknowns": {
        "rho": 1,
        "init_x": 1,
        "init_y": 1,
        "u": 1,
        "v": 1,
        "a_x": 1,
        "a_y": 1,
        "mass": 1,
    }
}

tracer_particles = project.add_tracer(
    name="Tracers", particle_attributes=tracer_attributes
)

generate_grid_particles(
    "InitTracers.dat", num_rows=5, num_cols=4, x_width=2.0, y_width=1.0, centre_x=0.0, centre_y=0.0
)
tracer_particles.insert_from_file(
    "InitTracers.dat", scale_factor=1.0, save_initial_position=True
)

rkdg_solver.trace(
    tracer_particles,
    forces=exahype2.particles.Forces(
        integration_scheme=integration_schemes[args.integration],
        interpolation_scheme=interpolation_methods[args.interpolation],
        force_method=0,
        # The interpolation overwrites rho, u and v of each tracer with the fluid state.
        fields_to_interpolate="0,3,4",
    ),
    initialisation=False,
)

project.generate_and_build()
