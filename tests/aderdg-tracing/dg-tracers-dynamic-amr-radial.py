# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import argparse

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

def parse_arguments():
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "-h",
        "--help",
        action="store_true",
        help="Show this help message and continue with project help.",
    )
    parser.add_argument("-et", "--end-time", type=float, default=0.4)
    parser.add_argument("-dofs", "--degrees-of-freedom", type=int, default=5)
    parser.add_argument("-l", "--grid-levels", type=int, default=None)
    parser.add_argument("-lmin", "--min-grid-levels", type=int, default=2)
    parser.add_argument("-lmax", "--max-grid-levels", type=int, default=4)
    args, _ = parser.parse_known_args()
    if args.grid_levels is not None:
        args.min_grid_levels = args.grid_levels
        args.max_grid_levels = args.grid_levels
    if args.help:
        parser.print_help()
        print()
    return args


args = parse_arguments()

project = exahype2.Project(
    namespace=["tests", "exahype2", "rkdg"],
    project_name="solution",
    dimensions=2,
    size=[2.0, 2.0, 2.0],
    offset=[-1.0, -1.0, -1.0],
    periodic_bc=[False, False, False],
    min_end_time=args.end_time,
    max_end_time=args.end_time,
    min_grid_levels=args.min_grid_levels,
    max_grid_levels=args.max_grid_levels,
)

initial_conditions = """
  const tarch::la::Vector<DIMENSIONS, double> centre(0.0);
  constexpr double BumpRadius = 0.3;
  const auto scaledRadiusSquared =
    tarch::la::norm2Squared(x - centre) / (BumpRadius * BumpRadius);
  Q[Shortcuts::rho] = 0.0;
  if (scaledRadiusSquared < 1.0) {
    const auto exponent = 1.0 - 1.0 / (1.0 - scaledRadiusSquared);
    const auto minimumNormalExponent = std::log(std::numeric_limits<double>::min());
    if (exponent >= minimumNormalExponent) {
      Q[Shortcuts::rho] = std::exp(exponent);
    }
  }
"""

boundary_conditions = """
  for (int i = 0; i < NumberOfUnknowns + NumberOfAuxiliaryVariables; i++) {
    Qoutside[i] = Qinside[i];
  }
"""

flux = """
  const tarch::la::Vector<DIMENSIONS, double> centre(0.0);
  const tarch::la::Vector<DIMENSIONS, double> radialOffset = x - centre;
  const double radius = tarch::la::norm2(radialOffset);
  F[Shortcuts::rho] = (radius > 0.0 ? radialOffset(normal) / radius : 0.0) * Q[Shortcuts::rho];
"""

# The unit offset bounds the relative detail where the density is zero.
refinement_criterion = """
  return ::exahype2::refineOnDetail(
    Q,
    Quadrature<double>::getQuadraturePoints1D(),
    NumberOfUnknowns + NumberOfAuxiliaryVariables,
    1.0e-2,
    ::exahype2::DefaultCoarsenSafetyFactor,
    [](const double* const values) { return 1.0 + values[Shortcuts::rho]; }
  );
"""

solver = exahype2.dg.ader.Solver(
    name="ADERDGSolver",
    rk_order=2,
    dg_order=args.degrees_of_freedom - 1,
    unknowns={"rho": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.9),
)

solver.set_pde_terms(
    initial_condition=initial_conditions,
    boundary_conditions=boundary_conditions,
    refinement_criterion=refinement_criterion,
    max_eigenvalue="return 1.0;",
    flux=flux,
)

project.add_solver(solver)


tracer_attributes = {
    "unknowns": {
        "rho": 1,
        "x": 1,
        "y":1,
        "z":1,
    }
}

tracer_particles = project.add_tracer(
    name="Tracers", particle_attributes=tracer_attributes
)

generate_grid_particles(
    "InitTracers.dat", num_rows=7, num_cols=8, x_width=2.0, y_width=2.0, centre_x=0.0, centre_y=0.0
)
tracer_particles.insert_from_file(
    "InitTracers.dat", scale_factor=1.0, save_initial_position=True
)

solver.trace(
    tracer_particles,
    forces=exahype2.particles.Forces(
        integration_scheme=1,
        interpolation_scheme=2,
        force_method=0,
        # The interpolation overwrites rho, u and v of each tracer with the fluid state.
        fields_to_interpolate="0",
    ),
    initialisation=False,
)


project.generate_and_build()
