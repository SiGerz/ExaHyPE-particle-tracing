# This file is part of the ExaHyPE2 project. For conditions of distribution and
# use, please see the copyright notice at www.peano-framework.org
import argparse

import exahype2

"""
This is a particle tracing simulation for volcanic plume rise using the Euler equations.
"""


def generate_particles_in_plume(
    file_path: str = "InitTracers.dat",
    num_x: int = 8,
    num_y: int = 8,
    num_z: int = 8,
    x_min: float = 0.46,
    x_max: float = 0.54,
    y_min: float = 0.11,
    y_max: float = 0.29,
    z_min: float = 0.46,
    z_max: float = 0.54,
    round_decimals: int = 6,
):
    """
    Generates a .dat file with particles on a uniform grid inside a rectangular plume.
    """
    if num_x < 1 or num_y < 1:
        raise ValueError("num_x and num_y must be >= 1")
    if x_min > x_max:
        raise ValueError("x_min must be <= x_max")
    if y_min > y_max:
        raise ValueError("y_min must be <= y_max")

    if num_x == 1:
        x_values = [round((x_min + x_max) / 2.0, round_decimals)]
    else:
        x_step = (x_max - x_min) / (num_x - 1)
        x_values = [round(x_min + i * x_step, round_decimals) for i in range(num_x)]

    if num_y == 1:
        y_values = [round((y_min + y_max) / 2.0, round_decimals)]
    else:
        y_step = (y_max - y_min) / (num_y - 1)
        y_values = [round(y_min + j * y_step, round_decimals) for j in range(num_y)]

    if num_z == 1:
            z_values = [round((z_min + z_max) / 2.0, round_decimals)]
    else:
        z_step = (z_max - z_min) / (num_z - 1)
        z_values = [round(z_min + j * z_step, round_decimals) for j in range(num_z)]

    with open(file_path, "w") as file:
        for x in x_values:
            for y in y_values:
                for z in z_values:
                    file.write(f"{x:.{round_decimals}f} {y:.{round_decimals}f} {z:.{round_decimals}f}\n")

    total = num_x * num_y * num_z
    print(f"Generated {total} particles inside plume region and saved to {file_path}")


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
        help="show this help message and continue with project help",
    )
    parser.add_argument(
        "-interp",
        "--interpolation",
        choices=interpolation_methods.keys(),
        default="Constant",
    )
    parser.add_argument(
        "-integ", "--integration", choices=integration_schemes.keys(), default="RK4"
    )
    parser.add_argument("-et", "--end-time", type=float, default=3.0)
    parser.add_argument("-dofs", "--degrees-of-freedom", type=int, default=16)
    parser.add_argument("-l", "--grid-levels", type=int, default=2)
    args, _ = parser.parse_known_args()
    if args.help:
        parser.print_help()
        print()
    return args


args = parse_arguments()

project = exahype2.Project(
    namespace=["tests", "exahype2", "fv"],
    project_name=".",
    dimensions=3,
    size=[1.0, 1.0, 1.0],
    offset=[0.0, 0.0, 0.0],
    periodic_bc=[False, False, False],
    min_end_time=args.end_time,
    max_end_time=args.end_time,
    min_grid_levels=args.grid_levels,
    max_grid_levels=args.grid_levels,
)

#if project.dimensions != 2:
#    raise ValueError("The volcanic plume tracer test runs in 2D only.")

initial_conditions = """
  const double plumeCentreX = DomainOffset[0] + 0.5 * DomainSize[0];
  const double plumeCentreZ = DomainOffset[2] + 0.5 * DomainSize[2];
  const double plumeRadius  = 0.05 * std::fmin(DomainSize[0], DomainSize[2]);

  const double plumeBaseY = DomainOffset[1] + 0.1 * DomainSize[1];
  const double plumeTopY  = DomainOffset[1] + 0.3 * DomainSize[1];

  const double dx = x(0) - plumeCentreX;
  const double dz = x(2) - plumeCentreZ;

  const bool inPlumeColumn =
      dx * dx + dz * dz <= plumeRadius * plumeRadius;

  const bool inPlume =
      inPlumeColumn && x(1) >= plumeBaseY && x(1) <= plumeTopY;

  double localPressure = 0.0;

  if (inPlumeColumn) {
    if (x(1) >= plumeTopY) {
      localPressure =
          P0 + GRAVITY * (DomainOffset[1] + DomainSize[1] - x(1))
                   * RHO_AMBIENT;
    } else if (x(1) < plumeBaseY) {
      localPressure =
          P0
          + GRAVITY * (DomainOffset[1] + DomainSize[1] - plumeTopY)
                        * RHO_AMBIENT
          + GRAVITY * (plumeTopY - plumeBaseY) * RHO_HOT
          + GRAVITY * (plumeBaseY - x(1)) * RHO_AMBIENT;
    } else {
      localPressure =
          P0
          + GRAVITY * (DomainOffset[1] + DomainSize[1] - plumeTopY)
                        * RHO_AMBIENT
          + GRAVITY * (plumeTopY - x(1)) * RHO_HOT;
    }
  } else {
    localPressure =
        P0 + GRAVITY * (DomainOffset[1] + DomainSize[1] - x(1))
                 * RHO_AMBIENT;
  }

  const double localEnergy = localPressure / (GAMMA - 1.0);

  Q[Shortcuts::rho]  = inPlume ? RHO_HOT : RHO_AMBIENT;
  Q[Shortcuts::rhoU] = 0.0;
  Q[Shortcuts::rhoV] = 0.0;
  Q[Shortcuts::rhoW] = 0.0;
  Q[Shortcuts::rhoE] = localEnergy;
"""


boundary_conditions = """
  const double bottomY = DomainOffset[1];
  const double topY    = DomainOffset[1] + DomainSize[1];
  const double tol     = 0.001;

  if (x(1) >= topY - tol) {
    Qoutside[Shortcuts::rho]  = RHO_AMBIENT;
    Qoutside[Shortcuts::rhoU] = 0.0;
    Qoutside[Shortcuts::rhoV] = 0.0;
    Qoutside[Shortcuts::rhoW] = 0.0;
    Qoutside[Shortcuts::rhoE] = P0 / (GAMMA - 1.0);
  } else if (x(1) <= bottomY + tol) {
    Qoutside[Shortcuts::rho]  = Qinside[Shortcuts::rho];
    Qoutside[Shortcuts::rhoU] = Qinside[Shortcuts::rhoU];
    Qoutside[Shortcuts::rhoV] = -Qinside[Shortcuts::rhoV];
    Qoutside[Shortcuts::rhoW] = Qinside[Shortcuts::rhoW];
    Qoutside[Shortcuts::rhoE] = Qinside[Shortcuts::rhoE];
  } else {
    Qoutside[Shortcuts::rho]  = Qinside[Shortcuts::rho];
    Qoutside[Shortcuts::rhoU] = Qinside[Shortcuts::rhoU];
    Qoutside[Shortcuts::rhoV] = Qinside[Shortcuts::rhoV];
    Qoutside[Shortcuts::rhoW] = Qinside[Shortcuts::rhoW];
    Qoutside[Shortcuts::rhoE] = Qinside[Shortcuts::rhoE];
  }
"""


compute_primitive_variables = r"""
  const auto irho = 1.0 / Q[Shortcuts::rho];

  const auto u0 = Q[Shortcuts::rhoU] * irho;
  const auto u1 = Q[Shortcuts::rhoV] * irho;
  const auto u2 = Q[Shortcuts::rhoW] * irho;

  const auto uSq =
      tarch::la::productRoundedOnce(u0, u0)
    + tarch::la::productRoundedOnce(u1, u1)
    + tarch::la::productRoundedOnce(u2, u2);

  const auto u_n =
      (normal == 0) ? u0 : ((normal == 1) ? u1 : u2);

  const auto internalE =
      Q[Shortcuts::rhoE] - 0.5 * Q[Shortcuts::rho] * uSq;
  const auto p = (GAMMA - 1.0) * internalE;
"""


max_eigenvalue = r"""
  {compute_primitive_variables}
  const auto speedOfSound = std::sqrt(GAMMA * p * irho);
  auto result = std::fmax(0.0, std::fabs(u_n - speedOfSound));
  result = std::fmax(result, std::fabs(u_n + speedOfSound));
  return result;
""".format(compute_primitive_variables=compute_primitive_variables)

flux = r"""
  {compute_primitive_variables}

  F[Shortcuts::rho] =
      (normal == 0) ? Q[Shortcuts::rhoU]
    : (normal == 1) ? Q[Shortcuts::rhoV]
                    : Q[Shortcuts::rhoW];

  F[Shortcuts::rhoU] = Q[Shortcuts::rhoU] * u_n;
  F[Shortcuts::rhoV] = Q[Shortcuts::rhoV] * u_n;
  F[Shortcuts::rhoW] = Q[Shortcuts::rhoW] * u_n;

  if (normal == 0) F[Shortcuts::rhoU] += p;
  if (normal == 1) F[Shortcuts::rhoV] += p;
  if (normal == 2) F[Shortcuts::rhoW] += p;

  F[Shortcuts::rhoE] = (Q[Shortcuts::rhoE] + p) * u_n;
""".format(compute_primitive_variables=compute_primitive_variables)


source_term = """
  S[Shortcuts::rho]  = 0.0;
  S[Shortcuts::rhoU] = 0.0;
  S[Shortcuts::rhoV] = -GRAVITY * Q[Shortcuts::rho];
  S[Shortcuts::rhoW] = 0.0;
  S[Shortcuts::rhoE] = -GRAVITY * Q[Shortcuts::rhoV];
"""


fv_solver = exahype2.fv.Solver(
    name="FVSolver",
    patch_size=args.degrees_of_freedom,
    unknowns={"rho": 1, "rhoU": 1, "rhoV": 1, "rhoW": 1, "rhoE": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.4),
    scheme=exahype2.fv.Godunov(),
)

fv_solver.set_pde_terms(
    initial_condition=initial_conditions,
    boundary_conditions=boundary_conditions,
    max_eigenvalue=max_eigenvalue,
    flux=flux,
    source=source_term,
)

project.add_solver(fv_solver)

tracer_attributes = {
    "unknowns": {
        "rho": 1,
        "fluid_u": 1,
        "fluid_v": 1,
        "fluid_w": 1,
        "u": 1,
        "v": 1,
        "w": 1,
        "a_x": 1,
        "a_y": 1,
        "a_z": 1,
        "d": 1,
    }
}

tracer_particles = project.add_tracer(
    name="Tracers", particle_attributes=tracer_attributes
)

generate_particles_in_plume(
    file_path="InitTracers.dat",
    num_x=8,
    num_y=8,
    num_z=8,
    x_min=0.46,
    x_max=0.54,
    y_min=0.11,
    y_max=0.29,
    z_min=0.46,
    z_max=0.54,
)
tracer_particles.insert_from_file("InitTracers.dat", scale_factor=1.0)

fv_solver.trace(
    tracer_particles,
    forces=exahype2.particles.Forces(
        integration_scheme=integration_schemes[args.integration],
        interpolation_scheme=interpolation_methods[args.interpolation],
    ),
)

project.output.makefile.add_compile_definition("GAMMA=1.4")
project.output.makefile.add_compile_definition("GRAVITY=2.0")
project.output.makefile.add_compile_definition("P0=1.0")
project.output.makefile.add_compile_definition("RHO_AMBIENT=1.0")
project.output.makefile.add_compile_definition("RHO_HOT=0.05")
project.generate_and_build()
