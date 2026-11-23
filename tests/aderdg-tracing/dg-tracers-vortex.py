import exahype2
import math


def generate_radial_particles(
    file_path,
    num_rays=8,
    particles_per_ray=1,
    r_min=0.5,
    r_max=1.0,
    centre_x=0.0,
    centre_y=0.0,
):
    """
    Generates a .dat file with particle positions on rays around (centre_x, centre_y).
    No ray lies on an axis. A tracer on an axis of a Hilbert grid reads only the
    cell below its face.
    """
    with open(file_path, "w") as file:
        for i in range(num_rays):
            theta = 2 * math.pi * (i + 0.5) / num_rays
            for j in range(particles_per_ray):
                r = r_min + (j * (r_max - r_min) / (particles_per_ray - 1) if particles_per_ray > 1 else 0)
                x = centre_x + r * math.cos(theta)
                y = centre_y + r * math.sin(theta)
                file.write(f"{x:.6f} {y:.6f}\n")

    total = num_rays * particles_per_ray
    print(
        f"Generated {total} radial particles in a 3x3 domain centred at ({centre_x},{centre_y}) and saved to {file_path}"
    )

project = exahype2.Project(
    namespace=["tests", "exahype2", "dg"],
    project_name="solution",
    dimensions=2,
    size=[3.0, 3.0, 3.0],
    offset=[-1.5, -1.5, -1.5],
    periodic_bc=[False, False, False],
    min_end_time=3.0,
    max_end_time=3.0,
    min_grid_levels=3,
    max_grid_levels=3,
)

initial_conditions = """
  constexpr double CentreX = 0.0;
  constexpr double CentreY = 0.0;
  constexpr double Gamma   = 1.0;  // rotation speed, sign determines direction of rotation
  constexpr double Epsilon = 1e-6;

  const double dx = x(0) - CentreX;
  const double dy = x(1) - CentreY;

  const double radiusSquared =
      tarch::la::productRoundedOnce(dx, dx)
    + tarch::la::productRoundedOnce(dy, dy)
    + Epsilon;

  Q[Shortcuts::rho] = 1.0;
  Q[Shortcuts::u] = -(Gamma / (2.0 * M_PI)) * dy / radiusSquared;
  Q[Shortcuts::v] =  (Gamma / (2.0 * M_PI)) * dx / radiusSquared;
  Q[Shortcuts::e] = 1.0;
"""

boundary_conditions = """
  Qoutside[0] = Qinside[0];
  Qoutside[1] = -Qinside[1];
  Qoutside[2] = -Qinside[2];
  Qoutside[3] = Qinside[3];
"""

flux = """
  for (int i = 0; i < NumberOfUnknowns; i++) {
    F[i] = 0.0;
   // if (normal == 0) {
   //     F[i] = Q[i + 1];
   // } else {
   //     F[i] = -Q[i-1]);
   // }
  }
"""

max_eigenvalue = """
  return 1.0;
"""

dg_solver = exahype2.dg.ader.Solver(
    name="ADERDGSolver",
    order = 4,
    unknowns={"rho": 1, "u": 1, "v": 1, "e": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.4),
)

solver = dg_solver

solver.set_pde_terms(
    initial_condition=initial_conditions,
    boundary_conditions=boundary_conditions,
    max_eigenvalue=max_eigenvalue,
    flux=flux,
)

project.add_solver(solver)

tracer_attributes = {
    "unknowns": {
        "rho": 1,
        "fluid_u": 1,
        "fluid_v": 1,
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

generate_radial_particles(
    "InitTracers.dat", num_rays=16, particles_per_ray=5, r_min=0.05, r_max=1.0
)
tracer_particles.insert_from_file(
    "InitTracers.dat", scale_factor=1.0, save_initial_position=True
)


solver.trace(
    tracer_particles,
    forces=exahype2.particles.Forces(
        integration_scheme=1,  # Explicit Euler
        interpolation_scheme=2,  # Piecewise linear
        force_method=0,
        fields_to_interpolate="0,3,4",
    ),
)

project.generate_and_build()