import exahype2


project = exahype2.Project(
    namespace=["tests", "exahype2", "dg"],
    project_name="solution",
    dimensions=2,
    size=[1.0, 1.0, 1.0],
    offset=[0.0, 0.0, 0.0],
    periodic_bc=[True, False, False],
    min_end_time=10.0,
    max_end_time=10.0,
    min_grid_levels=3,
    max_grid_levels=3,
)

initial_conditions = """
  Q[Shortcuts::rho] = 1.0;
  Q[Shortcuts::u]   = 0.5;
  Q[Shortcuts::v]   = 0.5;
  Q[Shortcuts::e]   = 1.0;
"""

boundary_conditions = """
  for (int i = 0; i < NumberOfUnknowns; i++) {
    Qoutside[i] = Qinside[i];
  }
"""

flux = """
  for (int i = 0; i < NumberOfUnknowns; i++) {
    F[i] = Q[i];
  }
"""

max_eigenvalue = """
  return 1.0;
"""

aderdg_solver = exahype2.dg.ader.Solver(
    name="ADERDGSolver",
    order = 4,
    unknowns={"rho": 1, "u": 1, "v": 1, "e": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.4),
)

solver = aderdg_solver

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

# tracer_particles.insert_by_coordinates(
#     [[0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]]
# )
tracer_particles.insert_with_data(
    "{0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}" #x, y, rho, fu, fv, u, v, ax, ay, mass
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