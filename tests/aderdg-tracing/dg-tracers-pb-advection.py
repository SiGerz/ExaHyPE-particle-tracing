import exahype2

project = exahype2.Project(
    namespace=["tests", "exahype2", "fv"],
    project_name="solution",
    dimensions=2,
    size=[1.0, 1.0, 1.0],
    offset=[0.0, 0.0, 0.0],
    periodic_bc=[True, True, False],
    min_end_time=10.0,
    max_end_time=10.0,
    min_grid_levels=3,
    max_grid_levels=3,
)

initial_conditions = """
  constexpr double Pi = 3.14159265358979323846;

  Q[0] = sin(2.0 * Pi * x[0]) * sin(2.0 * Pi * x[1]); // q
  Q[1] = 0.3;                                         // fluid_u
  Q[2] = 0.2;                                         // fluid_v
"""

boundary_conditions = """
  for (int i = 0; i < NumberOfUnknowns + NumberOfAuxiliaryVariables; i++) {
    Qoutside[i] = Qinside[i];
  }
"""

flux = """
  const double c = (normal == 0) ? 0.3 : 0.2;

  for (int i = 0; i < NumberOfUnknowns; i++) {
    F[i] = c * Q[i];
  }
"""

max_eigenvalue = """
  return (normal == 0) ? 0.3 : 0.2;
"""

source_term = """
  for (int i = 0; i < NumberOfUnknowns; i++) {
    S[i] = 0.0;
  }
"""


aderdg_solver = exahype2.dg.ader.Solver(
    name="ADERDGSolver",
    order = 4,
    unknowns={"q": 1, "fluid_u": 1, "fluid_v": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.4),
)

'''fv_solver = exahype2.fv.Solver(
    name="FVSolver",
    patch_size=16,
    unknowns={"q": 1, "fluid_u": 1, "fluid_v": 1},
    auxiliary_variables=0,
    time_stepping=exahype2.TimeStepping.Adaptive(0.4),
    scheme=exahype2.fv.Godunov(),
)'''

solver = aderdg_solver

solver.set_pde_terms(
    initial_condition=initial_conditions,
    boundary_conditions=boundary_conditions,
    max_eigenvalue=max_eigenvalue,
    flux=flux,
    source=source_term,
)

project.add_solver(solver)

tracer_attributes = {
    "unknowns": {
        "q": 1,
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
    "{0.2, 0.3, 0.9045084972, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0}"
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