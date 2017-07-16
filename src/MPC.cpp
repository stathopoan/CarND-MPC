#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 20;
double dt = 0.1;

size_t x_start = 0;
size_t y_start = N + x_start;
size_t psi_start = N + y_start;
size_t v_start = N + psi_start;
size_t cte_start = N + v_start;
size_t epsi_start = N + cte_start;
size_t delta_start = N + epsi_start;
size_t a_start = N-1 + delta_start;

// The reference speed
double ref_v = 50*0.44704; // MPH -> m/s

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
	// Define the cost function

	// The cost is stored is the first element of `fg`.

	fg[0] = 0;
	// The part of the cost based on the reference state.
	for (unsigned int t = 0; t < N; t++) {
		fg[0] += CppAD::pow(vars[cte_start + t], 2);
		fg[0] += CppAD::pow(vars[epsi_start + t], 2);
		fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
	}

	// Minimize the use of actuators.
	for (unsigned int t = 0; t < N - 1; t++) {
		fg[0] += 15000*CppAD::pow(vars[delta_start + t], 2);
		fg[0] += CppAD::pow(vars[a_start + t], 2);


	}

	// Minimize the value gap between sequential actuations.
	for (unsigned int t = 0; t < N - 2; t++) {
		fg[0] += 50000*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
		fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
	}

	// Add constraints

	// Initialization
	fg[1 + x_start] = vars[x_start];
	fg[1 + y_start] = vars[y_start];
	fg[1 + psi_start] = vars[psi_start];
	fg[1 + v_start] = vars[v_start];
	fg[1 + cte_start] = vars[cte_start];
	fg[1 + epsi_start] = vars[epsi_start];

	//  Rest of the constraints

	for (unsigned int t = 1; t < N; t++) {
	  // The state at time t+1 .
	  AD<double> x1 = vars[x_start + t];
	  AD<double> y1 = vars[y_start + t];
	  AD<double> psi1 = vars[psi_start + t];
	  AD<double> v1 = vars[v_start + t];
	  AD<double> cte1 = vars[cte_start + t];
	  AD<double> epsi1 = vars[epsi_start + t];

	  // The state at time t.
	  AD<double> x0 = vars[x_start + t - 1];
	  AD<double> y0 = vars[y_start + t - 1];
	  AD<double> psi0 = vars[psi_start + t - 1];
	  AD<double> v0 = vars[v_start + t - 1];
	  AD<double> cte0 = vars[cte_start + t - 1];
	  AD<double> epsi0 = vars[epsi_start + t - 1];

	  // Only consider the actuation at time t.
	  AD<double> delta0 = vars[delta_start + t - 1];
	  AD<double> a0 = vars[a_start + t - 1];

	  AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2]*x0*x0 + coeffs[3]*x0*x0*x0;
	  AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2]*x0 + 3 * coeffs[3]*x0*x0);

	  // Here's `x` to get you started.
	  // The idea here is to constraint this value to be 0.
	  //
	  // Recall the equations for the model:
	  // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
	  // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
	  // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
	  // v_[t] = v[t-1] + a[t-1] * dt
	  // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
	  // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
	  fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
	  fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
	  fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
	  fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
	  fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
	  fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
	}
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  //size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  //Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  //Set the number of constraints
  size_t n_constraints = N*6;

  // Initialize vars
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (unsigned int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set lower and upper limits for variables.

  // Set lower and upper limits for delta. The range of values δ is set to [-25, 25] in radians:
  for (unsigned int i = delta_start; i < a_start; i++) {
      vars_lowerbound[i] = -0.436332;
      vars_upperbound[i] = 0.436332;
  }
  // Set lower and upper limits for a.
  for (unsigned int i = a_start; i < n_vars; i++) {
      vars_lowerbound[i] = -1.0;
      vars_upperbound[i] = 1.0;
  }

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (unsigned int i = 0; i < delta_start; i++) {
	  vars_lowerbound[i] = -1.0e19;
	  vars_upperbound[i] = 1.0e19;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (unsigned int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set lower and upper limits for the constraints of the initial state
  constraints_lowerbound[x_start] = x;
  constraints_upperbound[x_start] = x;

  constraints_lowerbound[y_start] = y;
  constraints_upperbound[y_start] = y;

  constraints_lowerbound[psi_start] = psi;
  constraints_upperbound[psi_start] = psi;

  constraints_lowerbound[v_start] = v;
  constraints_upperbound[v_start] = v;

  constraints_lowerbound[cte_start] = cte;
  constraints_upperbound[cte_start] = cte;

  constraints_lowerbound[epsi_start] = epsi;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost is " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> resultSolution;
  cout<< "ok " << ok << endl;

  resultSolution.push_back(solution.x[delta_start]);
  resultSolution.push_back(solution.x[a_start]);

  // Add x and y values for optimal chosen trajectory
  for (unsigned int i = 0; i < N; ++i) {
	  resultSolution.push_back(solution.x[x_start + i]);
	  resultSolution.push_back(solution.x[y_start + i]);
  }

  return resultSolution;
}
