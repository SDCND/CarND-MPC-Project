#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
// size_t N = 25;
// double dt = 0.05;

size_t N = 11;
double dt = 0.15;

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

double cte_ref = 0;
double epsilon_ref = 0;
double v_ref = 80;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_ini = 0;
size_t y_ini = x_ini + N;
size_t psi_ini = y_ini + N;
size_t v_ini = psi_ini + N;
size_t cte_ini = v_ini + N;
size_t epsilon_ini = cte_ini + N;
size_t delta_ini = epsilon_ini + N;
size_t a_ini = delta_ini + N - 1;

class FG_eval 
{
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector& fg, const ADvector& vars) 
	{
    // TODO: implement MPC
    // fg a vector of constraints, x is a vector of constraints.
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    fg[0] = 0;

    for (int i = 0; i < N; i++) 
		{
      fg[0] += 1200 * CppAD::pow(vars[cte_ini + i] - cte_ref, 2);
      fg[0] += 1200 * CppAD::pow(vars[epsilon_ini + i] - epsilon_ref, 2);
      fg[0] += CppAD::pow(vars[v_ini + i] - v_ref, 2);
    }
   
    for (int i = 0; i < N - 1; i++) 
		{
      fg[0] += 20 * CppAD::pow(vars[delta_ini + i], 2);
      fg[0] += 50 * CppAD::pow(vars[a_ini + i], 2);
    }
 
    for (int i = 0; i < N - 2; i++) 
		{
      fg[0] += 250 * CppAD::pow(vars[delta_ini + i + 1] - vars[delta_ini + i], 2); 
      fg[0] += 100 * CppAD::pow(vars[a_ini + i + 1] - vars[a_ini + i], 2);
    }

    fg[x_ini + 1] = vars[x_ini];
    fg[y_ini + 1] = vars[y_ini];
    fg[psi_ini + 1] = vars[psi_ini];
    fg[v_ini + 1] = vars[v_ini];
    fg[cte_ini + 1] = vars[cte_ini];
    fg[epsilon_ini + 1] = vars[epsilon_ini];

    for (int i = 0; i < N - 1; i++) 
		{
      AD<double> x_1 = vars[x_ini + i + 1];
      AD<double> y_1 = vars[y_ini + i + 1];
      AD<double> psi_1 = vars[psi_ini + i + 1];
      AD<double> v_1 = vars[v_ini + i + 1];
      AD<double> cte_1 = vars[cte_ini + i + 1];
      AD<double> epsilon_1 = vars[epsilon_ini + i + 1];
  
      AD<double> x_0 = vars[x_ini + i];
      AD<double> y_0 = vars[y_ini + i];
      AD<double> psi_0 = vars[psi_ini + i];
      AD<double> v_0 = vars[v_ini + i];
      AD<double> cte_0 = vars[cte_ini + i];
      AD<double> epsilon_0 = vars[epsilon_ini + i];
 
      AD<double> delta_0 = vars[delta_ini + i];
      AD<double> a_0 = vars[a_ini + i];

      AD<double> f_0 = coeffs[0] + coeffs[1] * x_0 + coeffs[2] * x_0 * x_0 + coeffs[3] * x_0 * x_0 * x_0;
      AD<double> psi_d_0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x_0) + (3 * coeffs[3]* (x_0 * x_0)));

      fg[2 + x_ini + i] = x_1 - (x_0 + v_0 * CppAD::cos(psi_0) * dt);
      fg[2 + y_ini + i] = y_1 - (y_0 + v_0 * CppAD::sin(psi_0) * dt);
      fg[2 + psi_ini + i] = psi_1 - (psi_0 + v_0 * delta_0 / Lf * dt);
      fg[2 + v_ini + i] = v_1 - (v_0 + a_0 * dt);
      fg[2 + cte_ini + i] = cte_1 - ((f_0 - y_0) + (v_0 * CppAD::sin(epsilon_0) * dt));
      fg[2 + epsilon_ini + i] = epsilon_1 - ((psi_0 - psi_d_0) + v_0 * delta_0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) 
{
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsilon = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);

  for (int i = 0; i < n_vars; i++) 
	{
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

	vars[x_ini] = x;
	vars[y_ini] = y;
	vars[psi_ini] = psi;
	vars[v_ini] = v;
	vars[cte_ini] = cte;
	vars[epsilon_ini] = epsilon;
	
  for (int i = 0; i < delta_ini; i++) 
	{
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

	  for (int i = delta_ini; i < a_ini; i++) 
	{
    vars_lowerbound[i] = -0.5;
    vars_upperbound[i] = 0.5;
  }


  for (int i = a_ini; i < n_vars; i++) 
	{
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

	// TODO: Set lower and upper limits for variables.

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  for (int i = 0; i < n_constraints; i++) 
	{
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_ini] = x;
  constraints_lowerbound[y_ini] = y;
  constraints_lowerbound[psi_ini] = psi;
  constraints_lowerbound[v_ini] = v;
  constraints_lowerbound[cte_ini] = cte;
  constraints_lowerbound[epsilon_ini] = epsilon;

  constraints_upperbound[x_ini] = x;
  constraints_upperbound[y_ini] = y;
  constraints_upperbound[psi_ini] = psi;
  constraints_upperbound[v_ini] = v;
  constraints_upperbound[cte_ini] = cte;
  constraints_upperbound[epsilon_ini] = epsilon;

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
  CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound, constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {solution.x[delta_ini],solution.x[a_ini]};
}
