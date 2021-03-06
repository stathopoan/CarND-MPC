 CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

# Model Predictive Control with actuator delays

---

## Project description

The purpose of this project is the implementation of a Model Predictive Controller to drive a car safely around a track. The simulator provides the position of the car in the world, its speed, heading direction, steering angle and throttle. In addition, a trajectory of waypoints used as a reference  for the vehicle to follow is provided. The actuators (steering angle and throttle) have a latency of 100 ms that must be taken into account.

## The Model

The vehicle kinematic model used is the non-linear kinematic bicycle model. It ignores tire forces, gravity, and mass and is generally a simplification of a dynamic model. The model is described by the following equations:

    x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
    y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
    psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
    v_[t] = v[t-1] + a[t-1] * dt
    cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
    epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt

`x` and `y` are the coordinates of the car in 2D space, `v` its velocity, `cte` the cross track error and `epsi` the orientation error. The variable `Lf` is the distance between the center of mass of the vehicle and the front wheels.

The state of the vehicle can be represented in each timestep `dt` by a vector of these variables: `[x,y,psi,v,cte,epsi]`

The actuators in this model are two: The steering value and throttle. Using the actuators the implementation is capable of changing the vehicle state in each timestep resulting in controlling the car. 

## Timestep Length and Elapsed Duration (N & dt)

The prediction horizon `T` is a very important factor in this project. It can be derived by the formula: `T=N*dt` where `N` is the number of timesteps in the horizon and `dt` is how much time elapses between actuations. Generally `T` should be as large as possible, while `dt` should be as small as possible. `T` should be a few seconds at most, however this depends on the vehicle velocity. The greater the vehicle's velocity the larger the `T` must be in order to allow vehicle to adapt in new conditions. As initial values of `N` and `dt` i chose 25 and 0.05 s respectively resulting in `T=1.25 s`. Since my reference speed is around 50MPH~22.3 m/s a prediction of a max distance of `22.3*1.25 ~ 28m` is sufficient. Finally `dt` was set to the latency (100 ms) since i obtained the best results and `N` was set to 20. `N` shall not be too large since a long predicted horizon is useless in such a dynamic environment. With the parameters selected `T=20*0.1=2s`. This time window suffices for my reference speed giving emphasis on smooth driving (N=20) without demanding a lot of computational power (dt=0.1s). Since dt=latency the first timestep does not contribute to optimization (previous actuation values are still active).

## Polynomial Fitting and MPC Preprocessing

All computations are performed in reference to the vehicle coordinate system. The transformation functions are:

    X =   cos(-psi) * (ptsx[i] - x) - sin(-psi) * (ptsy[i] - y);
    Y =  sin(-psi) * (ptsx[i] - x) + cos(-psi) * (ptsy[i] - y);

where `ptsx` and `ptsy` are the waypoints of x and y axis respectively in global coordinates. `X` and `Y` are the coordinates of the waypoints in reference to the vehicle coordinate system. The transformation is conducted in function: `convertGlobalToLocalCoord`.

The initial state in the vehicle coordinate system is: `[0, 0, 0, v, cte, epsi]` since the car is always located at `(x,y)=(0,0)` with `psi=0`;

The velocity `v` of the vehicle is returned from the simulator in MPH scale. In order for latency equations to work correctly, a conversion to m/s is required: `v=v*0.44704; //mph -> m/s`

For fitting waypoints a 3rd order polynomial was chosen providing the degrees of freedom to adjust to the T horizon track ahead.

    // Try to fit a 3rd order polynomial to represent waypoints path and fetch coefficients
    Eigen::VectorXd coeffsFit = polyfit(carCoordptsX, carCoordptsY, 3);

## Model Predictive Control with Latency
The time between the control trigerring and its application to the vehicle is 100ms. This latency is actually introduced before the actuations are sent back to the simulator. When the latency is not taken account of, oscillations and unexpected behavior may occur. The problem originates from the fact that the optimization is conducted based on the vehicle current state and not on the state after the latency. 

My approach is to incorporate the latency model in the basic model by predicting the future car state after 100 ms time. The current state of the car is now the state after 100 ms and it is the state that the optimization will be conducted for. The state after 100ms is predicted below:

    double Lf = 2.67;
    // predict state in 100ms
    double latency = 0.1;
    px = px + v*cos(psi)*latency;
    py = py + v*sin(psi)*latency;
    psi = psi + v*delta/Lf*latency;
    v = v + acceleration*latency;

Next the transformation of waypoints from the global coordinate system to the vehicle coordinate system is required.

The NMPC trajectory is then determined by solving the control problem starting from that position with state:

`state << 0, 0, 0, v, cte, epsi;`

## Optimality

For every state value provided by the simulator an optimal trajectory for the parameters selected is computed by minimizing a cost function. The cost function is non-linear and includes the variables: cross-track error, error in heading direction, the current vehicle speed, a reference speed, the actuator values and more. The weights given in each term are chosen manually by trial and error.

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

The simulator receives back only the actuations corresponding to the first time step. At the next time step the entire optimization procedure is repeated.
## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.
* Fortran Compiler
  * Mac: `brew install gcc` (might not be required)
  * Linux: `sudo apt-get install gfortran`. Additionall you have also have to install gcc and g++, `sudo apt-get install gcc g++`. Look in [this Dockerfile](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile) for more info.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt`
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions


1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Tips

1. It's recommended to test the MPC on basic examples to see if your implementation behaves as desired. One possible example
is the vehicle starting offset of a straight line (reference). If the MPC implementation is correct, after some number of timesteps
(not too many) it should find and track the reference line.
2. The `lake_track_waypoints.csv` file has the waypoints of the lake track. You could use this to fit polynomials and points and see of how well your model tracks curve. NOTE: This file might be not completely in sync with the simulator so your solution should NOT depend on it.
3. For visualization this C++ [matplotlib wrapper](https://github.com/lava/matplotlib-cpp) could be helpful.

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!

More information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/b1ff3be0-c904-438e-aad3-2b5379f0e0c3/concepts/1a2255a0-e23c-44cf-8d41-39b8a3c8264a)
for instructions and the project rubric.

## Hints!

* You don't have to follow this directory structure, but if you do, your work
  will span all of the .cpp files here. Keep an eye out for TODOs.

## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to we ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./
